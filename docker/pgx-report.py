#!/usr/bin/env python3
"""
docker/pgx-report.py — Generate HTML clinical PGx reports

Produces:
  <output>/<sample>.pgx.html          — landing page (gene cards + BAM QC)
  <output>/<sample>.<gene>.pgx.html   — per-gene 17X4 detail table

Usage:
  pgx-report.py --sample <SAMPLE_ID> --output <OUTPUT_DIR> [--bam-stats <JSON>]

The script reads:
  <output>/all_genes_summary.tsv
  <output>/<GENE>/<GENE>_<SAMPLE>_detail.json  (one per gene)
  <bam_stats>   (optional; bam_stats.json written by pgx-bamstats.sh)
"""

import argparse
import glob
import json
import os
import sys
from collections import Counter
from datetime import date
from html import escape

# ── Field metadata (17 rows in the detail table) ──────────────────────────────
FIELDS = [
    ("diplotype",             "Diplotype"),
    ("haplotype1",            "Haplotype 1"),
    ("haplotype2",            "Haplotype 2"),
    ("sub_alleles",           "Sub-alleles"),
    ("alternative_diplotypes","Alternative diplotypes"),
    ("phenotype",             "Phenotype"),
    ("activity_score",        "Activity score"),
    ("sv_type",               "SV type"),
    ("copy_number",           "Copy number"),
    ("supporting_variants",   "Supporting variants"),
    ("functional_effects",    "Functional effects"),
    ("dbsnp_ids",             "dbSNP IDs"),
    ("allele_score",          "Allele score"),
    ("mean_af",               "Mean allele frequency"),
    ("phasing_method",        "Phasing method"),
    ("status",                "Status"),
    ("sv_mode",               "SV mode note"),
]

TOOLS = ["PyPGx", "Stargazer", "Aldy", "StellarPGx", "Cyrius", "PharmCAT", "OptiType", "mutserve", "VCF-Check"]

PHENOTYPE_COLORS = {
    "poor":         "#e53e3e",   # red
    "intermediate": "#ed8936",   # orange
    "normal":       "#38a169",   # green
    "ultrarapid":   "#3182ce",   # blue
    "rapid":        "#3182ce",
    "increased":    "#3182ce",
    "decreased":    "#ed8936",
    "indeterminate":"#718096",
    "uncertain":    "#718096",
    "favorable":    "#38a169",
    "unfavorable":  "#e53e3e",
}

def phenotype_color(pheno: str) -> str:
    if not pheno or pheno == "-":
        return "#718096"
    p = pheno.lower()
    for k, c in PHENOTYPE_COLORS.items():
        if k in p:
            return c
    return "#718096"


# ── CPIC / PharmVar clinical reference database ────────────────────────────────
# Static lookup for the 15 most actionable CPIC Level A/B gene-drug pairs.
# ── Gene loci (GRCh38, 1-based display coordinates) ──────────────────────────
# Used to show "chrN:start-end (GRCh38)" under the gene name on detail pages.
# HLA-A/HLA-B use the actual gene body coordinates, not the MHC extraction region.
# GSTT1 is on an alt contig.
GENE_LOCI: dict[str, str] = {
    "ABCG2":   "chr4:88,085,265-88,236,626",
    "CYP1A1":  "chr15:74,716,541-74,728,528",
    "CYP1A2":  "chr15:74,745,844-74,759,607",
    "CYP2A6":  "chr19:40,833,540-40,890,447",
    "CYP2B6":  "chr19:40,921,281-41,028,398",
    "CYP2C8":  "chr10:95,033,771-95,072,497",
    "CYP2C9":  "chr10:94,935,657-94,993,091",
    "CYP2C19": "chr10:94,759,680-94,858,547",
    "CYP2D6":  "chr22:42,116,498-42,155,810",
    "CYP2E1":  "chr10:133,517,362-133,549,123",
    "CYP3A4":  "chr7:99,753,966-99,787,184",
    "CYP3A5":  "chr7:99,645,193-99,682,996",
    "CYP4F2":  "chr19:15,863,022-15,913,074",
    "DPYD":    "chr1:97,074,742-97,924,034",
    "G6PD":    "chrX:154,528,389-154,550,018",
    "GSTM1":   "chr1:109,684,816-109,696,745",
    "GSTT1":   "chr22_KI270879v1_alt (alt contig)",
    "HLA-A":   "chr6:29,910,247-29,913,661",
    "HLA-B":   "chr6:31,321,649-31,324,666",
    "IFNL3":   "chr19:39,240,552-39,253,525",
    "NAT1":    "chr8:18,207,108-18,226,689",
    "NAT2":    "chr8:18,388,281-18,404,218",
    "NUDT15":  "chr13:48,034,725-48,050,221",
    "POR":     "chr7:75,912,154-75,989,855",
    "RYR1":    "chr19:38,430,690-38,590,564",
    "SLCO1B1": "chr12:21,128,193-21,242,796",
    "TPMT":    "chr6:18,125,310-18,158,169",
    "UGT1A1":  "chr2:233,754,269-233,779,300",
    "VKORC1":  "chr16:31,087,853-31,097,797",
    "CACNA1S": "chr1:201,006,956-201,083,927",
    "MT-RNR1": "chrM:648-1,601 (rCRS / GRCh38 chrM)",
}

# diplotype_check: lambda(diplotype_str) -> bool for genes where the diplotype
#   string must be inspected directly (VKORC1, DPYD non-reference, CYP4F2 *3,
#   RYR1 non-reference).
# Sources: cpicpgx.org, pharmvar.org (verify URLs at publication time).
# CPIC clinical knowledge base. Static data (drugs, phenotypes, dosing notes,
# CPIC tiers, URLs) lives in cpic.json so clinical staff can update it without
# editing code (F7). The dynamic per-gene diplotype predicates can't be JSON,
# so they live in DIPLOTYPE_CHECKS and are re-attached at load.
_CPIC_JSON = os.path.join(os.path.dirname(os.path.realpath(__file__)), "cpic.json")
with open(_CPIC_JSON, encoding="utf-8") as _fh:
    CPIC_DB: dict = json.load(_fh)

# diplotype_check: callable(diplotype_str) -> bool for genes whose diplotype
# string is inspected directly (VKORC1 rs9923231, DPYD non-reference, CYP4F2
# *3, ...). Verified equivalent to the pre-extraction lambdas by test_cpic.py.
DIPLOTYPE_CHECKS = {
    "DPYD":    lambda d: bool(d) and d not in ("-",) and any(v.strip() not in ("*1", "Reference", "") for v in d.split("/")),
    "VKORC1":  lambda d: any(x in (d or "") for x in ("rs9923231", "*S1", "*H1", "*H2")),
    "RYR1":    lambda d: bool(d) and d not in (None, "-", "Reference/Reference", "*1/*1", "*ference/*ference", ""),
    "CACNA1S": lambda d: bool(d) and d not in (None, "-", "Reference/Reference", "*1/*1", ""),
    "MT-RNR1": lambda d: bool(d) and d not in (None, "-", "Reference", ""),
    "CYP4F2":  lambda d: "*3" in (d or ""),
    "ABCG2":   lambda d: "rs2231142" in (d or ""),
    "HLA-A":   lambda d: "A*31:01" in (d or ""),
    "HLA-B":   lambda d: any(x in (d or "") for x in ("B*57:01", "B*15:02", "B*58:01", "B*13:01", "B*57:03")),
    "CYP3A4":  lambda d: "*22" in (d or ""),
    "IFNL3":   lambda d: "Reference/Reference" not in (d or "") and bool(d) and d not in ("-", ""),
}
for _g, _entry in CPIC_DB.items():
    _entry["diplotype_check"] = DIPLOTYPE_CHECKS.get(_g)


def _pheno_cat(phenotype: str) -> str:
    """Normalise a phenotype string to a lowercase keyword category."""
    p = (phenotype or "").lower().replace("-", " ").replace("_", " ").strip()
    if not p or p == "-":
        return "unknown"
    # Check ultrarapid before rapid (substring order matters)
    if "ultrarapid" in p.replace(" ", "") or "ultra rapid" in p:
        return "ultrarapid"
    if "poor" in p:       return "poor"
    if "rapid" in p:      return "rapid"      # rapid but not ultrarapid (caught above)
    if "intermediate" in p: return "intermediate"
    if "slow" in p:       return "slow"
    if "deficient" in p:  return "deficient"
    if "hemizygous" in p: return "hemizygous"
    if "homozygous" in p: return "homozygous"
    if "heterozygous" in p: return "heterozygous"
    if "variable" in p:   return "variable"
    if "decreased" in p or "reduced" in p: return "decreased"
    if "increased" in p:  return "increased"
    if "susceptible" in p or "malignant" in p: return "susceptible"
    if "non expresser" in p or "nonexpresser" in p: return "poor"  # CYP3A5 *3/*3
    if "normal" in p or "reference" in p or "wild" in p: return "normal"
    return "other"


def _get_tier(gene: str, pheno_cat: str, diplotype: str,
              all_diplotypes: list[str] | None = None) -> str | None:
    """Return 'high', 'moderate', 'informational', or None (not flagged).

    For diplotype_check genes, all_diplotypes (one per tool that called the gene)
    are checked individually so that nomenclature differences between tools
    (e.g. rs9923231 vs *S1 vs *H1 for VKORC1) do not block detection.
    The finding is raised only if ≥2 tool diplotypes trigger the lambda.
    """
    entry = CPIC_DB.get(gene)
    if not entry:
        return None
    dc = entry.get("diplotype_check")
    if dc:
        diplos = all_diplotypes if all_diplotypes else ([diplotype] if diplotype else [])
        min_t = entry.get("min_tools", 2)
        n_triggered = sum(1 for d in diplos if d and d not in ("-", "") and dc(d))
        if n_triggered >= min_t:
            if gene in ("DPYD", "RYR1", "TPMT", "NUDT15", "HLA-A", "HLA-B"):
                return "high"
            elif gene in ("VKORC1", "ABCG2"):
                return "moderate"
            else:
                return "informational"
    # Phenotype-based tier
    if pheno_cat in ("normal", "unknown", "other"):
        return None
    if pheno_cat in entry.get("high_pheno", []):
        return "high"
    if pheno_cat in entry.get("moderate_pheno", []):
        return "moderate"
    return None


def _get_landing_note(gene: str, pheno_cat: str, diplotype: str,
                      all_diplotypes: list[str] | None = None) -> str:
    """Return the brief TL;DR clinical note for a finding."""
    entry = CPIC_DB.get(gene, {})
    notes = entry.get("landing_notes", {})
    dc = entry.get("diplotype_check")
    diplos = all_diplotypes if all_diplotypes else ([diplotype] if diplotype else [])

    # HLA genes: generate note from whichever risk alleles are present
    risk_alleles = entry.get("risk_alleles")
    if risk_alleles and dc:
        found = [note for allele, note in risk_alleles.items()
                 if any(allele in (d or "") for d in diplos if d not in ("-", ""))]
        if found:
            return " | ".join(found)
        return entry.get("no_risk_note", "")

    if dc and any(dc(d) for d in diplos if d and d not in ("-", "")) and "_diplotype" in notes:
        return notes["_diplotype"]
    return notes.get(pheno_cat, notes.get("_diplotype", ""))


def build_clinical_findings_section(genes_data: list, sample: str = "",
                                     genes_rel_prefix: str = "",
                                     embedded: bool = False) -> str:
    """Build the Key Clinical Findings TL;DR section for the landing page.

    embedded — True when gene detail panels are inlined in the same HTML file;
               links use onclick='pgxShowGene()' instead of separate file hrefs.
    """
    findings: dict[str, list] = {"high": [], "moderate": [], "informational": []}

    for gd in genes_data:
        gene             = gd["gene"]
        pheno            = gd["consensus_phenotype"]
        diplo            = gd["consensus_diplotype"]
        n_agree          = gd["n_agree"]
        n_called         = gd["n_called"]
        all_tool_diplos  = gd.get("all_tool_diplotypes", [])

        entry = CPIC_DB.get(gene, {})
        dc    = entry.get("diplotype_check")
        pheno_cat = _pheno_cat(pheno)

        if dc:
            # Diplotype-check gene: count how many tools trigger the lambda.
            # Different tools use different nomenclature for the same variant
            # (e.g. rs9923231 / *S1 / *H1 for VKORC1), so we do NOT use the
            # consensus diplotype n_agree filter here — we count lambda hits.
            min_t = entry.get("min_tools", 2)
            n_triggered = sum(1 for d in all_tool_diplos if d and d not in ("-", "") and dc(d))
            if n_triggered < min_t:
                continue
            tier = _get_tier(gene, pheno_cat, diplo, all_tool_diplos)
            note = _get_landing_note(gene, pheno_cat, diplo, all_tool_diplos)
            concordance_warn = n_triggered < n_called
        else:
            # Phenotype-based gene: require ≥2 tools to agree on diplotype.
            if n_agree < 2:
                continue
            tier = _get_tier(gene, pheno_cat, diplo)
            note = _get_landing_note(gene, pheno_cat, diplo)
            concordance_warn = n_agree < n_called

        if tier is None:
            continue
        findings[tier].append({
            "gene":             gene,
            "diplotype":        diplo,
            "phenotype":        pheno,
            "note":             note,
            "concordance_warn": concordance_warn,
            "n_agree":          n_agree,
            "n_called":         n_called,
        })

    all_findings = findings["high"] + findings["moderate"] + findings["informational"]
    if not all_findings:
        return """
    <div class="section">
        <h2>Key Clinical Findings</h2>
        <div class="cf-normal">
            <span style="font-size:1.4rem">&#10003;</span>
            <div><strong>No clinically actionable pharmacogenomic variants identified.</strong><br>
            All tested genes show normal/reference function based on concordant tool results.</div>
        </div>
    </div>"""

    def render_tier(tier_findings, tier_css, tier_icon, tier_label):
        if not tier_findings:
            return ""
        html = f'<div class="cf-tier-label {tier_css}-label">{tier_icon}&nbsp; {tier_label}</div>'
        for f in tier_findings:
            pill_color = phenotype_color(f["phenotype"])
            pheno_str  = f["phenotype"] if f["phenotype"] not in ("-", "") else "—"
            diplo_str  = f["diplotype"] if f["diplotype"] not in ("-", "") else "—"
            warn_html  = (
                f'<span class="cf-concordance-warn" title="{f["n_agree"]}/{f["n_called"]} tools agree">'
                f'&#9888; Discordant ({f["n_agree"]}/{f["n_called"]})</span>'
                if f["concordance_warn"] else ""
            )
            if embedded:
                detail_href  = "#"
                detail_click = f" onclick=\"pgxShowGene('{f['gene']}'); return false;\""
            elif sample and genes_rel_prefix:
                detail_href  = f"{genes_rel_prefix}/{f['gene']}/{sample}.{f['gene']}.pgx.html"
                detail_click = ""
            elif sample:
                detail_href  = f"{sample}.{f['gene']}.pgx.html"
                detail_click = ""
            else:
                detail_href  = "#"
                detail_click = ""
            html += f"""
            <a href="{detail_href}"{detail_click} class="cf-finding-link">
            <div class="cf-finding {tier_css}">
                <div class="cf-finding-header">
                    <span class="cf-gene">{f["gene"]}</span>
                    <span class="cf-diplo">{diplo_str}</span>
                    <span class="gene-pheno-pill" style="background:{pill_color};font-size:0.72rem;padding:0.15rem 0.5rem">{pheno_str}</span>
                    {warn_html}
                    <span class="cf-detail-arrow">&#8594; Detail</span>
                </div>
                <div class="cf-note">{f["note"]}</div>
            </div>
            </a>"""
        return html

    body  = render_tier(findings["high"],         "cf-high",     "&#9888;", "High Priority")
    body += render_tier(findings["moderate"],      "cf-moderate", "&#9679;", "Moderate")
    body += render_tier(findings["informational"], "cf-info",     "&#8505;", "Informational")

    return f"""
    <div class="section">
        <h2>Key Clinical Findings</h2>
        <p class="cf-disclaimer">Based on concordant results (&#8805;2/4 tools). For research and clinical decision support only — not a standalone diagnostic.
        Consult CPIC guidelines and a clinical pharmacist before prescribing.</p>
        <div class="cf-container">
{body}
        </div>
    </div>"""


def build_gene_cpic_section(gene: str, phenotype: str, diplotype: str,
                             all_tool_diplos: list[str] | None = None,
                             id_prefix: str = "") -> str:
    """Build the comprehensive CPIC reference section for a gene detail page."""
    entry = CPIC_DB.get(gene)
    if not entry:
        return ""

    pheno_cat = _pheno_cat(phenotype)
    tier      = _get_tier(gene, pheno_cat, diplotype, all_tool_diplos)
    note      = _get_landing_note(gene, pheno_cat, diplotype, all_tool_diplos)

    tier_badge = {
        "high":          '<span class="badge badge-red" style="font-size:0.85rem;padding:0.25rem 0.75rem">High priority</span>',
        "moderate":      '<span class="badge badge-amber" style="font-size:0.85rem;padding:0.25rem 0.75rem">Moderate</span>',
        "informational": '<span class="badge badge-blue" style="font-size:0.85rem;padding:0.25rem 0.75rem">Informational</span>',
    }.get(tier or "", '<span class="badge badge-grey" style="font-size:0.85rem;padding:0.25rem 0.75rem">Normal / Not flagged</span>')

    # M2: escalate the note's visual weight by tier so a high-consequence warning
    # (e.g. "contraindicated") no longer looks identical to an informational note.
    _note_mod = {"high": " cpic-patient-note--high",
                 "moderate": " cpic-patient-note--moderate"}.get(tier or "", "")
    patient_note_html = (
        f'<div class="cpic-patient-note{_note_mod}">&#128203; {note}</div>' if note else ""
    )

    drug_rows = ""
    for d in entry.get("drugs", []):
        lvl_cls = f'cpic-level-{d["level"].lower()}'
        drug_rows += f"""
            <tr>
                <td><a href="{d['url']}" target="_blank" class="cpic-link">{d['name']}</a></td>
                <td><span class="cpic-level-badge {lvl_cls}">{d['level']}</span></td>
                <td>{d['rec']}</td>
            </tr>"""

    pharmvar_url = entry.get("pharmvar_url")
    pharmvar_html = (
        f'<a href="{pharmvar_url}" target="_blank" class="cpic-link">PharmVar — {gene} allele nomenclature &#8599;</a>'
        if pharmvar_url else ""
    )
    sep = "&ensp;|&ensp;" if pharmvar_html else ""

    return f"""
    <div class="cpic-section" id="{id_prefix}cpic-reference">
        <div class="cpic-header">
            <div>
                <h2>CPIC Clinical Reference — {gene}</h2>
                <p class="cpic-desc">{entry.get("desc", "")}</p>
            </div>
            <div style="text-align:right;flex-shrink:0">
                <div class="qc-label">This patient</div>
                <div style="margin-top:0.35rem">{tier_badge}</div>
            </div>
        </div>
        {patient_note_html}
        <div class="cpic-drug-wrap table-scroll">
            <table class="cpic-drug-table">
                <thead><tr>
                    <th>Drug / Drug class</th>
                    <th>CPIC Level</th>
                    <th>Recommendation</th>
                </tr></thead>
                <tbody>{drug_rows}</tbody>
            </table>
        </div>
        <div class="cpic-footer">
            <a href="https://cpicpgx.org" target="_blank" class="cpic-link">cpicpgx.org &#8599;</a>
            {sep}{pharmvar_html}
        </div>
    </div>"""


def concordance_color(n_agree: int, n_called: int) -> tuple[str, str]:
    """Return (card_class, badge_text) based on tool agreement."""
    if n_called == 0:
        return ("card-no-data",  "No data")
    frac = n_agree / n_called
    if frac == 1.0:
        return ("card-green",  f"{n_agree}/{n_called}")
    elif frac >= 0.75:
        return ("card-amber",  f"{n_agree}/{n_called}")
    elif frac >= 0.5:
        return ("card-orange", f"{n_agree}/{n_called}")
    else:
        return ("card-red",    f"{n_agree}/{n_called}")


def normalize_diplotype(d: str) -> str:
    """Sort alleles alphabetically so *2/*4 == *4/*2."""
    if not d or d == "-":
        return "-"
    sep = "/" if "/" in d else None
    if sep is None:
        return d
    parts = [p.strip() for p in d.split(sep)]
    return "/".join(sorted(parts))


def compute_concordance(tool_data: dict) -> tuple[str, str, str, int]:
    """
    Returns (consensus_diplotype, consensus_phenotype, card_class, n_agree).
    """
    diplotypes = []
    phenotypes = []
    for tool in TOOLS:
        td = tool_data.get(tool, {})
        d = normalize_diplotype(td.get("diplotype", "-"))
        p = td.get("phenotype", "-")
        if d and d != "-":
            diplotypes.append(d)
        if p and p != "-":
            phenotypes.append(p)

    n_called = len(diplotypes)
    if n_called == 0:
        return ("-", "-", "card-no-data", 0)

    most_common_dip, n_agree = Counter(diplotypes).most_common(1)[0]
    most_common_pheno = Counter(phenotypes).most_common(1)[0][0] if phenotypes else "-"

    card_class, _ = concordance_color(n_agree, n_called)
    return (most_common_dip, most_common_pheno, card_class, n_agree)


def verdict_card(verdict: dict) -> tuple[str, str, str, int, int]:
    """Map the authoritative detail.json verdict to
    (diplotype, phenotype, card_class, n_agree, n_called).

    This is the consumer of pgx-compare.py's single-source verdict. NO_CALL and
    DISCORDANT deliberately assert NO phenotype, so they cannot drive a clinical
    recommendation downstream. Replaces the old per-page concordance recompute
    that ignored the coverage gate and resolved ties as a silent winner.
    """
    st = verdict.get("status", "no_data")
    n_agree = int(verdict.get("n_agree", 0) or 0)
    n_called = int(verdict.get("n_called", 0) or 0)
    if st == "no_call":
        return ("NO_CALL", "-", "card-no-data", n_agree, n_called)
    if st == "discordant":
        return ("DISCORDANT", "-", "card-red", n_agree, n_called)
    if st == "no_data":
        return ("-", "-", "card-no-data", 0, 0)
    # Colour by the VERDICT, not the raw tool fraction: an authority/phenotype
    # resolved call (e.g. PharmCAT at 1/4 tools) is a confident green, not red.
    # The raw "n/N tools" detail is still shown in the per-card concordance badge.
    card_class = "card-green" if st == "concordant" else "card-amber"
    return (verdict.get("consensus_diplotype", "-"),
            verdict.get("consensus_phenotype", "-"), card_class, n_agree, n_called)


# ── CSS / JS shared across all pages ──────────────────────────────────────────
SHARED_CSS = """
:root {
    --primary: #1a365d;
    --primary-light: #2b6cb0;
    --bg: #f3f6fa;
    --card-bg: #ffffff;
    --border: #dde3ec;
    --text: #1e2a3a;
    --muted: #64748b;
    --green: #14532d;   /* darkened to clear WCAG AAA (>=7:1) on --green-bg */
    --green-bg: #dcfce7;
    --amber: #78350f;   /* AAA on --amber-bg */
    --amber-bg: #fef3c7;
    --orange: #7c2d12;  /* AAA on --orange-bg */
    --orange-bg: #ffedd5;
    --red: #7f1d1d;     /* AAA on --red-bg */
    --red-bg: #fee2e2;
    --blue: #1e3a5f;
    --blue-bg: #dbeafe;
    --grey: #334155;    /* AAA on --grey-bg */
    --grey-bg: #f1f5f9;
}
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: 'Inter', 'Segoe UI', system-ui, -apple-system, sans-serif;
    background: var(--bg);
    color: var(--text);
    line-height: 1.6;
    font-size: 15px;
}
header {
    background: var(--primary);
    color: white;
    padding: 0.9rem 2rem;
    display: flex;
    align-items: center;
    gap: 1rem;
    border-bottom: 3px solid #2b6cb0;
}
header .logo {
    font-size: 1.25rem;
    font-weight: 700;
    letter-spacing: -0.3px;
}
header .subtitle {
    font-size: 0.78rem;
    opacity: 0.8;
    margin-top: 0.1rem;
}
header .spacer { flex: 1; }
header .report-meta {
    text-align: right;
    font-size: 0.75rem;
    opacity: 0.9;
    line-height: 1.75;
}
.container { max-width: 1200px; margin: 0 auto; padding: 2rem 1.5rem; }
h2 {
    font-family: 'Source Serif 4', 'Georgia', serif;
    font-size: 1.1rem;
    font-weight: 700;
    color: var(--primary);
    margin-bottom: 0.85rem;
    padding-bottom: 0.35rem;
    border-bottom: 1px solid var(--border);
}
h3 { font-size: 0.95rem; font-weight: 600; color: var(--primary); margin-bottom: 0.75rem; }
.section {
    margin-bottom: 2.25rem;
    padding-top: 1.75rem;
    border-top: 1px solid var(--border);
}
.section:first-child { border-top: none; padding-top: 0; }
.badge {
    display: inline-block;
    padding: 0.15rem 0.55rem;
    border-radius: 4px;
    font-size: 0.72rem;
    font-weight: 600;
    letter-spacing: 0.02em;
}
.badge-green  { background: var(--green-bg);  color: var(--green);  }
.badge-amber  { background: var(--amber-bg);  color: var(--amber);  }
.badge-orange { background: var(--orange-bg); color: var(--orange); }
.badge-red    { background: var(--red-bg);    color: var(--red);    }
.badge-grey   { background: var(--grey-bg);   color: var(--grey);   }
.badge-blue   { background: var(--blue-bg);   color: var(--blue);   }
/* At-a-glance headline summary strip (clinician 1-second scan) */
.summary-strip { display: flex; flex-wrap: wrap; gap: 0.45rem 0.7rem; align-items: center;
    margin-bottom: 1.75rem; padding: 0.6rem 0.9rem; background: var(--card-bg);
    border: 1px solid var(--border); border-radius: 8px; }
.summary-total { font-weight: 700; font-size: 0.95rem; color: var(--primary); }
.summary-flagged { font-weight: 700; font-size: 0.8rem; color: var(--red); }
.summary-flagged.summary-ok { color: var(--green); }
.summary-pill { font-size: 0.76rem; font-weight: 600; padding: 0.15rem 0.55rem; border-radius: 999px; }
.summary-pill b { font-weight: 800; }
.summary-pill.card-red     { background: var(--red-bg);    color: var(--red); }
.summary-pill.card-no-data { background: var(--grey-bg);   color: var(--grey); }
.summary-pill.card-orange  { background: var(--orange-bg); color: var(--orange); }
.summary-pill.card-amber   { background: var(--amber-bg);  color: var(--amber); }
.summary-pill.card-green   { background: var(--green-bg);  color: var(--green); }
table { width: 100%; border-collapse: collapse; font-size: 0.84rem; }
th {
    background: var(--primary);
    color: white;
    padding: 0.55rem 0.75rem;
    text-align: left;
    font-weight: 600;
    font-size: 0.73rem;
    text-transform: uppercase;
    letter-spacing: 0.06em;
}
thead tr th { border-bottom: 2px solid #2b6cb0; }
td { padding: 0.5rem 0.75rem; border-bottom: 1px solid var(--border); vertical-align: top; font-variant-numeric: tabular-nums; }
tr:nth-child(even) td { background: #f8fafc; }
tr:hover td { background: #eef2f7; }
tr:last-child td { border-bottom: none; }
.dash { color: var(--muted); }
footer {
    text-align: center;
    padding: 1.5rem 2rem;
    font-size: 0.73rem;
    color: var(--muted);
    border-top: 2px solid var(--border);
    margin-top: 3rem;
    line-height: 1.9;
}
.back-link {
    display: inline-flex;
    align-items: center;
    gap: 0.4rem;
    color: var(--primary-light);
    text-decoration: none;
    font-size: 0.875rem;
    font-weight: 500;
    margin-bottom: 1.5rem;
}
.back-link:hover { text-decoration: underline; }
@media print {
    header .report-meta { opacity: 1; }
    .gene-panel[hidden] { display: block !important; page-break-before: always; }
    .gene-panel { display: block !important; }
    nav, .gene-page-nav, .back-link { display: none; }
    body { font-size: 10.5pt; background: white; color: #000; }
    * { box-shadow: none !important; }
    a[href]::after { content: none; }
    .gene-grid { grid-template-columns: repeat(3, 1fr) !important; gap: 0.5rem !important; }
    .qc-bar { flex-wrap: wrap; }
    .container { padding: 0.5rem; }
    th { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
    .gene-card { border: 1px solid #ccc; break-inside: avoid; }
    .report-status-banner { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
    /* C1: preserve severity colour in print/PDF. Without this the archived report
       drops all status backgrounds (a discordant/NO_CALL prints like a normal
       call, and white-text pheno pills become invisible on white). */
    .badge, .gene-card, .cpic-level-badge, .gene-pheno-pill, .cpic-patient-note,
    .status-final, .status-prelim {
        -webkit-print-color-adjust: exact; print-color-adjust: exact;
    }
}

/* M1: responsive reflow — previously the ONLY @media block was print. */
.table-scroll { overflow-x: auto; -webkit-overflow-scrolling: touch; }
@media (max-width: 768px) {
    .container { padding: 1rem 0.75rem; }
    header { padding: 0.75rem 1rem; flex-wrap: wrap; }
    .qc-bar { flex-wrap: wrap; }
    .gene-grid { grid-template-columns: repeat(auto-fill, minmax(150px, 1fr)); }
    .table-scroll table, .table-scroll .cpic-drug-table { min-width: 560px; }
}
/* ── CAP/CLIA compliance elements ─────────────────────────────────── */
.report-status-banner {
    text-align: center;
    font-size: 0.72rem;
    font-weight: 700;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    padding: 0.3rem 1rem;
}
.status-final    { background: #166534; color: #fff; }
.status-prelim   { background: #92400e; color: #fff; }
.status-research { background: #475569; color: #fff; }
.lab-block {
    text-align: right;
    border-right: 1px solid rgba(255,255,255,0.25);
    margin-right: 1rem;
    padding-right: 1rem;
}
.lab-block .lab-name   { font-size: 0.82rem; font-weight: 600; color: #fff; }
.lab-block .lab-detail { font-size: 0.72rem; opacity: 0.85; color: #fff; line-height: 1.6; }
.accession-row { font-size: 0.78rem; color: var(--muted); margin-top: 0.15rem; }
.accession-row .accession-id { font-family: monospace; font-weight: 600; color: var(--primary); }
.auth-block {
    font-size: 0.75rem;
    color: var(--muted);
    margin-bottom: 0.5rem;
    padding-bottom: 0.5rem;
    border-bottom: 1px solid var(--border);
}
.auth-block .auth-label { font-weight: 600; }
.auth-block .auth-name  { color: var(--primary); font-weight: 500; }
"""

# ── Landing page ──────────────────────────────────────────────────────────────
LANDING_EXTRA_CSS = """
.sample-banner {
    background: white;
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1.25rem 1.5rem;
    margin-bottom: 2rem;
    display: flex;
    align-items: flex-start;
    gap: 2rem;
    flex-wrap: wrap;
}
.sample-id { font-size: 1.5rem; font-weight: 700; color: var(--primary); }
.sample-meta { font-size: 0.85rem; color: var(--muted); margin-top: 0.2rem; }
.qc-bar {
    background: white;
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 0.7rem 1.25rem;
    display: flex;
    align-items: stretch;
    flex-wrap: wrap;
    row-gap: 0.5rem;
    margin-bottom: 1.75rem;
}
.qc-bar-item {
    display: flex;
    flex-direction: column;
    align-items: flex-start;
    padding: 0.2rem 1rem;
    min-width: 88px;
}
.qc-bar-item:first-child { padding-left: 0.25rem; }
.qc-bar-sep {
    width: 1px;
    background: var(--border);
    align-self: stretch;
    margin: 0.2rem 0;
    flex-shrink: 0;
}
.qc-label { font-size: 0.72rem; text-transform: uppercase; letter-spacing: 0.07em; color: var(--muted); font-weight: 500; }
.qc-value { font-size: 1.05rem; font-weight: 700; color: var(--primary); margin-top: 0.1rem; font-variant-numeric: tabular-nums; }
.qc-unit  { font-size: 0.7rem; color: var(--muted); font-weight: 400; }
.gene-grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
    gap: 1rem;
}
.gene-card {
    background: white;
    border-radius: 8px;
    padding: 1rem 1.1rem;
    cursor: pointer;
    text-decoration: none;
    color: var(--text);
    transition: transform 0.12s, box-shadow 0.12s;
    border-left: 5px solid #ccc;
    display: block;
}
.gene-card:hover { transform: translateY(-2px); box-shadow: 0 4px 16px rgba(0,0,0,0.12); }
.card-green  { border-left-color: #38a169; }
.card-amber  { border-left-color: #d69e2e; }
.card-orange { border-left-color: #ed8936; }
.card-red    { border-left-color: #e53e3e; }
.card-no-data{ border-left-color: #a0aec0; }
.gene-name   { font-size: 1rem; font-weight: 700; }
.gene-status { font-size: 0.72rem; font-weight: 700; margin-top: 0.15rem; color: var(--text); }
.grid-divider { grid-column: 1 / -1; font-size: 0.74rem; font-weight: 700; text-transform: uppercase; letter-spacing: 0.06em; padding: 0.6rem 0 0.15rem; border-bottom: 1px solid var(--border); }
.grid-divider-flagged { color: var(--red); }
.grid-divider-normal  { color: var(--green); margin-top: 0.5rem; }
.gene-diplo  { font-size: 0.8rem; color: var(--muted); margin-top: 0.25rem; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.gene-pheno-pill {
    display: inline-block;
    margin-top: 0.5rem;
    padding: 0.2rem 0.55rem;
    border-radius: 999px;
    font-size: 0.7rem;
    font-weight: 600;
    color: white;
}
.gene-concord { font-size: 0.7rem; margin-top: 0.4rem; }
.gene-authority { margin-top: 0.3rem; }
.gene-authority .badge { font-size: 0.62rem; padding: 0.1rem 0.4rem; }
.gene-xcheck { margin-top: 0.25rem; }
.gene-xcheck .badge { font-size: 0.6rem; padding: 0.1rem 0.4rem; opacity: 0.92; }
.xcheck-note { margin: 0.75rem 0 0; padding: 0.5rem 0.75rem; background: #f1f5f9; border-left: 4px solid var(--primary-light); border-radius: 4px; font-size: 0.8rem; color: var(--muted); }
.gene-depth   { font-size: 0.7rem; margin-top: 0.25rem; color: var(--muted); }
/* ── Sequencing depth coverage flags ── */
.depth-ok       { color: #2d7d46; font-weight: 700; }   /* ≥80% at ≥30× */
.depth-caution  { color: #b38600; font-weight: 700; }   /* 50–79% */
.depth-poor     { color: #d46b08; font-weight: 700; }   /* 10–49% */
.depth-critical { color: #c0392b; font-weight: 700; }   /* <10%  */
/* ── Key Clinical Findings section ── */
.cf-disclaimer { font-size: 0.75rem; color: var(--muted); margin-bottom: 1rem; font-style: italic; }
.cf-container { display: flex; flex-direction: column; gap: 0.4rem; }
.cf-tier-label { font-size: 0.72rem; font-weight: 700; text-transform: uppercase; letter-spacing: 0.08em; padding: 0.6rem 0 0.2rem; }
.cf-high-label     { color: #742a2a; }
.cf-moderate-label { color: #744210; }
.cf-info-label     { color: #2a4365; }
.cf-finding-link { display: block; text-decoration: none; color: inherit; border-radius: 8px; transition: transform 0.12s ease, box-shadow 0.12s ease; }
.cf-finding-link:hover { transform: translateY(-2px); box-shadow: 0 4px 12px rgba(0,0,0,0.12); }
.cf-finding-link:hover .cf-detail-arrow { opacity: 1; }
.cf-detail-arrow { font-size: 0.72rem; color: var(--primary-light); opacity: 0; transition: opacity 0.12s ease; margin-left: auto; white-space: nowrap; }
.cf-finding { background: white; border-radius: 8px; padding: 0.85rem 1.1rem; border-left: 5px solid #ccc; }
.cf-high     { border-left-color: #e53e3e; }
.cf-moderate { border-left-color: #d69e2e; }
.cf-info     { border-left-color: #3182ce; }
.cf-finding-header { display: flex; align-items: center; gap: 0.65rem; flex-wrap: wrap; margin-bottom: 0.35rem; }
.cf-gene  { font-weight: 700; font-size: 0.95rem; color: var(--primary); min-width: 70px; }
.cf-diplo { font-size: 0.8rem; color: var(--muted); font-family: monospace; }
.cf-note  { font-size: 0.82rem; color: var(--text); line-height: 1.5; }
.cf-concordance-warn { font-size: 0.72rem; color: #744210; background: #fefcbf; border-radius: 4px; padding: 0.1rem 0.4rem; }
.cf-normal { background: #c6f6d5; color: #276749; border-radius: 8px; padding: 1rem 1.25rem; font-size: 0.88rem; display: flex; align-items: center; gap: 0.75rem; }
.legend-row {
    display: flex;
    gap: 1.2rem;
    flex-wrap: wrap;
    font-size: 0.78rem;
    margin-bottom: 1.25rem;
    align-items: center;
}
.depth-details > summary { cursor: pointer; font-size: 1.15rem; font-weight: 600; color: var(--primary); list-style: none; }
.depth-details > summary::-webkit-details-marker { display: none; }
.depth-details > summary::before { content: "▸ "; color: var(--muted); }
.depth-details[open] > summary::before { content: "▾ "; }
.depth-hint { font-size: 0.78rem; font-weight: 400; color: var(--muted); }
.legend-authority { font-size: 0.78rem; margin: -0.5rem 0 1.25rem; }
.legend-auth-title { color: var(--muted); font-weight: 600; margin-bottom: 0.4rem; }
.legend-auth-row { display: flex; align-items: baseline; gap: 0.5rem; margin: 0.25rem 0; }
.legend-auth-row .badge { flex: 0 0 auto; }
.legend-auth-row > span:last-child { color: var(--muted); }
.legend-item { display: flex; align-items: center; gap: 0.4rem; }
.legend-swatch { width: 14px; height: 14px; border-radius: 3px; }
"""

def fmt_num(n, suffix=""):
    """Format large numbers with comma separators."""
    try:
        return f"{int(n):,}{suffix}"
    except Exception:
        return str(n)

def bam_stats_cards(bs: dict) -> str:
    items = [
        ("Total Reads",  fmt_num(bs.get("total_reads", 0)),      ""),
        ("Mapped",       f"{bs.get('mapped_pct','—')}",           "%"),
        ("Duplicates",   f"{bs.get('duplicate_pct','—')}",        "%"),
        ("Depth",        f"{bs.get('mean_depth_genome','—')}",    "×"),
        ("Read Length",  f"{bs.get('read_length','—')}",          "\u202fbp"),
        ("Insert Size",  f"{bs.get('insert_size_mean','—')}",     "\u202fbp"),
        ("MAPQ\u226520", f"{bs.get('mapq20_pct','—')}",           "%"),
        ("Sex",          f"{bs.get('inferred_sex','—')}",         ""),
        ("X/Y Ratio",    f"{bs.get('xy_depth_ratio','—')}",       ""),
        ("Error Rate",   f"{bs.get('error_rate','—')}",           ""),
    ]
    parts = []
    for i, (label, val, unit) in enumerate(items):
        if i > 0:
            parts.append('<div class="qc-bar-sep"></div>')
        parts.append(
            f'<div class="qc-bar-item">'
            f'<span class="qc-label">{label}</span>'
            f'<span class="qc-value">{val}<span class="qc-unit">{unit}</span></span>'
            f'</div>'
        )
    return f'<div class="qc-bar">{"".join(parts)}</div>'


# Genes where low mean gene depth (&lt;30×) materially reduces calling sensitivity.
# Validation study (TTSH 2026-03): DPYD rare variants and RYR1 het variants were
# missed in samples with gene depth &lt;30×; all tools performed adequately above that.
_DEPTH_SENSITIVE_GENES: dict[str, int] = {
    "DPYD": 30,   # rare het variants (*5, *9A) missed below 30× in TTSH cohort
    "RYR1": 30,   # het MH-susceptibility variants missed below 30× in TTSH cohort
}


def _depth_css(pct_ge_30x) -> str:
    """Return CSS class for a ≥30× coverage fraction value.

    Thresholds align with CAP/AMP PGx guidance and tool recommendations:
      ≥80%  → depth-ok       (green)  — reliable calling
      50–79% → depth-caution (amber)  — SNV OK, SV uncertain
      10–49% → depth-poor    (orange) — results may be unreliable
       <10%  → depth-critical (red)   — insufficient for calling
    """
    try:
        v = float(pct_ge_30x)
    except (TypeError, ValueError):
        return ""
    if v >= 80:
        return "depth-ok"
    if v >= 50:
        return "depth-caution"
    if v >= 10:
        return "depth-poor"
    return "depth-critical"


def gene_depth_table(bs: dict) -> str:
    gd = bs.get("gene_depth", {})
    if not gd:
        return ""
    rows = ""
    for gene in sorted(gd.keys()):
        d = gd[gene]
        if d.get("note") == "alt_contig":
            rows += f"""
        <tr>
            <td>{gene}</td>
            <td colspan="3" style="color:var(--muted);font-style:italic">Alt contig (chr22_KI270879v1_alt) — not measurable on primary assembly</td>
        </tr>"""
        else:
            mean  = d.get('mean', '—')
            p20   = d.get('pct_ge_20x', '—')
            p30   = d.get('pct_ge_30x', '—')
            p30_cls = _depth_css(p30)
            p30_cell = f'<span class="{p30_cls}">{p30}%</span>' if p30_cls else f'{p30}%'
            rows += f"""
        <tr>
            <td>{gene}</td>
            <td>{mean}X</td>
            <td>{p20}%</td>
            <td>{p30_cell}</td>
        </tr>"""
    return f"""
    <details class="section depth-details">
        <summary>Per-gene sequencing depth <span class="depth-hint">(click to expand)</span></summary>
        <div class="table-scroll" style="margin-top:0.75rem">
        <table>
            <thead><tr>
                <th>Gene</th><th>Mean depth</th><th>≥20X</th><th>≥30X</th>
            </tr></thead>
            <tbody>{rows}</tbody>
        </table>
        </div>
    </details>"""


def build_landing(sample: str, bam: str, genes_data: list, bs: dict | None, out_dir: str,
                  genes_rel_prefix: str = "",
                  gene_fragments: dict | None = None,
                  lab_info: dict | None = None,
                  provenance: dict | None = None,
                  axiom: list | None = None):
    """Build <sample>_pgx_report.html landing page with embedded gene detail panels."""

    gene_cards_html = ""
    gene_depth_map = (bs or {}).get("gene_depth", {})
    # M4: float criticals (discordant / no-call) to the top-left so they are scanned
    # first. C2: a non-colour status token (glyph + word) so severity survives
    # greyscale, colour-blindness, and print where backgrounds may drop.
    _SEV_RANK = {"card-red": 0, "card-no-data": 1, "card-orange": 2, "card-amber": 3, "card-green": 4}
    _STATUS_TOKEN = {
        "card-green":   "✓ Concordant",
        "card-amber":   "≈ Majority",
        "card-orange":  "≈ Partial",
        "card-red":     "✗ Discordant",
        "card-no-data": "▢ No call",
    }
    # Actionability-first ordering: lead with abnormal/actionable phenotypes
    # (reusing the clinical-findings tier) regardless of concordance, so a
    # concordant poor-metaboliser is not buried below a benign majority call.
    # Groups: ★ Actionable findings -> ⚠ Needs review -> ✓ Normal.
    _NEEDS_REVIEW = ("card-red", "card-no-data")
    _TIER_RANK = {"high": 0, "moderate": 1, "informational": 2}
    for g in genes_data:
        g["_tier"] = _get_tier(g["gene"], _pheno_cat(g["consensus_phenotype"]),
                               g["consensus_diplotype"], g.get("all_tool_diplotypes"))

    def _grp(g):
        if g["_tier"]:                              return 0   # actionable finding
        if g["card_class"] in _NEEDS_REVIEW:        return 1   # no confident verdict
        return 2                                               # normal

    _sorted_genes = sorted(genes_data, key=lambda g: (
        _grp(g), _TIER_RANK.get(g["_tier"], 9), _SEV_RANK.get(g["card_class"], 5), g["gene"]))
    _grp_counts = Counter(_grp(g) for g in genes_data)
    _GRP_LABEL = {0: ("★ Actionable findings", "grid-divider-flagged"),
                  1: ("⚠ Needs review",        "grid-divider-flagged"),
                  2: ("✓ Normal",              "grid-divider-normal")}
    _prev_grp = None
    for _i, gd in enumerate(_sorted_genes):
        gene = gd["gene"]
        dip  = gd["consensus_diplotype"]
        pheno = gd["consensus_phenotype"]
        card_class = gd["card_class"]
        n_agree = gd["n_agree"]
        n_called = gd["n_called"]
        badge_text = f"{n_agree}/{n_called} tools" if n_called else "No data"
        pill_color = phenotype_color(pheno)
        pheno_short = pheno if pheno != "-" else "—"

        card_class_map = {
            "card-green":   "card-green",
            "card-amber":   "card-amber",
            "card-orange":  "card-orange",
            "card-red":     "card-red",
            "card-no-data": "card-no-data",
        }

        badge_cls_map = {
            "card-green":   "badge-green",
            "card-amber":   "badge-amber",
            "card-orange":  "badge-orange",
            "card-red":     "badge-red",
            "card-no-data": "badge-grey",
        }
        badge_cls = badge_cls_map.get(card_class, "badge-grey")
        css_cls    = card_class_map.get(card_class, "")

        # Authority / resolution-tier provenance (why this is concordant despite
        # differing diplotype strings) — visible for CAP/CLIA review.
        _auth = gd.get("authority")
        if gd.get("phenotype_tier"):
            auth_badge = '<div class="gene-authority"><span class="badge badge-blue">Resolved by phenotype concordance</span></div>'
        elif _auth and _auth != "-":
            auth_badge = f'<div class="gene-authority"><span class="badge badge-blue">{_auth} authoritative</span></div>'
        else:
            auth_badge = ""

        # PharmCAT cross-check (informational; never changes the verdict). Green
        # when its phenotype agrees with the verdict, amber when it differs (review).
        _xc = (gd.get("cross_check") or {}).get("PharmCAT")
        if _xc:
            _xc_cls, _xc_mark = (("badge-green", "✓") if _xc.get("agrees")
                                 else ("badge-amber", "⚠ differs"))
            xcheck_badge = (f'<div class="gene-xcheck"><span class="badge {_xc_cls}">'
                            f'PharmCAT {_xc["diplotype"]} {_xc_mark}</span></div>')
        else:
            xcheck_badge = ""

        depth_info = gene_depth_map.get(gene, {})
        if depth_info:
            if depth_info.get("note") == "alt_contig":
                depth_html = '<div class="gene-depth">Depth: Alt contig</div>'
            else:
                mean  = depth_info.get("mean", "—")
                pct30 = depth_info.get("pct_ge_30x", "—")
                p30_cls = _depth_css(pct30)
                p30_str = f'<span class="{p30_cls}">{pct30}%</span>' if p30_cls else f'{pct30}%'
                depth_html = f'<div class="gene-depth">Depth: {mean}X &nbsp;|&nbsp; &#8805;30X: {p30_str}</div>'
        else:
            depth_html = ""

        if gene_fragments is not None:
            card_onclick = f'onclick="pgxShowGene(\'{gene}\'); return false;"'
            card_href = "#"
        else:
            card_onclick = ""
            if genes_rel_prefix:
                card_href = f"{genes_rel_prefix}/{gene}/{sample}.{gene}.pgx.html"
            else:
                card_href = f"{sample}.{gene}.pgx.html"

        _g = _grp(gd)
        if _g != _prev_grp:
            _lbl, _cls = _GRP_LABEL[_g]
            gene_cards_html += f'<div class="grid-divider {_cls}">{_lbl} ({_grp_counts[_g]})</div>'
            _prev_grp = _g
        gene_cards_html += f"""
            <a href="{card_href}" {card_onclick} class="gene-card {css_cls}">
                <div class="gene-name">{gene}</div>
                <div class="gene-status">{_STATUS_TOKEN.get(card_class, "")}</div>
                <div class="gene-diplo" title="{dip}">{dip if dip != '-' else '—'}</div>
                <div>
                    <span class="gene-pheno-pill" style="background:{pill_color}">{pheno_short}</span>
                </div>
                <div class="gene-concord">
                    Concordance: <span class="badge {badge_cls}">{badge_text}</span>
                </div>
                {auth_badge}
                {xcheck_badge}
                {depth_html}
            </a>"""

    bam_basename = os.path.basename(bam) if bam else "—"

    # At-a-glance headline counts — the 1-second clinician scan before the cards.
    _counts = Counter(gd["card_class"] for gd in genes_data)
    _total = len(genes_data)
    # need review = genes with no confident verdict (discordant / no-call)
    _flagged = sum(_counts.get(c, 0) for c in ("card-red", "card-no-data"))
    _flag_html = (f'<span class="summary-flagged">{_flagged} need review</span>'
                  if _flagged else '<span class="summary-flagged summary-ok">all resolved</span>')
    _pills = "".join(
        f'<span class="summary-pill {_cc}"><b>{_counts[_cc]}</b> {_g} {_lbl}</span>'
        for _cc, _g, _lbl in (("card-red", "✗", "Discordant"), ("card-no-data", "▢", "No call"),
                              ("card-amber", "≈", "Majority"), ("card-green", "✓", "Concordant"))
        if _counts.get(_cc, 0)
    )
    summary_strip_html = (f'<div class="summary-strip"><span class="summary-total">'
                          f'{_total} genes assessed</span>{_flag_html}{_pills}</div>')

    clinical_html = build_clinical_findings_section(
        genes_data, sample, genes_rel_prefix,
        embedded=(gene_fragments is not None))

    qc_html = ""
    gene_depth_html = ""
    if bs:
        qc_html = f"""
        <div class="section">
            <h2>BAM Quality Control</h2>
            {bam_stats_cards(bs)}
        </div>"""
        gene_depth_html = gene_depth_table(bs)

    today = date.today().isoformat()

    # ── CAP/CLIA lab identification blocks ────────────────────────────────────
    li               = lab_info or {}
    lab_name         = li.get("lab_name", "")
    clia_number      = li.get("clia_number", "")
    cap_number       = li.get("cap_number", "")
    medical_director = li.get("medical_director", "")
    accession_id     = li.get("accession_id", "")
    authorized_by    = li.get("authorized_by", "")
    report_status    = li.get("report_status", "RESEARCH USE ONLY")

    if lab_name or clia_number:
        _rows = ""
        if lab_name:          _rows += f'<div class="lab-name">{lab_name}</div>'
        if clia_number:       _rows += f'<div class="lab-detail">CLIA: {clia_number}</div>'
        if cap_number:        _rows += f'<div class="lab-detail">CAP: {cap_number}</div>'
        if medical_director:  _rows += f'<div class="lab-detail">Director: {medical_director}</div>'
        lab_header_html = f'<div class="lab-block">{_rows}</div>'
    else:
        lab_header_html = ""

    _rs = report_status.upper()
    _status_cls = "status-final" if _rs == "FINAL" else ("status-prelim" if _rs == "PRELIMINARY" else "status-research")
    status_banner_html = f'<div class="report-status-banner {_status_cls}">{report_status}</div>'

    accession_html = (
        f'<div class="accession-row"><span class="qc-label">Accession:</span>'
        f' <span class="accession-id">{accession_id}</span></div>'
    ) if accession_id else ""

    auth_html = (
        f'<div class="auth-block"><span class="auth-label">Authorized by:</span>'
        f' <span class="auth-name">{authorized_by}</span>'
        f' &nbsp;&mdash;&nbsp; <span class="auth-date">{today}</span></div>'
    ) if authorized_by else ""

    # Provenance footer (F2) — proves what produced this report. Optional: only
    # rendered when a provenance.json was passed in.
    provenance_html = ""
    if provenance:
        ref = provenance.get("reference", {}) or {}
        tools = provenance.get("tools", {}) or {}
        tool_bits = " &bull; ".join(f"{k} {v}" for k, v in tools.items() if v)
        rows = [
            ("Pipeline",   provenance.get("pipeline", "")),
            ("Generated",  provenance.get("generated_utc", "")),
            ("Reference",  f'{ref.get("path","")} ({ref.get("n_contigs","?")} contigs, {ref.get("bytes","?")} bytes)' if ref else ""),
            ("Tools",      tool_bits),
        ]
        cells = "".join(
            f'<div><span class="prov-k">{k}:</span> <span class="prov-v">{escape(str(v))}</span></div>'
            for k, v in rows if v
        )
        provenance_html = (
            '<div class="provenance" style="margin-top:0.6rem;font-size:0.78rem;'
            'color:#475569;text-align:left;display:inline-block">'
            '<div style="font-weight:600;color:#334155">Provenance</div>'
            f'{cells}</div>'
        )

    # Orthogonal validation panel — optional. Rendered only when an Axiom (or
    # other array) concordance table was passed in. Sample-specific external truth.
    axiom_html = ""
    if axiom:
        _glyph = {"MATCH": ("&#10003;", "Concordant", "#1a7f37"),
                  "PARTIAL": ("&#9680;", "Partial", "#9a6700"),
                  "MISMATCH": ("&#9888;", "Review", "#cf222e")}
        # Adjudication direction takes priority over raw array agreement: an
        # array false-negative is the array's gap, not a pipeline error, so it
        # must NOT render as a red "Review". Glyph chosen for colourblind safety
        # (shape + word, not hue alone).
        _dir_glyph = {
            "concordant":          ("&#10003;", "Concordant",          "#1a7f37"),
            "array_fn":            ("&#9650;",  "NGS — array gap", "#0969da"),
            "array_authoritative": ("&#9670;",  "Array (single-SNP)",  "#6e7781"),
            "array_fp_or_ngs_fn":  ("&#9888;",  "Confirm",             "#bf3989"),
            "ngs_unresolved":      ("&#9888;",  "Review",              "#cf222e"),
            "review":              ("&#9680;",  "Review",              "#9a6700"),
        }
        # Optional second array-software column (e.g. Axiom P6 vs P9). Detected from
        # the TSV header: if Axiom_P9 is present, render two array columns, else fall
        # back to the single "Axiom" column (EQ_2017-style TSVs stay valid).
        two_arrays = any(("Axiom_P9" in r) for r in axiom)
        body = ""
        for r in axiom:
            g, lbl, col = _dir_glyph.get((r.get("Direction") or "").lower()) or \
                _glyph.get(r.get("Match", "").upper(),
                           ("&bull;", r.get("Match", "?"), "#57606a"))
            if two_arrays:
                arr_cells = (f'<td>{escape(r.get("Axiom_P6",""))}</td>'
                             f'<td>{escape(r.get("Axiom_P9",""))}</td>')
            else:
                arr_cells = f'<td>{escape(r.get("Axiom",""))}</td>'
            reason = escape(r.get("Reason", ""))
            reason_html = (f'<div style="font-size:0.78rem;color:#57606a">{reason}</div>'
                           if reason else "")
            body += (
                f'<tr><td style="font-weight:600">{escape(r.get("Gene",""))}</td>'
                f'{arr_cells}'
                f'<td>{escape(r.get("Consensus",""))}</td>'
                f'<td style="color:{col};min-width:9rem">{g} {lbl}{reason_html}</td></tr>'
            )
        arr_hdr = ('<th>Axiom P6</th><th>Axiom P9</th>' if two_arrays
                   else '<th>Axiom</th>')
        axiom_html = (
            '<div class="container" style="margin-top:1.5rem">'
            '<h2 style="margin-bottom:0.3rem">Orthogonal Validation &mdash; Axiom Array</h2>'
            '<p style="font-size:0.85rem;color:#57606a;margin-top:0">Pipeline diplotypes adjudicated '
            'against the Thermo Fisher Axiom microarray (a comparator, not ground truth). '
            '&#10003; Concordant; &#9650; NGS &mdash; array gap = off-panel variant the array cannot '
            'score, NGS authoritative; &#9670; Array (single-SNP) = locus the array genotypes directly, '
            'trust the array; &#9888; Confirm = array and NGS disagree, orthogonal confirmation needed; '
            '&#9680; Review = nomenclature or panel difference.</p>'
            '<table class="axiom-table" style="width:100%;border-collapse:collapse;font-size:0.9rem">'
            '<thead><tr style="text-align:left;border-bottom:2px solid #d0d7de">'
            f'<th>Gene</th>{arr_hdr}<th>Pipeline consensus</th><th>Agreement</th></tr></thead>'
            f'<tbody>{body}</tbody></table></div>'
        )

    # Build embedded gene panels HTML
    gene_panels_html = ""
    if gene_fragments:
        for gfragment_gene, fragment in gene_fragments.items():
            gene_panels_html += (
                f'<div id="gene-panel-{gfragment_gene}" class="gene-panel" hidden>\n'
                f'<div class="container">\n{fragment}\n</div>\n</div>\n'
            )

    # Include DETAIL_EXTRA_CSS when embedding gene panels
    extra_style = f"\n{DETAIL_EXTRA_CSS}" if gene_fragments else ""

    # JS show/hide for embedded mode
    embedded_js = ""
    if gene_fragments:
        embedded_js = """
<script>
function pgxShowGene(gene) {
    document.getElementById('main-view').hidden = true;
    document.querySelectorAll('.gene-panel').forEach(function(el) { el.hidden = true; });
    var panel = document.getElementById('gene-panel-' + gene);
    if (panel) { panel.hidden = false; }
    window.scrollTo(0, 0);
}
function pgxShowMain() {
    document.querySelectorAll('.gene-panel').forEach(function(el) { el.hidden = true; });
    document.getElementById('main-view').hidden = false;
    window.scrollTo(0, 0);
}
// Wide-table horizontal scrolling: drag-to-pan anywhere + a floating bottom
// scrollbar synced to whichever overflowing table is on screen.
// Runs on DOMContentLoaded — this <script> lives in <head>, so the body/tables
// don't exist yet at parse time.
document.addEventListener('DOMContentLoaded', function () {
    document.querySelectorAll('.table-scroll').forEach(function (el) {
        var down = false, startX = 0, startLeft = 0, moved = false;
        el.addEventListener('mousedown', function (e) {
            if (e.button !== 0 || e.target.closest('a, button, input, select')) return;
            down = true; moved = false; startX = e.pageX; startLeft = el.scrollLeft;
            el.classList.add('dragging');
        });
        window.addEventListener('mousemove', function (e) {
            if (!down) return;
            var dx = e.pageX - startX;
            if (Math.abs(dx) > 3) moved = true;
            el.scrollLeft = startLeft - dx;
            if (moved) e.preventDefault();
        });
        window.addEventListener('mouseup', function () { down = false; el.classList.remove('dragging'); });
        el.addEventListener('click', function (e) {
            if (moved) { e.preventDefault(); e.stopPropagation(); }
        }, true);
    });

    var bar = document.createElement('div');
    bar.className = 'floating-xscroll';
    var inner = document.createElement('div');
    bar.appendChild(inner);
    document.body.appendChild(bar);
    var active = null, syncing = false;

    function usable(el) {
        if (el.offsetParent === null) return false;            // hidden gene panel
        if (el.scrollWidth - el.clientWidth < 4) return false; // no overflow
        var r = el.getBoundingClientRect();
        // on screen, but its own bottom scrollbar is below the fold
        return r.top < window.innerHeight && r.bottom > 60 && r.bottom > window.innerHeight;
    }
    function refresh() {
        var cands = document.querySelectorAll('.table-scroll');
        active = null;
        for (var i = 0; i < cands.length; i++) { if (usable(cands[i])) { active = cands[i]; break; } }
        if (!active) { bar.style.display = 'none'; return; }
        var r = active.getBoundingClientRect();
        bar.style.display = 'block';
        bar.style.left = r.left + 'px';
        bar.style.width = r.width + 'px';
        inner.style.width = active.scrollWidth + 'px';
        syncing = true; bar.scrollLeft = active.scrollLeft; syncing = false;
    }
    bar.addEventListener('scroll', function () {
        if (syncing || !active) return;
        syncing = true; active.scrollLeft = bar.scrollLeft; syncing = false;
    });
    document.addEventListener('scroll', function () {
        if (active && !syncing) { syncing = true; bar.scrollLeft = active.scrollLeft; syncing = false; }
        refresh();
    }, true);  // capture: also fires for inner table scroll (drag/wheel)
    window.addEventListener('resize', refresh);
    setTimeout(refresh, 0);
});
</script>"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PGx Report &#8212; {sample}</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600&family=Source+Serif+4:opsz,wght@8..60,600;8..60,700&display=swap" rel="stylesheet">
    <style>
{SHARED_CSS}
{LANDING_EXTRA_CSS}{extra_style}
    </style>
{embedded_js}
</head>
<body>
{status_banner_html}
<header onclick="if(window.pgxShowMain)pgxShowMain();return false;" style="cursor:pointer" title="Back to sample summary">
    <div>
        <div class="logo">PGx Suite</div>
        <div class="subtitle">Pharmacogenomics Star-Allele Report</div>
    </div>
    <div class="spacer"></div>
    {lab_header_html}
    <div class="report-meta">
        <div><strong>Sample:</strong> {sample}</div>
        <div><strong>Date:</strong> {today}</div>
        <div>GRCh38 &bull; CPIC Guidelines</div>
    </div>
</header>

<div id="main-view">
<div class="container">

    <div class="sample-banner">
        <div>
            <div class="qc-label">Sample ID</div>
            <div class="sample-id">{sample}</div>
            {accession_html}
            <div class="sample-meta">BAM: {bam_basename}</div>
        </div>
        <div>
            <div class="qc-label">Tools</div>
            <div style="margin-top:0.3rem; display:flex; gap:0.5rem; flex-wrap:wrap;">
                <span class="badge badge-blue">PyPGx 0.26</span>
                <span class="badge badge-blue">Stargazer 2.0.3</span>
                <span class="badge badge-blue">Aldy 4.8</span>
                <span class="badge badge-blue">StellarPGx 1.2.7</span>
                <span class="badge badge-blue">Cyrius</span>
                <span class="badge badge-blue">PharmCAT 3.2.0</span>
                <span class="badge badge-blue">OptiType 1.3.5</span>
                <span class="badge badge-blue">mutserve 2.0.3</span>
            </div>
        </div>
    </div>

{summary_strip_html}

{clinical_html}

{qc_html}

    <div class="section">
        <h2>Gene Diplotype Summary</h2>
        <div class="legend-row">
            <span style="font-weight:600;color:var(--muted)">Verdict:</span>
            <span class="legend-item">
                <span class="legend-swatch" style="background:#38a169"></span>Concordant
            </span>
            <span class="legend-item">
                <span class="legend-swatch" style="background:#d69e2e"></span>Majority
            </span>
            <span class="legend-item">
                <span class="legend-swatch" style="background:#e53e3e"></span>Discordant
            </span>
            <span class="legend-item">
                <span class="legend-swatch" style="background:#a0aec0"></span>No call
            </span>
            <span class="legend-item" style="color:var(--muted)">— per-card <b>n/N</b> shows how many callers agreed</span>
        </div>
        <div class="legend-authority">
            <div class="legend-auth-title">When star-allele callers can't agree, a CPIC reference method sets the verdict (shown as a badge on the gene card):</div>
            <div class="legend-auth-row"><span class="badge badge-blue">PharmCAT authoritative</span><span>diplotype set by PharmCAT, the CPIC reference star-allele caller (UGT1A1, CYP2B6, CYP4F2).</span></div>
            <div class="legend-auth-row"><span class="badge badge-blue">VCF-Check authoritative</span><span>genotype read directly from the CPIC variant table (RYR1, VKORC1, G6PD, CACNA1S, IFNL3).</span></div>
            <div class="legend-auth-row"><span class="badge badge-blue">Cyrius authoritative</span><span>CYP2D6 structural variant / copy-number resolved by Cyrius.</span></div>
            <div class="legend-auth-row"><span class="badge badge-blue">Resolved by phenotype concordance</span><span>callers gave different diplotype strings but agree on the CPIC phenotype.</span></div>
        </div>
        <div class="gene-grid">
{gene_cards_html}
        </div>
    </div>

{gene_depth_html}

</div>
</div>

{axiom_html}

<!-- Gene detail panels (embedded, shown/hidden by JS) -->
{gene_panels_html}

<footer>
    {auth_html}
    <div>PGx Suite &bull; GRCh38 Reference &bull; CPIC Guidelines (cpicpgx.org) &bull; PharmVar (pharmvar.org)</div>
    <div>PyPGx 0.26 &bull; Stargazer 2.0.3 &bull; Aldy 4.8.3 &bull; StellarPGx 1.2.7 &bull; OptiType 1.3.5</div>
    <div style="margin-top:0.4rem;font-style:italic">For clinical decision support only &#8212; not a standalone diagnostic report. Results require clinical validation and interpretation by a qualified clinician.</div>
    {provenance_html}
</footer>
</body>
</html>"""

    out_path = os.path.join(out_dir, f"{sample}_pgx_report.html")
    with open(out_path, "w") as fh:
        fh.write(html)
    print(f"Landing page: {out_path}")


# ── Per-gene detail page ───────────────────────────────────────────────────────
DETAIL_EXTRA_CSS = """
.gene-header {
    background: white;
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1.25rem 1.5rem;
    margin-bottom: 2rem;
    display: flex;
    align-items: flex-start;
    gap: 2rem;
    flex-wrap: wrap;
}
.gene-title { font-size: 1.8rem; font-weight: 700; color: var(--primary); }
.detail-table-wrap {
    background: white;
    border: 1px solid var(--border);
    border-radius: 8px;
    overflow-x: auto;              /* horizontal scrollbar when many tool columns */
    -webkit-overflow-scrolling: touch;
}
/* Wide tables: grab-and-drag to pan horizontally (see wideScroll JS) + a scrollbar
   that floats at the viewport bottom so you never scroll to the table's end to use it. */
.detail-table-wrap, .var-table-wrap { cursor: grab; }
.table-scroll.dragging { cursor: grabbing; user-select: none; }
.floating-xscroll {
    position: fixed; bottom: 0; left: 0; height: 15px;
    overflow-x: auto; overflow-y: hidden;
    z-index: 50; display: none;
    background: rgba(247,250,252,0.96);
    border-top: 1px solid var(--border);
    box-shadow: 0 -2px 6px rgba(0,0,0,0.08);
}
.floating-xscroll > div { height: 1px; }
/* per-column floor so columns stay readable; the table then overflows the wrap
   (one column per caller — now up to 9) and a scrollbar appears at the bottom */
.detail-table th, .detail-table td { min-width: 120px; vertical-align: top; }
.detail-table th:first-child,
.detail-table td:first-child {
    position: sticky; left: 0; z-index: 1;   /* keep the field-name column in view */
    background: #f7fafc;
    font-weight: 600;
    color: var(--primary);
    min-width: 180px;
    border-right: 2px solid var(--border);
}
.detail-table th { font-size: 0.85rem; }
.variant-list { font-family: monospace; font-size: 0.78rem; }
.sv-note {
    background: #ebf8ff;
    border-left: 4px solid #3182ce;
    border-radius: 4px;
    padding: 0.6rem 0.9rem;
    font-size: 0.83rem;
    color: #2a4365;
    margin-top: 1.25rem;
}
/* ── Variant evidence subtable ─────────────────────────────────────── */
.variant-section {
    background: white;
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1.25rem 1.5rem;
    margin-top: 2rem;
}
.var-meta { display: flex; align-items: center; gap: 1.5rem; flex-wrap: wrap; margin-bottom: 0.75rem; }
.var-summary { font-size: 0.82rem; color: var(--muted); }
.var-legend { display: flex; gap: 1rem; font-size: 0.75rem; align-items: center; }
.vleg-swatch { display: inline-block; width: 12px; height: 12px; border-radius: 2px; margin-right: 3px; vertical-align: middle; }
.vleg-all   { background: #c6f6d5; }
.vleg-multi { background: #fefcbf; }
.vleg-single{ background: #edf2f7; }
.var-table-wrap { overflow-x: auto; border: 1px solid var(--border); border-radius: 6px; }
.var-table { min-width: 750px; font-size: 0.78rem; border-collapse: collapse; width: 100%; }
.var-table thead th {
    background: var(--primary);
    color: white;
    padding: 0.45rem 0.65rem;
    font-size: 0.72rem;
    text-transform: uppercase;
    letter-spacing: 0.05em;
    text-align: left;
    font-weight: 600;
    white-space: nowrap;
}
.var-table thead th.vtool-col { background: #2c5282; }
.var-table td { padding: 0.38rem 0.65rem; border-bottom: 1px solid var(--border); vertical-align: top; }
.var-table tr:last-child td { border-bottom: none; }
.vrow-all   td { background: #f0fff4; }
.vrow-multi td { background: #fffff0; }
.vrow-single td { background: inherit; }
.vpos-col { font-family: 'Consolas', 'Cascadia Code', monospace; white-space: nowrap; color: #4a5568; }
.vchange-col { font-family: 'Consolas', 'Cascadia Code', monospace; font-weight: 700; white-space: nowrap; }
.vrsid-col { white-space: nowrap; }
.veff-col { font-size: 0.73rem; max-width: 200px; }
.vtool-data { font-size: 0.73rem; font-family: 'Consolas', 'Cascadia Code', monospace; line-height: 1.6; }
.vnot-detected { color: #cbd5e0; text-align: center; font-size: 0.82rem; }
.vallele-tag {
    display: inline-block;
    background: #ebf8ff; color: #2b6cb0;
    padding: 0.05rem 0.35rem;
    border-radius: 4px;
    font-size: 0.68rem;
    font-weight: 700;
    margin-bottom: 0.1rem;
    font-family: 'Segoe UI', sans-serif;
}
.vcopy-tag {
    display: inline-block;
    background: #f0e6ff; color: #553c9a;
    padding: 0.05rem 0.35rem;
    border-radius: 4px;
    font-size: 0.68rem;
    font-weight: 600;
    margin-bottom: 0.1rem;
    font-family: 'Segoe UI', sans-serif;
}
.vgt-tag {
    display: inline-block;
    background: #fff5f5; color: #742a2a;
    padding: 0.05rem 0.35rem;
    border-radius: 4px;
    font-size: 0.68rem;
    font-weight: 600;
    font-family: 'Segoe UI', sans-serif;
}
.rsid-link { color: var(--primary-light); text-decoration: none; font-size: 0.75rem; }
.rsid-link:hover { text-decoration: underline; }
.var-count-cell { font-size: 0.78rem; color: var(--muted); }
.var-anchor-link { font-size: 0.72rem; color: var(--primary-light); margin-left: 0.4rem; }
.gene-page-nav { display: flex; gap: 1rem; margin-bottom: 1.5rem; flex-wrap: wrap; }
.gene-page-nav a { font-size: 0.85rem; color: var(--primary); text-decoration: none; padding: 0.35rem 0.8rem; border: 1px solid var(--primary-light); border-radius: 6px; transition: background 0.12s ease; }
.gene-page-nav a:hover { background: var(--primary-light); color: white; }
.gene-depth-detail { font-size: 0.82rem; color: var(--muted); margin-top: 0.4rem; }
.depth-sensitivity-warn { font-size: 0.78rem; color: #92400e; background: #fef3c7; border: 1px solid #f59e0b; border-radius: 4px; padding: 0.25rem 0.5rem; margin-top: 0.35rem; }
.gene-locus { font-size: 0.78rem; color: var(--muted); font-family: monospace; margin-top: 0.2rem; }
/* ── CPIC reference section (gene detail page) ── */
.cpic-section { background: white; border: 1px solid var(--border); border-radius: 8px; padding: 1.5rem; margin-top: 2rem; }
.cpic-header { display: flex; justify-content: space-between; align-items: flex-start; gap: 1.5rem; margin-bottom: 0.9rem; flex-wrap: wrap; }
.cpic-desc { font-size: 0.85rem; color: var(--muted); margin-top: 0.4rem; max-width: 720px; line-height: 1.5; }
.cpic-patient-note { background: #ebf8ff; border-left: 4px solid #3182ce; border-radius: 4px; padding: 0.6rem 0.9rem; font-size: 0.83rem; color: #2a4365; margin-bottom: 1rem; line-height: 1.5; }
.cpic-patient-note--high     { background: #fff5f5; border-left-color: #c53030; color: #742a2a; }
.cpic-patient-note--moderate { background: #fffaf0; border-left-color: #dd6b20; color: #7b341e; }
.cpic-drug-wrap { overflow-x: auto; border: 1px solid var(--border); border-radius: 6px; margin-bottom: 1rem; }
.cpic-drug-table { width: 100%; border-collapse: collapse; font-size: 0.83rem; }
.cpic-drug-table th { background: var(--primary); color: white; padding: 0.5rem 0.75rem; text-align: left; font-size: 0.78rem; text-transform: uppercase; letter-spacing: 0.05em; font-weight: 600; }
.cpic-drug-table td { padding: 0.5rem 0.75rem; border-bottom: 1px solid var(--border); vertical-align: top; font-variant-numeric: tabular-nums; }
.cpic-drug-table tr:last-child td { border-bottom: none; }
.cpic-level-badge { display: inline-block; padding: 0.1rem 0.45rem; border-radius: 4px; font-weight: 700; font-size: 0.75rem; }
.cpic-level-a { background: #c6f6d5; color: #276749; }
.cpic-level-b { background: #bee3f8; color: #2a4365; }
.cpic-footer { font-size: 0.78rem; color: var(--muted); margin-top: 0.5rem; }
.cpic-link { color: var(--primary-light); text-decoration: none; }
.cpic-link:hover { text-decoration: underline; }
"""


def fmt_value(key: str, val) -> str:
    """Format a field value for HTML display."""
    if val is None or val == "" or val == "-":
        return '<span class="dash">—</span>'
    if isinstance(val, list):
        if not val:
            return '<span class="dash">—</span>'
        return "<br>".join(f'<span class="variant-list">{v}</span>' for v in val)
    s = str(val)
    if s == "-":
        return '<span class="dash">—</span>'
    return s


# ── Variant harmonization ──────────────────────────────────────────────────────

def harmonize_variants(tools_data: dict) -> list[dict]:
    """Merge supporting_variants from all tools into a position-sorted union.

    Variants within ±1 bp with the same REF and ALT are treated as the same
    event (handles VCF coordinate convention differences between tools).

    Returns a list of cluster dicts, each with:
        pos    – canonical GRCh38 position (int, lowest seen)
        ref    – reference allele (str)
        alt    – alternate allele (str)
        rsid   – best rsID across tools (str, or "-")
        effect – best functional effect label across tools (str, or "-")
        tools  – dict {tool_name: variant_dict} for tools that detected this variant
    """
    collected: list[tuple[int, str, str, str, dict]] = []
    for tool in TOOLS:
        for v in tools_data.get(tool, {}).get("supporting_variants", []):
            try:
                pos = int(v.get("pos", 0))
            except (ValueError, TypeError):
                continue
            ref = (v.get("ref") or "").upper().strip()
            alt = (v.get("alt") or "").upper().strip()
            if not ref or not alt or ref == "-" or alt == "-":
                continue
            collected.append((pos, ref, alt, tool, v))

    if not collected:
        return []

    collected.sort(key=lambda x: (x[0], x[1], x[2]))

    clusters: list[dict] = []
    for pos, ref, alt, tool, v in collected:
        # Search for a matching cluster (same ref+alt, position within ±1 bp)
        merged = False
        for cl in reversed(clusters):
            if cl["pos"] < pos - 2:
                break          # sorted — no older cluster can be within ±1
            if cl["ref"] == ref and cl["alt"] == alt and abs(cl["pos"] - pos) <= 1:
                cl["pos"] = min(cl["pos"], pos)   # keep the lower coordinate
                # Accumulate best rsid and effect from any tool
                if cl["rsid"] in ("-", "", None) and (v.get("rsid") or "-") != "-":
                    cl["rsid"] = v["rsid"]
                if cl["effect"] in ("-", "", None) and (v.get("effect") or "-") != "-":
                    cl["effect"] = v["effect"]
                # First reporter wins within a tool (avoids duplicate Aldy rows)
                cl["tools"].setdefault(tool, v)
                merged = True
                break
        if not merged:
            clusters.append({
                "pos":    pos,
                "ref":    ref,
                "alt":    alt,
                "rsid":   (v.get("rsid") or "-"),
                "effect": (v.get("effect") or "-"),
                "tools":  {tool: v},
            })

    clusters.sort(key=lambda c: c["pos"])
    return clusters


def _tool_cell(tool: str, v: dict | None) -> str:
    """Return the <td> HTML for one tool's evidence at a variant locus."""
    if v is None:
        return '<td class="vnot-detected">—</td>'
    parts: list[str] = []
    if tool == "PyPGx":
        allele = (v.get("allele") or "-")
        af     = (v.get("af")     or "-")
        if allele != "-":
            parts.append(f'<span class="vallele-tag">{allele}</span>')
        if af != "-":
            parts.append(f"AF&nbsp;{af}")
    elif tool == "Stargazer":
        af    = (v.get("af")    or "-")
        depth = (v.get("depth") or "-")
        if af    != "-": parts.append(f"AF&nbsp;{af}")
        if depth != "-": parts.append(f"Dp&nbsp;{depth}")
    elif tool == "Aldy":
        allele = (v.get("allele") or "-")
        depth  = (v.get("depth")  or "-")
        rsid   = (v.get("rsid")   or "-")
        if allele != "-":
            parts.append(f'<span class="vcopy-tag">{allele}</span>')
        if depth  != "-": parts.append(f"Dp&nbsp;{depth}")
        if rsid   != "-":
            parts.append(
                f'<a href="https://www.ncbi.nlm.nih.gov/snp/{rsid}" '
                f'target="_blank" class="rsid-link">{rsid}</a>'
            )
    elif tool == "StellarPGx":
        gt = (v.get("gt") or "-")
        if gt != "-":
            parts.append(f'<span class="vgt-tag">GT&nbsp;{gt}</span>')
    inner = "<br>".join(parts) if parts else "✓"
    return f'<td class="vtool-data">{inner}</td>'


def render_variant_subtable(clusters: list[dict], n_tools_per_gene: dict | None = None,
                            id_prefix: str = "") -> str:
    """Return the full HTML for the harmonised cross-tool variant evidence section."""
    if not clusters:
        return ""

    n_total  = len(clusters)
    n_shared = sum(1 for c in clusters if len(c["tools"]) >= 2)

    summary = (
        f'<span class="var-summary">'
        f'{n_total} unique variant{"s" if n_total != 1 else ""} '
        f'&nbsp;·&nbsp; {n_shared} detected by ≥2 tools'
        f'</span>'
    )
    legend = """
        <span class="var-legend">
            <span class="vleg-swatch vleg-all"></span>All tools&ensp;
            <span class="vleg-swatch vleg-multi"></span>≥2 tools&ensp;
            <span class="vleg-swatch vleg-single"></span>1 tool only
        </span>"""

    tool_header_cells = "".join(
        f'<th class="vtool-col">{t}</th>' for t in TOOLS
    )

    rows = ""
    for cl in clusters:
        n_det = len(cl["tools"])
        row_cls = (
            "vrow-all"    if n_det == len(TOOLS) else
            "vrow-multi"  if n_det >= 2          else
            "vrow-single"
        )
        rsid = cl.get("rsid") or "-"
        rsid_html = (
            f'<a href="https://www.ncbi.nlm.nih.gov/snp/{rsid}" '
            f'target="_blank" class="rsid-link">{rsid}</a>'
            if rsid != "-" else '<span class="dash">—</span>'
        )
        eff = cl.get("effect") or "-"
        eff_html = eff if eff != "-" else '<span class="dash">—</span>'

        tool_cells = "".join(_tool_cell(t, cl["tools"].get(t)) for t in TOOLS)

        rows += (
            f'<tr class="{row_cls}">'
            f'<td class="vpos-col">{cl["pos"]}</td>'
            f'<td class="vchange-col">{cl["ref"]}&rarr;{cl["alt"]}</td>'
            f'<td class="vrsid-col">{rsid_html}</td>'
            f'<td class="veff-col">{eff_html}</td>'
            f'{tool_cells}'
            f'</tr>\n'
        )

    return f"""
    <div class="variant-section" id="{id_prefix}variant-evidence">
        <h2>Supporting Variant Evidence</h2>
        <div class="var-meta">{summary}{legend}</div>
        <div class="var-table-wrap table-scroll">
        <table class="var-table">
            <thead><tr>
                <th>Position (GRCh38)</th>
                <th>REF&rarr;ALT</th>
                <th>rsID</th>
                <th>Functional effect</th>
                {tool_header_cells}
            </tr></thead>
            <tbody>
{rows}            </tbody>
        </table>
        </div>
    </div>"""


def _build_gene_inner(sample: str, gene: str, detail: dict, gene_depth: dict | None,
                      back_href: str, id_prefix: str = "") -> str:
    """Return the inner HTML fragment for a gene detail panel (no html/head/body wrapper).

    id_prefix — prepended to all section IDs (e.g. 'CYP2D6-') to avoid collisions
                when multiple gene panels are embedded in the same document.
    back_href — href for the <- Back link; use 'javascript:void(0)' for embedded mode.
    """
    tools_data = detail.get("tools", {})
    sv_note = detail.get("sv_mode", "")

    # Fill the PharmCAT column for genes where it ran as an informational
    # cross-check (not the authority). Display-only: a synthesized copy keeps
    # the cross-check out of the vote math (which uses the verdict, or the
    # cross-check-free `tools_data` fallback below). ponytail: shallow copy.
    _ccp = (detail.get("cross_check") or {}).get("PharmCAT")
    if _ccp and "PharmCAT" not in tools_data:
        display_tools = dict(tools_data)
        display_tools["PharmCAT"] = {
            "diplotype": _ccp.get("diplotype", "-"),
            "phenotype": _ccp.get("phenotype", "-"),
        }
    else:
        display_tools = tools_data

    clusters = harmonize_variants(tools_data)
    var_subtable_html = render_variant_subtable(clusters, id_prefix=id_prefix)

    # Consensus comes from the authoritative verdict in detail.json (single
    # source). NO_CALL / DISCORDANT yield no asserted phenotype, so the CPIC
    # clinical section below cannot fire off a non-call.
    _v = detail.get("verdict")
    if _v:
        consensus_diplo, consensus_pheno, _cc, _na, _nc = verdict_card(_v)
        consensus_pheno = "" if consensus_pheno == "-" else consensus_pheno
        consensus_diplo = "" if consensus_diplo == "-" else consensus_diplo
    else:
        from collections import Counter as _Counter
        pheno_list = [tools_data.get(t, {}).get("phenotype", "-") for t in TOOLS
                      if tools_data.get(t, {}).get("phenotype", "-") not in ("-", "")]
        consensus_pheno = _Counter(pheno_list).most_common(1)[0][0] if pheno_list else ""
        diplo_list = [normalize_diplotype(tools_data.get(t, {}).get("diplotype", "-")) for t in TOOLS
                      if normalize_diplotype(tools_data.get(t, {}).get("diplotype", "-")) != "-"]
        consensus_diplo = _Counter(diplo_list).most_common(1)[0][0] if diplo_list else ""
    all_tool_diplos_page = [tools_data.get(t, {}).get("diplotype", "-") for t in TOOLS
                            if tools_data.get(t, {}).get("diplotype", "-") not in ("-", "")]
    cpic_section_html = build_gene_cpic_section(
        gene, consensus_pheno, consensus_diplo, all_tool_diplos_page, id_prefix=id_prefix)

    header_cells = "".join(f"<th>{t}</th>" for t in TOOLS)
    rows_html = ""
    for field_key, field_label in FIELDS:
        if field_key == "sv_mode":
            val_html = f'<td colspan="{len(TOOLS)}" style="color:var(--muted);font-style:italic">{fmt_value("sv_mode", sv_note)}</td>'
            rows_html += f"<tr><td>{field_label}</td>{val_html}</tr>\n"
        elif field_key == "supporting_variants":
            cells = ""
            for tool in TOOLS:
                vlist = display_tools.get(tool, {}).get("supporting_variants", [])
                n = len(vlist) if isinstance(vlist, list) else 0
                cells += f'<td><span class="var-count-cell">{n} variant{"s" if n != 1 else ""}</span></td>'
            anchor = f'<a href="#{id_prefix}variant-evidence" class="var-anchor-link">&#8595; see table below</a>'
            rows_html += f"<tr><td>{field_label}{anchor}</td>{cells}</tr>\n"
        else:
            cells = ""
            for tool in TOOLS:
                td = display_tools.get(tool, {})
                raw = td.get(field_key, "-")
                cells += f"<td>{fmt_value(field_key, raw)}</td>"
            rows_html += f"<tr><td>{field_label}</td>{cells}</tr>\n"

    diplotypes = []
    for tool in TOOLS:
        d = normalize_diplotype(tools_data.get(tool, {}).get("diplotype", "-"))
        if d and d != "-":
            diplotypes.append(d)
    n_called = len(diplotypes)
    n_agree = Counter(diplotypes).most_common(1)[0][1] if diplotypes else 0
    # Prefer the single-source verdict so the detail page matches the landing card
    # exactly (failed / no-call callers like Cyrius are excluded the same way,
    # and the colour follows the verdict status, not the raw tool fraction).
    if _v:
        n_agree, n_called, card_class = _na, _nc, _cc
        badge_text = f"{n_agree}/{n_called}" if n_called else "No data"
    else:
        card_class, badge_text = concordance_color(n_agree, n_called)
    badge_cls_map = {
        "card-green":   "badge-green",
        "card-amber":   "badge-amber",
        "card-orange":  "badge-orange",
        "card-red":     "badge-red",
        "card-no-data": "badge-grey",
    }
    badge_cls = badge_cls_map.get(card_class, "badge-grey")

    if gene_depth and gene_depth.get("note") == "alt_contig":
        depth_detail_html = '<div class="gene-depth-detail">Depth: Alt contig (chr22_KI270879v1_alt)</div>'
    elif gene_depth:
        mean  = gene_depth.get("mean", "—")
        pct30 = gene_depth.get("pct_ge_30x", "—")
        pct20 = gene_depth.get("pct_ge_20x", "—")
        p30_cls = _depth_css(pct30)
        p30_str = (f'<span class="{p30_cls}">{pct30}%</span>'
                   if p30_cls else f'<strong>{pct30}%</strong>')
        depth_detail_html = (
            f'<div class="gene-depth-detail">'
            f'Mean depth: <strong>{mean}X</strong>'
            f' &nbsp;|&nbsp; &#8805;20X: <strong>{pct20}%</strong>'
            f' &nbsp;|&nbsp; &#8805;30X: {p30_str}'
            f'</div>'
        )
        # Low-depth clinical warning for depth-sensitive genes
        _min_depth = _DEPTH_SENSITIVE_GENES.get(gene)
        try:
            _mean_val = float(mean)
        except (TypeError, ValueError):
            _mean_val = None
        if _min_depth and _mean_val is not None and _mean_val < _min_depth:
            depth_detail_html += (
                f'<div class="depth-sensitivity-warn">'
                f'&#9888; Mean gene depth {mean}X is below the recommended {_min_depth}X for '
                f'{gene}. Rare heterozygous variants may be missed — consider confirmatory testing.'
                f'</div>'
            )
    else:
        depth_detail_html = ""

    locus = GENE_LOCI.get(gene, "")
    locus_html = (
        f'<div class="gene-locus">&#128205; {locus} (GRCh38)</div>'
        if locus else ""
    )

    sv_note_html = f'<div class="sv-note">&#128202; {sv_note}</div>' if sv_note and sv_note != "-" else ""

    # PharmCAT cross-check (informational; not part of the verdict).
    _xc = (detail.get("cross_check") or {}).get("PharmCAT")
    if _xc:
        _xc_cls = "badge-green" if _xc.get("agrees") else "badge-amber"
        _xc_word = "agrees with verdict" if _xc.get("agrees") else "differs — review"
        xcheck_detail_html = (
            '<div class="xcheck-note">Cross-check &mdash; PharmCAT (CPIC reference, '
            'informational, not part of the verdict): '
            f'<span class="badge {_xc_cls}">{_xc["diplotype"]}</span> '
            f'{_xc.get("phenotype","-")} &middot; {_xc_word}</div>'
        )
    else:
        xcheck_detail_html = ""

    # back_onclick is set when back_href signals embedded mode
    if back_href.startswith("javascript:"):
        back_el = f'<a href="#" onclick="pgxShowMain(); return false;" class="back-link">&#8592; Back to sample summary</a>'
    else:
        back_el = f'<a href="{back_href}" class="back-link">&#8592; Back to sample summary</a>'

    return f"""
    {back_el}

    <nav class="gene-page-nav">
        <a href="#{id_prefix}tool-results">Tool Results</a>
        <a href="#{id_prefix}variant-evidence">Variant Evidence</a>
        <a href="#{id_prefix}cpic-reference">CPIC Reference</a>
    </nav>

    <div class="gene-header">
        <div>
            <div class="qc-label">Gene</div>
            <div class="gene-title">{gene}</div>
            {locus_html}
            <div class="sample-meta" style="margin-top:0.3rem">Sample: {sample} &bull; GRCh38</div>
        </div>
        <div>
            <div class="qc-label">Tool concordance</div>
            <div style="margin-top:0.4rem">
                <span class="badge {badge_cls}" style="font-size:0.9rem;padding:0.3rem 0.9rem">{badge_text}</span>
            </div>
            {depth_detail_html}
        </div>
    </div>
    {xcheck_detail_html}

    <div class="section" id="{id_prefix}tool-results">
        <h2>Tool Results &#8212; {gene}</h2>
        <div class="detail-table-wrap table-scroll">
            <table class="detail-table">
                <thead>
                    <tr>
                        <th>Field</th>
                        {header_cells}
                    </tr>
                </thead>
                <tbody>
{rows_html}                </tbody>
            </table>
        </div>
{sv_note_html}
    </div>

{var_subtable_html}

{cpic_section_html}
"""


def build_gene_page(sample: str, gene: str, detail: dict, landing_file: str, out_dir: str,
                    gene_depth: dict | None = None):
    """Build <sample>.<gene>.pgx.html detail page (standalone file, not used in default pipeline)."""
    today = date.today().isoformat()
    inner = _build_gene_inner(sample, gene, detail, gene_depth,
                              back_href=landing_file, id_prefix="")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PGx Report &#8212; {sample} &#8212; {gene}</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600&family=Source+Serif+4:opsz,wght@8..60,600;8..60,700&display=swap" rel="stylesheet">
    <style>
{SHARED_CSS}
{DETAIL_EXTRA_CSS}
    </style>
</head>
<body>
<header>
    <div>
        <div class="logo">PGx Suite</div>
        <div class="subtitle">Pharmacogenomics Star-Allele Report</div>
    </div>
    <div class="spacer"></div>
    <div class="report-meta">
        <div><strong>Sample:</strong> {sample} &bull; <strong>Gene:</strong> {gene}</div>
        <div><strong>Date:</strong> {today} &bull; GRCh38</div>
    </div>
</header>

<div class="container">
{inner}
</div>

<footer>
    <div>PGx Suite &bull; GRCh38 Reference &bull; CPIC Guidelines (cpicpgx.org)</div>
    <div>PyPGx 0.26 &bull; Stargazer 2.0.3 &bull; Aldy 4.8.3 &bull; StellarPGx 1.2.7 &bull; OptiType 1.3.5</div>
    <div style="margin-top:0.4rem;font-style:italic">For clinical decision support only &#8212; not a standalone diagnostic report.</div>
</footer>
</body>
</html>"""

    out_path = os.path.join(out_dir, f"{sample}.{gene}.pgx.html")
    with open(out_path, "w") as fh:
        fh.write(html)
    print(f"  {gene}: {out_path}")



# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description="Generate PGx HTML reports")
    ap.add_argument("--sample",           required=True, help="Sample ID")
    ap.add_argument("--output",           required=True, help="Root output directory (landing page written here)")
    ap.add_argument("--genes-dir",        default="",    help="Directory containing per-gene subdirs (default: <output>/Genes)")
    ap.add_argument("--bam",              default="",    help="BAM file path (for display)")
    ap.add_argument("--bam-stats",        default="",    help="Path to bam_stats.json (default: <output>/log/bam_stats.json)")
    # CAP/CLIA Phase 1 compliance fields (all optional)
    ap.add_argument("--lab-name",         default="",                  help="Laboratory name")
    ap.add_argument("--clia-number",      default="",                  help="CLIA certificate number")
    ap.add_argument("--cap-number",       default="",                  help="CAP accreditation number")
    ap.add_argument("--medical-director", default="",                  help="Medical/laboratory director name and credentials")
    ap.add_argument("--accession-id",     default="",                  help="Unique report/accession number")
    ap.add_argument("--authorized-by",    default="",                  help="Authorizing pathologist or clinical geneticist")
    ap.add_argument("--report-status",    default="RESEARCH USE ONLY", help="Report status: FINAL, PRELIMINARY, or RESEARCH USE ONLY (default)")
    ap.add_argument("--provenance",       default="",                  help="Path to provenance.json (tool versions/reference/timestamp; rendered in the report footer)")
    ap.add_argument("--axiom",            default="",                  help="Path to an orthogonal-validation TSV (cols: Gene, Axiom, Consensus, Match); renders an Axiom-array concordance panel")
    args = ap.parse_args()

    sample  = args.sample
    out_dir = args.output
    genes_dir = args.genes_dir or os.path.join(out_dir, "Genes")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(genes_dir, exist_ok=True)

    # Load BAM stats
    bs = None
    bam_stats_path = args.bam_stats
    if not bam_stats_path:
        # bamstats rule writes <output>/bam_stats.json; older layouts used log/.
        for _cand in (os.path.join(out_dir, "bam_stats.json"),
                      os.path.join(out_dir, "log", "bam_stats.json")):
            if os.path.isfile(_cand):
                bam_stats_path = _cand
                break
        else:
            bam_stats_path = os.path.join(out_dir, "bam_stats.json")
    if os.path.isfile(bam_stats_path):
        with open(bam_stats_path) as fh:
            bs = json.load(fh)
    else:
        print(f"[warn] bam_stats.json not found at {bam_stats_path}; skipping BAM QC section",
              file=sys.stderr)

    bam_path = args.bam or (bs.get("bam", "") if bs else "")
    gene_depth_map = (bs or {}).get("gene_depth", {})

    # Discover genes from summary TSV (written to <output>/log/ by pgx-all-genes.sh)
    summary_tsv = os.path.join(out_dir, "log", "all_genes_summary.tsv")
    if not os.path.isfile(summary_tsv):
        # Fallback: legacy location at output root
        summary_tsv = os.path.join(out_dir, "all_genes_summary.tsv")
    if not os.path.isfile(summary_tsv):
        print(f"ERROR: all_genes_summary.tsv not found in {out_dir}/log/ or {out_dir}/", file=sys.stderr)
        sys.exit(1)

    seen_genes: list[str] = []
    with open(summary_tsv) as fh:
        for line in fh:
            if line.startswith("Gene"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0] and parts[0] not in seen_genes:
                seen_genes.append(parts[0])

    # Build gene fragments and collect landing data
    print(f"Collecting gene detail panels for {len(seen_genes)} genes …")
    genes_data = []
    gene_fragments: dict = {}

    for gene in seen_genes:
        # Find detail JSON — new layout: genes_dir/<gene>/; fallback: out_dir/<gene>/
        detail_json = os.path.join(genes_dir, gene, f"{gene}_{sample}_detail.json")
        if not os.path.isfile(detail_json):
            matches = glob.glob(os.path.join(genes_dir, gene, "*_detail.json"))
            if not matches:
                # Legacy fallback
                matches = glob.glob(os.path.join(out_dir, gene, "*_detail.json"))
            detail_json = matches[0] if matches else None

        if detail_json and os.path.isfile(detail_json):
            with open(detail_json) as fh:
                detail = json.load(fh)
        else:
            print(f"  [warn] no detail JSON for {gene}", file=sys.stderr)
            detail = {"gene": gene, "sample": sample, "tools": {}}

        tools_data = detail.get("tools", {})
        _v = detail.get("verdict")
        if _v:
            consensus_dip, consensus_pheno, card_class, n_agree, n_called = verdict_card(_v)
        else:
            # Backward-compat: detail.json predates the verdict authority.
            consensus_dip, consensus_pheno, card_class, n_agree = compute_concordance(tools_data)
            n_called = sum(
                1 for t in TOOLS
                if normalize_diplotype(tools_data.get(t, {}).get("diplotype", "-")) != "-"
            )

        # Collect raw per-tool diplotypes (unnormalised) for diplotype_check genes
        all_tool_diplotypes = [
            tools_data.get(t, {}).get("diplotype", "-")
            for t in TOOLS
            if tools_data.get(t, {}).get("diplotype", "-") not in ("-", "")
        ]

        genes_data.append({
            "gene": gene,
            "consensus_diplotype": consensus_dip,
            "consensus_phenotype": consensus_pheno,
            "card_class": card_class,
            "n_agree": n_agree,
            "n_called": n_called,
            "all_tool_diplotypes": all_tool_diplotypes,
            "authority": (_v or {}).get("authority"),
            "phenotype_tier": bool((_v or {}).get("phenotype_tier")),
            "cross_check": detail.get("cross_check", {}),
        })

        fragment = _build_gene_inner(
            sample, gene, detail, gene_depth_map.get(gene),
            back_href="javascript:void(0)", id_prefix=f"{gene}-")
        gene_fragments[gene] = fragment
        print(f"  {gene}: fragment built")

    # Build landing page at <out_dir>/<sample>_pgx_report.html
    print(f"Generating standalone HTML report …")
    lab_info = {
        "lab_name":         args.lab_name,
        "clia_number":      args.clia_number,
        "cap_number":       args.cap_number,
        "medical_director": args.medical_director,
        "accession_id":     args.accession_id,
        "authorized_by":    args.authorized_by,
        "report_status":    args.report_status,
    }
    prov = None
    if args.provenance and os.path.isfile(args.provenance):
        try:
            with open(args.provenance) as fh:
                prov = json.load(fh)
        except (OSError, json.JSONDecodeError) as e:
            print(f"  warning: could not read provenance {args.provenance}: {e}")
    axiom_rows = None
    if args.axiom and os.path.isfile(args.axiom):
        import csv
        with open(args.axiom, newline="") as fh:
            axiom_rows = list(csv.DictReader(fh, delimiter="\t"))
    build_landing(sample, bam_path, genes_data, bs, out_dir,
                  genes_rel_prefix="Genes",
                  gene_fragments=gene_fragments,
                  lab_info=lab_info, provenance=prov, axiom=axiom_rows)
    print("Done.")


if __name__ == "__main__":
    main()
