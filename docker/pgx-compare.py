#!/usr/bin/env python3
"""
pgx-compare.py — Parse outputs from all 4 PGx callers and produce a comparison
table (TSV) and a full-detail JSON used by pgx-report.py for HTML reporting.

Usage: pgx-compare.py --gene GENE --sample SAMPLE --output-dir DIR
"""
import argparse
import csv
import glob
import io
import json
import os
import re
import sys
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any


@dataclass
class CallerResult:
    tool: str
    # ── Fields 3-7: allele calls ──────────────────────────────────────────────
    diplotype: str = "-"
    haplotype1: str = "-"
    haplotype2: str = "-"
    sub_alleles: str = "-"
    alternative_diplotypes: str = "-"
    # ── Fields 8-9: clinical ──────────────────────────────────────────────────
    phenotype: str = "-"
    activity_score: str = "-"
    # ── Fields 10-11: structural variants ─────────────────────────────────────
    sv_type: str = "-"
    copy_number: str = "-"
    # ── Fields 12-16: supporting evidence ────────────────────────────────────
    supporting_variants: list = field(default_factory=list)
    functional_effects: str = "-"
    dbsnp_ids: str = "-"
    allele_score: str = "-"
    mean_af: str = "-"
    # ── Field 17: method ──────────────────────────────────────────────────────
    phasing_method: str = "-"
    # ── Operational ───────────────────────────────────────────────────────────
    status: str = "not_run"   # not_run | ok | failed

    def to_dict(self, sample: str, gene: str) -> dict[str, Any]:
        """Serialise to the 17-field dict written into the detail JSON."""
        return {
            "sample_id":              sample,
            "gene":                   gene,
            "diplotype":              self.diplotype,
            "haplotype1":             self.haplotype1,
            "haplotype2":             self.haplotype2,
            "sub_alleles":            self.sub_alleles,
            "alternative_diplotypes": self.alternative_diplotypes,
            "phenotype":              self.phenotype,
            "activity_score":         self.activity_score,
            "sv_type":                self.sv_type,
            "copy_number":            self.copy_number,
            "supporting_variants":    self.supporting_variants,
            "functional_effects":     self.functional_effects,
            "dbsnp_ids":              self.dbsnp_ids,
            "allele_score":           self.allele_score,
            "mean_af":                self.mean_af,
            "phasing_method":         self.phasing_method,
            "status":                 self.status,
        }


# ── SV gene sets ──────────────────────────────────────────────────────────────
PYPGX_SV_GENES = frozenset(
    {"CYP2A6", "CYP2B6", "CYP2D6", "CYP2E1", "CYP4F2", "G6PD", "GSTM1", "GSTT1"}
)
STARGAZER_SV_GENES = frozenset({"CYP2A6", "CYP2B6", "CYP2D6"})

# ── Gene support matrix ───────────────────────────────────────────────────────
GENE_SUPPORT: dict[str, dict[str, bool]] = {
    "CYP2D6":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP2C19":  {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP2C9":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP2B6":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP2C8":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP3A4":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP3A5":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP4F2":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "NUDT15":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "TPMT":     {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "UGT1A1":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "SLCO1B1":  {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "DPYD":     {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": False},
    "NAT1":     {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "NAT2":     {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "G6PD":     {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": False},
    "GSTM1":    {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "GSTT1":    {"pypgx": True,  "stargazer": False, "aldy": False, "stellarpgx": True},
    "POR":      {"pypgx": True,  "stargazer": True,  "aldy": False, "stellarpgx": True},
    "CYPOR":    {"pypgx": True,  "stargazer": True,  "aldy": False, "stellarpgx": True},
    "VKORC1":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": False},
    "CYP1A1":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP1A2":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP2A6":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "CYP2E1":   {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": True},
    "IFNL3":    {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": False},
    "RYR1":     {"pypgx": True,  "stargazer": True,  "aldy": True,  "stellarpgx": False},
    "ABCG2":    {"pypgx": False, "stargazer": False, "aldy": True,  "stellarpgx": True,  "optitype": False},
    "HLA-A":    {"pypgx": False, "stargazer": False, "aldy": False, "stellarpgx": False, "optitype": True},
    "HLA-B":    {"pypgx": False, "stargazer": False, "aldy": False, "stellarpgx": False, "optitype": True},
    "CACNA1S":  {"pypgx": True,  "stargazer": False, "aldy": False, "stellarpgx": True},
    "MT-RNR1":  {"pypgx": False, "stargazer": False, "aldy": False, "stellarpgx": False, "mutserve": True},
}


# ── Helpers ───────────────────────────────────────────────────────────────────
def _dash(v: str | None) -> str:
    return (v or "").strip() or "-"


def _parse_pypgx_variant_data(vdata: str) -> tuple[list[dict], str, str]:
    """Parse PyPGx VariantData field.

    Format: *allele:chrom-pos-ref-alt[,pos-ref-alt,...]:AF[,AF,...];...
    Returns (variants_list, functional_effects, dbsnp_ids).
    PyPGx does not report rsIDs or effect labels in data.tsv.
    """
    variants: list[dict] = []
    if not vdata or vdata == "-":
        return variants, "-", "-"
    for entry in vdata.rstrip(";").split(";"):
        entry = entry.strip()
        if not entry:
            continue
        m = re.match(r"^(\S+?):(.+)$", entry)
        if not m:
            continue
        allele = m.group(1)
        rest = m.group(2)
        # rest = "chrom-pos-ref-alt[,...]:AF[,...]"
        parts = rest.rsplit(":", 1)
        loci_str = parts[0]
        afs_str = parts[1] if len(parts) > 1 else ""
        afs = [a.strip() for a in afs_str.split(",")]
        loci = [l.strip() for l in loci_str.split(",")]
        for i, locus in enumerate(loci):
            lm = re.match(r"(\w+)-(\d+)-([A-Z]+)-([A-Z<>.]+)", locus)
            if lm:
                variants.append({
                    "allele": allele,
                    "chrom":  lm.group(1),
                    "pos":    lm.group(2),
                    "ref":    lm.group(3),
                    "alt":    lm.group(4),
                    "af":     afs[i] if i < len(afs) else "-",
                    "depth":  "-",
                    "effect": "-",
                    "rsid":   "-",
                })
    return variants, "-", "-"


def _parse_stargazer_core(core_str: str) -> tuple[list[dict], list[str]]:
    """Parse Stargazer hap_main_core field.

    Format: <pos:ref>alt:reads_ref/reads_alt:AF:effect_type:impact:protein_change>,...
    Returns (variants, effect_labels).
    """
    variants: list[dict] = []
    effects: list[str] = []
    if not core_str or core_str in ("-", "."):
        return variants, effects
    for token in core_str.split(">,"):
        token = token.strip().lstrip("<").rstrip(">")
        if not token:
            continue
        # <42126611:C>G:37/37:1.00:missense_variant:low_impact:S486T>
        m = re.match(
            r"(\d+):([A-Z]+)>([A-Z.]+):(\d+)/(\d+):([\d.]+)"
            r"(?::([^:]+))?(?::([^:]+))?(?::([^:>]+))?",
            token,
        )
        if m:
            depth_ref, depth_alt = int(m.group(4)), int(m.group(5))
            total = depth_ref + depth_alt
            effect = m.group(7) or "-"
            protein = m.group(9) or "-"
            if effect != "-":
                effects.append(f"{protein} ({effect})" if protein != "-" else effect)
            variants.append({
                "allele": "-",
                "chrom":  "-",
                "pos":    m.group(1),
                "ref":    m.group(2),
                "alt":    m.group(3),
                "af":     m.group(6),
                "depth":  str(total),
                "effect": f"{protein} ({effect})" if protein != "-" else effect,
                "rsid":   "-",
            })
    return variants, effects


# ── Parser: PyPGx ─────────────────────────────────────────────────────────────
def parse_pypgx(output_dir: str, gene: str, sample: str) -> CallerResult:
    result = CallerResult(tool="PyPGx", phasing_method="Beagle (1KGP reference panel)")

    candidates = [os.path.join(output_dir, "pypgx", "results.zip")]
    candidates += glob.glob(
        os.path.join(output_dir, "pypgx", "**", "results.zip"), recursive=True
    )
    zip_path = next((p for p in candidates if os.path.exists(p)), None)
    if zip_path is None:
        return result

    try:
        with zipfile.ZipFile(zip_path) as zf:
            tsv_name = next(
                (n for n in zf.namelist() if n.endswith("data.tsv")), None
            )
            if tsv_name is None:
                result.status = "failed"
                result.diplotype = "no data.tsv in results.zip"
                return result

            with zf.open(tsv_name) as f:
                reader = csv.DictReader(io.TextIOWrapper(f), delimiter="\t")
                for row in reader:
                    result.diplotype = _dash(
                        row.get("Genotype") or row.get("Diplotype") or row.get("Haplotype")
                    )
                    result.activity_score = _dash(
                        row.get("ActivityScore") or row.get("Activity_Score")
                    )
                    result.phenotype = _dash(row.get("Phenotype"))
                    result.haplotype1 = _dash(row.get("Haplotype1"))
                    result.haplotype2 = _dash(row.get("Haplotype2"))

                    # Sub-alleles: Haplotype2 often contains phased sub-alleles (;-separated)
                    h2 = row.get("Haplotype2") or ""
                    parts = [p for p in h2.split(";") if p.strip()]
                    if len(parts) > 1:
                        result.sub_alleles = "; ".join(parts)

                    result.alternative_diplotypes = _dash(row.get("AlternativePhase"))
                    result.sv_type = _dash(row.get("CNV"))
                    result.copy_number = "2" if result.sv_type == "Normal" else "-"

                    vdata = row.get("VariantData") or ""
                    result.supporting_variants, result.functional_effects, result.dbsnp_ids = \
                        _parse_pypgx_variant_data(vdata)

                    result.allele_score = "-"
                    result.mean_af = "-"
                    if result.supporting_variants:
                        afs = [v["af"] for v in result.supporting_variants if v["af"] != "-"]
                        if afs:
                            result.mean_af = f"per-variant AF: {', '.join(afs[:4])}" + \
                                             ("…" if len(afs) > 4 else "")

                    result.status = "ok"
                    break
    except Exception as exc:
        result.status = "failed"
        result.diplotype = f"parse error: {exc}"

    return result


# ── Parser: Stargazer ─────────────────────────────────────────────────────────
def parse_stargazer(output_dir: str, gene: str, sample: str) -> CallerResult:
    result = CallerResult(
        tool="Stargazer",
        phasing_method="Beagle statistical phasing + SSR CN marker",
    )
    gene_lower = gene.lower()

    candidates = [
        os.path.join(output_dir, "stargazer", "genotype-calls.tsv"),
        os.path.join(output_dir, "stargazer", "genotype.txt"),
        os.path.join(output_dir, "stargazer", f"{gene_lower}.genotype.txt"),
        os.path.join(output_dir, "stargazer", f"{gene}.genotype.txt"),
    ]
    candidates += glob.glob(
        os.path.join(output_dir, "stargazer", "**", "genotype*.tsv"), recursive=True
    )
    candidates += glob.glob(
        os.path.join(output_dir, "stargazer", "**", "genotype*.txt"), recursive=True
    )
    geno_path = next((p for p in candidates if os.path.exists(p)), None)
    if geno_path is None:
        return result

    try:
        with open(geno_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                row = {k.lstrip("#").strip(): v for k, v in row.items()}

                hap1 = _dash(row.get("hap1_main"))
                hap2 = _dash(row.get("hap2_main"))
                if hap1 != "-" and hap2 != "-":
                    result.diplotype = f"{hap1}/{hap2}"
                else:
                    result.diplotype = _dash(
                        row.get("Genotype") or row.get("genotype") or row.get("Diplotype")
                    )

                raw_score = row.get("dip_score") or ""
                try:
                    if float(raw_score) < 0:
                        result.status = "not_run"
                        break
                except (ValueError, TypeError):
                    pass

                result.haplotype1 = hap1
                result.haplotype2 = hap2
                result.activity_score = _dash(raw_score)
                result.phenotype = _dash(row.get("phenotype") or row.get("Phenotype"))

                # Sub-alleles from candidate lists
                h1c = _dash(row.get("hap1_cand"))
                h2c = _dash(row.get("hap2_cand"))
                if h1c != "-" or h2c != "-":
                    result.sub_alleles = f"Hap1 candidates: {h1c} | Hap2 candidates: {h2c}"

                result.alternative_diplotypes = _dash(row.get("dip_cand"))

                sv1 = _dash(row.get("hap1_sv"))
                sv2 = _dash(row.get("hap2_sv"))
                result.sv_type = _dash(row.get("dip_sv")) or f"{sv1} / {sv2}"
                result.copy_number = "2" if "no_sv" in result.sv_type else "-"

                # Supporting variants from core fields
                variants: list[dict] = []
                effects: list[str] = []
                for hap_field in ("hap1_main_core", "hap2_main_core",
                                  "hap1_main_tag", "hap2_main_tag"):
                    v, e = _parse_stargazer_core(row.get(hap_field) or "")
                    variants.extend(v)
                    effects.extend(e)
                result.supporting_variants = variants
                result.functional_effects = "; ".join(dict.fromkeys(effects)) or "-"

                ssr = _dash(row.get("ssr"))
                h1s = _dash(row.get("hap1_score"))
                h2s = _dash(row.get("hap2_score"))
                result.allele_score = (
                    f"dip={raw_score}, hap1={h1s}, hap2={h2s}, SSR={ssr}"
                )

                af1g = _dash(row.get("hap1_af_mean_gene"))
                af2g = _dash(row.get("hap2_af_mean_gene"))
                af1m = _dash(row.get("hap1_af_mean_main"))
                af2m = _dash(row.get("hap2_af_mean_main"))
                result.mean_af = (
                    f"Hap1 gene={af1g}, main={af1m} | Hap2 gene={af2g}, main={af2m}"
                )

                result.status = "ok"
                break
    except Exception as exc:
        result.status = "failed"
        result.diplotype = f"parse error: {exc}"

    return result


# ── Parser: Aldy ──────────────────────────────────────────────────────────────
def parse_aldy(output_dir: str, gene: str, sample: str) -> CallerResult:
    result = CallerResult(
        tool="Aldy",
        phasing_method="Integer Linear Programming (joint CN + allele decomposition)",
    )

    candidates = [
        os.path.join(output_dir, "aldy", f"{gene}.aldy"),
        os.path.join(output_dir, "aldy", f"{gene.lower()}.aldy"),
        os.path.join(output_dir, "aldy", f"{sample}.{gene}.aldy"),
    ]
    candidates += glob.glob(
        os.path.join(output_dir, "aldy", "**", "*.aldy"), recursive=True
    )
    aldy_path = next((p for p in candidates if os.path.exists(p)), None)
    if aldy_path is None:
        return result

    try:
        with open(aldy_path) as f:
            lines = f.readlines()

        header: list[str] = []
        cpic_phenotype = "-"
        cpic_score = "-"
        diplotype = "-"
        minor = "-"
        variants: list[dict] = []
        effects: list[str] = []
        rsids: list[str] = []
        copies_seen: set[str] = set()

        for line in lines:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith("#"):
                # Extract cpic metadata from solution comment
                m_cpic = re.search(r"cpic=(\S+?)(?:;|$)", line)
                m_score = re.search(r"cpic_score=([\d.]+)", line)
                if m_cpic:
                    cpic_phenotype = m_cpic.group(1)
                if m_score:
                    cpic_score = m_score.group(1)
                # Parse column header
                candidate_header = [c.lstrip("#").strip() for c in line.split("\t")]
                if "Major" in candidate_header or "Sample" in candidate_header:
                    header = candidate_header
                continue

            if not header:
                continue

            parts = line.split("\t")
            if len(parts) < 4:   # need at least Sample, Gene, SolutionID, Major
                continue
            row = dict(zip(header, parts))

            # First data row: extract diplotype / minor
            if diplotype == "-":
                major = row.get("Major") or ""
                diplotype = major.strip() or "-"
                minor = _dash(row.get("Minor"))
                result.sv_type = "-"

            # Per-variant rows
            loc = row.get("Location") or row.get("Allele") or ""
            typ = row.get("Type") or ""
            cov = row.get("Coverage") or "-"
            eff = row.get("Effect") or "-"
            rsid = row.get("dbSNP") or "-"
            copy = row.get("Copy") or ""

            if loc and typ and re.search(r"[A-Z]>[A-Z]", typ):
                ref_alt = typ.split(">")
                variant = {
                    "allele": f"Copy {copy}" if copy != "" else "-",
                    "chrom":  "-",
                    "pos":    loc,
                    "ref":    ref_alt[0] if len(ref_alt) > 0 else "-",
                    "alt":    ref_alt[1] if len(ref_alt) > 1 else "-",
                    "af":     "-",
                    "depth":  cov,
                    "effect": eff if eff != "none" else "-",
                    "rsid":   rsid,
                }
                variants.append(variant)
                if eff and eff not in ("none", "-"):
                    effects.append(eff)
                if rsid and rsid != "-":
                    rsids.append(rsid)

        result.diplotype = diplotype
        result.haplotype1 = diplotype.split("/")[0].strip() if "/" in diplotype else diplotype
        result.haplotype2 = diplotype.split("/")[1].strip() if "/" in diplotype else "-"
        result.sub_alleles = minor
        result.phenotype = cpic_phenotype.replace("_", " ").title() if cpic_phenotype != "-" else "-"
        result.activity_score = cpic_score
        result.copy_number = "2"   # Aldy always solves to diploid by default
        result.supporting_variants = variants
        result.functional_effects = "; ".join(dict.fromkeys(effects)) or "-"
        result.dbsnp_ids = ", ".join(dict.fromkeys(rsids)) or "-"
        result.allele_score = f"SolutionID=1 (single best solution)"
        result.mean_af = f"{len(variants)} variants, depth in Coverage column"
        if diplotype != "-":
            result.status = "ok"

    except Exception as exc:
        result.status = "failed"
        result.diplotype = f"parse error: {exc}"

    return result


# ── Parser: StellarPGx ────────────────────────────────────────────────────────
def parse_stellarpgx(output_dir: str, gene: str, sample: str) -> CallerResult:
    result = CallerResult(
        tool="StellarPGx",
        phasing_method="Genome graph variant calling (graphtyper / vg)",
    )
    gene_lower = gene.lower()

    candidates = (
        glob.glob(os.path.join(output_dir, "stellarpgx", gene_lower, "alleles", "*.alleles"))
        + glob.glob(os.path.join(output_dir, "stellarpgx", gene, "alleles", "*.alleles"))
    )
    star_dir = os.path.join(output_dir, "stellarpgx", "star_allele_calls")
    if os.path.isdir(star_dir):
        candidates += (
            glob.glob(os.path.join(star_dir, "**", f"*{gene_lower}*"), recursive=True)
            + glob.glob(os.path.join(star_dir, "**", "*.txt"), recursive=True)
        )
    call_path = next((p for p in candidates if os.path.isfile(p)), None)
    if call_path is None:
        return result

    try:
        with open(call_path) as f:
            lines = f.readlines()

        variants: list[dict] = []
        cn = "-"

        for i, line in enumerate(lines):
            stripped = line.strip()

            if stripped.startswith("Initially computed CN"):
                m = re.search(r"=\s*(\d+)", stripped)
                if m:
                    cn = m.group(1)
                    result.copy_number = cn

            elif stripped == "Sample core variants:" and i + 1 < len(lines):
                var_str = lines[i + 1].strip()
                for tok in var_str.rstrip(";").split(";"):
                    tok = tok.strip()
                    if not tok:
                        continue
                    # format: pos~ref>alt~GT
                    m = re.match(r"(\d+)~([A-Z]+)>([A-Z]+)~([\d/|]+)", tok)
                    if m:
                        variants.append({
                            "allele": "-",
                            "chrom":  "-",
                            "pos":    m.group(1),
                            "ref":    m.group(2),
                            "alt":    m.group(3),
                            "af":     "-",
                            "depth":  "-",
                            "effect": "-",
                            "rsid":   "-",
                            "gt":     m.group(4),
                        })

            elif stripped == "Candidate alleles:" and i + 1 < len(lines):
                result.sub_alleles = lines[i + 1].strip()

            elif stripped == "Result:" and i + 1 < len(lines):
                dip = lines[i + 1].strip()
                if re.match(r"^\*", dip):
                    result.diplotype = dip
                    alleles = dip.split("/")
                    result.haplotype1 = alleles[0].strip() if alleles else "-"
                    result.haplotype2 = alleles[1].strip() if len(alleles) > 1 else "-"
                    result.status = "ok"

            elif stripped == "Activity score:" and i + 1 < len(lines):
                result.activity_score = _dash(lines[i + 1].strip())

            elif stripped == "Metaboliser status:" and i + 1 < len(lines):
                result.phenotype = _dash(lines[i + 1].strip())

            elif result.status != "ok":
                m = re.search(r"(\*\w+/\*\w+)", stripped)
                if m:
                    result.diplotype = m.group(1)
                    result.status = "ok"

        result.supporting_variants = variants
        result.sv_type = f"CN={cn}" if cn != "-" else "-"
        result.allele_score = "-"
        result.mean_af = "GT per variant (0/1 or 1/1)"

    except Exception as exc:
        result.status = "failed"
        result.diplotype = f"parse error: {exc}"

    return result


# ── Parser: OptiType (HLA-A, HLA-B) ─────────────────────────────────────────
def parse_optitype(output_dir: str, gene: str, sample: str) -> CallerResult:
    """Parse OptiType result TSV for HLA-A or HLA-B.

    OptiType outputs a TSV with columns:
        (index)  A1  A2  B1  B2  C1  C2  Reads  Objective
    Values are like "A*01:01" (without "HLA-" prefix).
    """
    result = CallerResult(
        tool="OptiType",
        phasing_method="ILP optimisation on HLA-reference-filtered MHC reads",
    )

    optitype_dir = os.path.join(output_dir, "optitype")
    # Look for <sample>_result.tsv or any *_result.tsv
    candidates = [os.path.join(optitype_dir, f"{sample}_result.tsv")]
    candidates += glob.glob(os.path.join(optitype_dir, "*_result.tsv"))
    tsv_path = next((p for p in candidates if os.path.exists(p)), None)
    if tsv_path is None:
        return result

    try:
        with open(tsv_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if gene == "HLA-A":
                    a1 = (row.get("A1") or "").strip()
                    a2 = (row.get("A2") or "").strip()
                    if a1 and a2:
                        # OptiType outputs "A*01:01"; prepend "HLA-" → "HLA-A*01:01"
                        result.haplotype1 = f"HLA-{a1}" if not a1.startswith("HLA") else a1
                        result.haplotype2 = f"HLA-{a2}" if not a2.startswith("HLA") else a2
                        result.diplotype  = f"{result.haplotype1}/{result.haplotype2}"
                        result.status = "ok"
                elif gene == "HLA-B":
                    b1 = (row.get("B1") or "").strip()
                    b2 = (row.get("B2") or "").strip()
                    if b1 and b2:
                        # OptiType outputs "B*57:01"; prepend "HLA-" → "HLA-B*57:01"
                        result.haplotype1 = f"HLA-{b1}" if not b1.startswith("HLA") else b1
                        result.haplotype2 = f"HLA-{b2}" if not b2.startswith("HLA") else b2
                        result.diplotype  = f"{result.haplotype1}/{result.haplotype2}"
                        result.status = "ok"
                reads = _dash(row.get("Reads"))
                obj   = _dash(row.get("Objective"))
                result.allele_score = f"Reads={reads}, Objective={obj}"
                break
    except Exception as exc:
        result.status = "failed"
        result.diplotype = f"parse error: {exc}"

    return result


# ── Parser: mutserve (MT-RNR1) ────────────────────────────────────────────────
def parse_mutserve(output_dir: str, gene: str, sample: str) -> CallerResult:
    """Parse pgx-mt.sh JSON result for MT-RNR1."""
    result = CallerResult(
        tool="mutserve",
        phasing_method="AF-based mitochondrial variant calling (mutserve v2)",
    )
    json_path = os.path.join(output_dir, "mt-rnr1", f"{sample}_mtrna1_result.json")
    if not os.path.exists(json_path):
        return result
    try:
        import json as _json
        data = _json.load(open(json_path))
        result.diplotype  = data.get("diplotype", "-")
        result.phenotype  = data.get("phenotype", "-")
        result.status     = "ok"
        variants = data.get("variants", [])
        if variants:
            af_strs = [f"{v['label']}(AF={v['af']:.2f},{v['type']})" for v in variants]
            result.allele_score = "; ".join(af_strs)
        else:
            result.allele_score = "No CPIC Level A variants detected"
    except Exception as exc:
        result.status    = "failed"
        result.diplotype = f"parse error: {exc}"
    return result


# ── SV mode note ──────────────────────────────────────────────────────────────
def _sv_note(gene: str) -> str:
    notes = []
    if gene in PYPGX_SV_GENES:
        notes.append("PyPGx: depth-of-coverage + VDR control stats")
    if gene in STARGAZER_SV_GENES:
        notes.append("Stargazer: GDF from BAM (paralog CN normalisation)")
    return ("SV mode — " + "; ".join(notes)) if notes else "no SVs expected"


# ── Output ────────────────────────────────────────────────────────────────────
def print_table(
    gene: str,
    sample: str,
    results: list[CallerResult],
    output_dir: str,
) -> None:
    W = 72
    SEP = "─" * W

    print()
    print("=" * W)
    print(f" PGx Star Allele Results")
    print(f" Gene:   {gene}")
    print(f" Sample: {sample}")
    print(f" Build:  GRCh38")
    print(f" SV:     {_sv_note(gene)}")
    print("=" * W)
    print(f"{'Tool':<14}{'Diplotype':<18}{'Activity Score':<18}{'Phenotype'}")
    print(SEP)

    called_diplotypes: list[str] = []
    for r in results:
        print(f"{r.tool:<14}{r.diplotype:<18}{r.activity_score:<18}{r.phenotype}")
        if r.status == "ok" and r.diplotype not in ("-", ""):
            called_diplotypes.append(r.diplotype)

    print(SEP)

    if called_diplotypes:
        counts = Counter(called_diplotypes)
        top_dip, top_count = counts.most_common(1)[0]
        total = len(results)
        print(f"Concordance: {top_count}/{total} tools agree on {top_dip}")
    else:
        print("Concordance: no calls available")

    print("=" * W)
    print()

    # ── TSV ───────────────────────────────────────────────────────────────────
    tsv_path = os.path.join(output_dir, f"{gene}_{sample}_comparison.tsv")
    with open(tsv_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            ["Gene", "Sample", "Build", "Tool", "Diplotype",
             "ActivityScore", "Phenotype", "Status", "SVMode"]
        )
        sv_mode = _sv_note(gene)
        for r in results:
            writer.writerow(
                [gene, sample, "GRCh38", r.tool, r.diplotype,
                 r.activity_score, r.phenotype, r.status, sv_mode]
            )
    print(f"Full results saved to: {tsv_path}")

    # ── Detail JSON ───────────────────────────────────────────────────────────
    json_path = os.path.join(output_dir, f"{gene}_{sample}_detail.json")
    detail = {
        "gene":    gene,
        "sample":  sample,
        "build":   "GRCh38",
        "sv_mode": _sv_note(gene),
        "tools":   {r.tool: r.to_dict(sample, gene) for r in results},
    }
    with open(json_path, "w") as fh:
        json.dump(detail, fh, indent=2)
    print(f"Detail JSON saved to:  {json_path}")


# ── Main ──────────────────────────────────────────────────────────────────────
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse PGx caller outputs and produce comparison table + detail JSON."
    )
    parser.add_argument("--gene",       required=True, help="Gene name (e.g. CYP2D6)")
    parser.add_argument("--sample",     required=True, help="Sample name")
    parser.add_argument("--output-dir", required=True, help="Results directory")
    args = parser.parse_args()

    gene = args.gene.upper()
    support = GENE_SUPPORT.get(gene, {})
    if not support:
        print(f"ERROR: Gene '{gene}' not found in support matrix.", file=sys.stderr)
        sys.exit(1)

    results: list[CallerResult] = []
    if support.get("pypgx"):
        results.append(parse_pypgx(args.output_dir, gene, args.sample))
    if support.get("stargazer"):
        results.append(parse_stargazer(args.output_dir, gene, args.sample))
    if support.get("aldy"):
        results.append(parse_aldy(args.output_dir, gene, args.sample))
    if support.get("stellarpgx"):
        results.append(parse_stellarpgx(args.output_dir, gene, args.sample))
    if support.get("optitype"):
        results.append(parse_optitype(args.output_dir, gene, args.sample))
    if support.get("mutserve"):
        results.append(parse_mutserve(args.output_dir, gene, args.sample))

    print_table(gene, args.sample, results, args.output_dir)


if __name__ == "__main__":
    main()
