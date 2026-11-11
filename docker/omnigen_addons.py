#!/usr/bin/env python3
"""
omnigen_addons.py — deterministic helpers for the OmniGen add-on features that
pgx-suite calls/types from the aligned GRCh38 BAM:

  1. HLA class II         (arcasHLA)                   -> genes/HLA-*/*_comparison.tsv
  2. ABO blood type       (SNP+indel, VCF-Check-style) -> abo/abo_type.tsv

(mtDNA + Y-chromosome haplogroups were moved to OmniGen — they classify directly
from the 30X WGS VCFs downstream, so pgx-suite no longer emits them.)

EVERY function here is pure and stdlib-only so it is unit-testable in-session
WITHOUT a BAM, a container, or a full pipeline run (see test_omnigen_addons.py).
The wet-lab / caller wiring lives in the shell scripts + Snakefile rules; this
module only parses tool output and assigns/serialises calls.

ABO is typed from BAM pileup evidence (abo_pileup.py) at the c.261delG O site and
four A/B SNPs; GRCh38 carries the O allele, so the non-O allele is the insertion.
An uncovered O site gives "Not typed", never a default type.
"""
from __future__ import annotations

import csv
import json
import os

# ─────────────────────────────────────────────────────────────────────────────
# 1. HLA class II — arcasHLA genotype.json parser
# ─────────────────────────────────────────────────────────────────────────────
# arcasHLA `genotype` emits <sample>.genotype.json:
#   {"A":["A*01:01:01","A*11:01:01"], "DQB1":["DQB1*06:02:01","DQB1*03:01:01"], ...}
# A homozygous locus may carry a single entry.
def _gene_to_arcas_key(gene: str) -> str:
    """HLA-DQB1 -> DQB1 ; DQB1 -> DQB1 (arcasHLA keys are unprefixed)."""
    return gene[4:] if gene.upper().startswith("HLA-") else gene


def parse_arcashla_genotype(genotype: dict, gene: str) -> tuple[str, str] | None:
    """Return (haplotype1, haplotype2) formatted 'HLA-DQB1*06:02:01', or None if
    the gene is absent. Homozygous single-entry loci duplicate the allele."""
    key = _gene_to_arcas_key(gene)
    alleles = genotype.get(key) or genotype.get(key.upper())
    if not alleles:
        return None
    def _fmt(a: str) -> str:
        a = a.strip()
        return a if a.upper().startswith("HLA-") else f"HLA-{a}"
    h1 = _fmt(alleles[0])
    h2 = _fmt(alleles[1]) if len(alleles) > 1 else h1
    return h1, h2


# HLA class-II disease-association risk logic (susceptibility, NOT diagnostic).
# Celiac DQ2.5 / DQ8 are DQA1+DQB1 *haplotype* combinations, so evaluated across
# both loci. Returns a list of finding dicts consumable by the report.
def hla2_disease_findings(calls: dict[str, str]) -> list[dict]:
    """calls: {"HLA-DQA1": diplo, "HLA-DQB1": diplo, "HLA-DRB1": diplo} (or bare keys).
    Each diplo is a '/'-joined string of alleles. Returns risk findings."""
    def _has(gene, *needles):
        d = calls.get(gene) or calls.get(f"HLA-{gene}") or ""
        return any(n in d for n in needles)

    findings: list[dict] = []
    # Celiac disease — HLA-DQ2.5 (DQA1*05 + DQB1*02:01) or HLA-DQ8 (DQA1*03 + DQB1*03:02)
    dq25 = _has("DQA1", "DQA1*05") and _has("DQB1", "DQB1*02:01", "DQB1*02:02", "DQB1*02")
    dq8 = _has("DQA1", "DQA1*03") and _has("DQB1", "DQB1*03:02")
    if dq25 or dq8:
        which = " + ".join([w for w, f in (("DQ2.5", dq25), ("DQ8", dq8)) if f])
        findings.append({
            "condition": "Celiac disease susceptibility",
            "risk_haplotype": which,
            "tier": "moderate",
            "note": ("Celiac-permissive HLA-DQ haplotype present. HLA is necessary "
                     "but not sufficient for celiac disease (present in ~30-40% of "
                     "the general population). Susceptibility marker only — not diagnostic."),
        })
    # Narcolepsy type 1 — DQB1*06:02
    if _has("DQB1", "DQB1*06:02"):
        findings.append({
            "condition": "Narcolepsy type 1 association",
            "risk_haplotype": "DQB1*06:02",
            "tier": "informational",
            "note": ("DQB1*06:02 present (~98% of narcolepsy-cataplexy cases carry it, "
                     "but also ~12-25% of the general population). Susceptibility marker only."),
        })
    # Type 1 diabetes — DR3-DQ2 and/or DR4-DQ8 risk arms; DQB1*06:02 is protective.
    # Report ONLY the arm(s) the sample actually carries so the label can never imply
    # a haplotype that is absent (e.g. a DR4-DQ8-only sample must not read as DR3-DQ2).
    dr3_dq2 = _has("DRB1", "DRB1*03:01") and _has("DQB1", "DQB1*02")
    dr4_dq8 = _has("DRB1", "DRB1*04") and _has("DQB1", "DQB1*03:02")
    t1d_risk = dr3_dq2 or dr4_dq8
    if t1d_risk:
        present = " + ".join(w for w, f in (("DR3-DQ2", dr3_dq2), ("DR4-DQ8", dr4_dq8)) if f)
        absent = " + ".join(w for w, f in (("DR3-DQ2", dr3_dq2), ("DR4-DQ8", dr4_dq8)) if not f)
        note = (f"HLA class-II T1D risk haplotype present: {present}. "
                "Susceptibility marker only — most carriers never develop T1D.")
        if absent:
            note += f" ({absent} not detected.)"
        findings.append({
            "condition": "Type 1 diabetes susceptibility",
            "risk_haplotype": f"{present} (present)",
            "tier": "informational",
            "note": note,
        })
    if _has("DQB1", "DQB1*06:02") and t1d_risk:
        findings.append({
            "condition": "Type 1 diabetes (protective allele co-present)",
            "risk_haplotype": "DQB1*06:02 (protective)",
            "tier": "informational",
            "note": "DQB1*06:02 is strongly protective against T1D and dominantly modifies risk.",
        })
    return findings


# ─────────────────────────────────────────────────────────────────────────────
# 2. ABO blood type — read-level assignment from BAM pileup evidence
# ─────────────────────────────────────────────────────────────────────────────
# GRCh38 carries the O haplotype at c.261delG (rs8176719): the reference LACKS the
# G, and the non-O allele is a one-base C insertion after chr9:133257521. So a
# sample whose reads all match the reference there has TWO O alleles; it is never
# A/A. Among non-O haplotypes, B is set by four SNPs whose reference base is the
# A-like base (c.526C>G, c.703G>A, c.796C>A, c.803G>C). The sites are read from the
# BAM (abo_pileup.py): an uncovered site yields "Not typed", never a default type.
ABO_MIN_DEPTH = 10          # reads required to genotype a site
ABO_HOM_ALT = 0.80          # alt fraction >= this -> 2 copies
ABO_HET_LOW = 0.20          # alt fraction in [LOW, HOM) -> 1 copy; < LOW -> 0
ABO_VALIDATION_STATUS = (
    "Research use. Read-level ABO call at the c.261delG O site and four A/B SNPs "
    "(GRCh38 reference = O haplotype). A1/A2 and O subtypes, cis-AB and weak/Bombay "
    "phenotypes are not resolved. Validation: docs/ToolsDocumentation.md 5e (GIAB HG001-HG007)."
)


def _site_copies(var: dict, obs: dict | None) -> dict:
    """Genotype one site from pileup counts -> {copies, depth, alt, af, status}."""
    if not obs or not obs.get("depth"):
        return {"rsid": var["rsid"], "status": "no_data", "copies": None,
                "depth": 0, "n_alt": 0, "n_ref": 0, "af": None}
    depth = int(obs["depth"])
    if var["type"] == "insertion":
        ins = obs.get("ins_seqs", {}).get(var["non_o_insertion"], 0)
        n_alt = ins
        n_ref = max(obs.get("ref_count", 0) + sum(obs.get("alt_counts", {}).values()) - ins, 0)
    else:
        n_alt = int(obs.get("alt_counts", {}).get(var["b_allele"], 0))
        n_ref = int(obs.get("ref_count", 0)) + int(obs.get("alt_counts", {}).get(var["ref_allele"], 0))
    informative = n_alt + n_ref
    if informative < ABO_MIN_DEPTH:
        return {"rsid": var["rsid"], "status": "no_data", "copies": None,
                "depth": depth, "n_alt": n_alt, "n_ref": n_ref, "af": None}
    af = n_alt / informative
    copies = 2 if af >= ABO_HOM_ALT else (1 if af >= ABO_HET_LOW else 0)
    return {"rsid": var["rsid"], "status": "ok", "copies": copies, "depth": depth,
            "n_alt": n_alt, "n_ref": n_ref, "af": round(af, 3)}


def abo_assign_type(pileup: dict, defs: dict) -> dict:
    """Assign ABO genotype/serology from per-site pileup evidence.

    pileup: {"sites": {rsid: {depth, ref, ref_count, alt_counts, ins_seqs, ...}}}
            (output of abo_pileup.py). A missing/empty pileup, or an uncovered
            O-defining site, gives ABO_type "Not typed" (ERROR is never NEGATIVE).

    Returns {ABO_type, alleles, confidence, per_site, validation_status,
             n_sites_observed, flags, evidence}.
    """
    psites = (pileup or {}).get("sites", {}) or {}
    var_by_rs = {v["rsid"]: v for v in defs["variants"]}
    o_rs = defs["o_defining_rsid"]
    ab_rs = list(defs["ab_discriminating_rsids"])
    extra = [r for r in (list(defs.get("supporting_rsids", [])) + [defs.get("o2_rsid")]) if r]
    calls = {r: _site_copies(var_by_rs[r], psites.get(r)) for r in [o_rs, *ab_rs, *extra]
             if r in var_by_rs}
    flags = []

    def _result(abo, alleles, conf, why=None):
        if why:
            flags.append(why)
        per_site = [{"rsid": r, "role": var_by_rs[r]["role"], "cdna": var_by_rs[r]["cdna"],
                     "pos": var_by_rs[r]["pos_grch38"], **c} for r, c in calls.items()]
        # legacy keys kept for existing readers
        for ps in per_site:
            ps["alt_copies"] = ps["copies"] if ps["copies"] is not None else 0
            ps["gt"] = {None: "./.", 0: "0/0", 1: "0/1", 2: "1/1"}[ps["copies"]]
        return {"ABO_type": abo, "alleles": alleles, "confidence": conf,
                "per_site": per_site, "validation_status": ABO_VALIDATION_STATUS,
                "n_sites_observed": sum(1 for c in calls.values() if c["status"] == "ok"),
                "flags": flags,
                "evidence": {"min_depth_required": ABO_MIN_DEPTH,
                             "o_site_depth": calls[o_rs]["depth"],
                             "o_site_reads_non_o": calls[o_rs]["n_alt"],
                             "o_site_reads_o": calls[o_rs]["n_ref"]}}

    o_call = calls[o_rs]
    if o_call["status"] != "ok":
        return _result("Not typed", "-", "none",
                       f"O-defining site c.261delG not covered (depth {o_call['depth']}, "
                       f"need >= {ABO_MIN_DEPTH}); no ABO call made")

    non_o = o_call["copies"]            # copies of the C insertion (non-O)
    o_copies = 2 - non_o
    if non_o == 0:                       # reference-homozygous at c.261 -> O/O
        haps = ["O", "O"]
    else:
        b_calls = [calls[r] for r in ab_rs]
        ok_b = [c for c in b_calls if c["status"] == "ok"]
        if len(ok_b) < 2:
            return _result("Not typed", "-", "none",
                           "non-O allele present but fewer than 2 of 4 A/B SNPs covered; "
                           "A vs B cannot be called")
        b_copies = min(max(c["copies"] for c in ok_b), non_o)
        # Sites should agree: every covered B SNP should carry the same copy number
        if len({c["copies"] for c in ok_b}) > 1:
            flags.append("A/B SNPs disagree (" +
                         ", ".join(f"{c['rsid']}={c['copies']}" for c in ok_b) +
                         "); possible recombinant/cis-AB/subtype, review")
        if len(ok_b) < len(b_calls):
            flags.append(f"{len(b_calls) - len(ok_b)} of 4 A/B SNPs not covered")
        a_copies = non_o - b_copies
        haps = ["O"] * o_copies + ["B"] * b_copies + ["A"] * a_copies
    o2 = calls.get(defs.get("o2_rsid"))
    if o2 and o2["status"] == "ok" and o2["copies"]:
        flags.append("c.802G>A (O2 marker) present: an O allele without c.261delG may be "
                     "present, so a non-O call may under-count O alleles")
    haps.sort(key=lambda h: {"A": 0, "B": 1, "O": 2}[h])
    alleles = "/".join(f"ABO*{h}" for h in haps)
    present = set(haps)
    if present == {"O"}:
        abo = "O"
    elif {"A", "B"} <= present:
        abo = "AB"
    elif "A" in present:
        abo = "A"
    else:
        abo = "B"
    ok_depths = [c["depth"] for c in calls.values() if c["status"] == "ok"
                 and c["rsid"] in [o_rs, *ab_rs]]
    conf = "high" if (not flags and min(ok_depths) >= 20) else (
        "medium" if min(ok_depths) >= ABO_MIN_DEPTH and len(flags) <= 1 else "low")
    return _result(abo, alleles, conf)


def write_abo_tsv(path: str, sample: str, result: dict) -> None:
    """Write results/<S>/abo/abo_type.tsv (stable OmniGen contract)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "ABO_type", "alleles", "confidence", "o_site_depth",
                    "n_sites_observed", "flags"])
        w.writerow([sample, result.get("ABO_type", "Not typed"),
                    result.get("alleles", "-"), result.get("confidence", "none"),
                    (result.get("evidence") or {}).get("o_site_depth", ""),
                    result.get("n_sites_observed", 0),
                    "; ".join(result.get("flags", []))])


def read_abo_tsv(path: str) -> dict | None:
    if not os.path.isfile(path):
        return None
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows[0] if rows else None


def load_abo_defs(path: str | None = None) -> dict:
    path = path or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "abo_allele_defs.json")
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)
