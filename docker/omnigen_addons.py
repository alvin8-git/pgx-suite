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

ABO IS PROVISIONAL. The GRCh38 reference-allele identity at the ABO locus and
the exact variant coordinates require empirical confirmation (samtools faidx)
plus validation against a known-type control before any clinical use. Until then
every ABO result carries validation_status = "UNVALIDATED" and confidence "low".
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
    # Type 1 diabetes — DR3-DQ2 / DR4-DQ8 risk; DQB1*06:02 is protective
    t1d_risk = (_has("DRB1", "DRB1*03:01") and _has("DQB1", "DQB1*02")) or \
               (_has("DRB1", "DRB1*04") and _has("DQB1", "DQB1*03:02"))
    if t1d_risk:
        findings.append({
            "condition": "Type 1 diabetes susceptibility",
            "risk_haplotype": "DR3-DQ2 / DR4-DQ8",
            "tier": "informational",
            "note": ("HLA class-II T1D risk haplotype present. Susceptibility marker "
                     "only — most carriers never develop T1D."),
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
# 2. ABO blood type — SNP+indel assignment (VCF-Check-style)  ** PROVISIONAL **
# ─────────────────────────────────────────────────────────────────────────────
# Model (classic 3-site ABO*O.01 / A1 / B): O is defined by the rs8176719
# frameshift deletion; among non-O haplotypes A vs B is set by rs8176746 +
# rs8176747. This is diploid and unphased, so we assign genotype from zygosity
# counts. Reference-allele identity at these positions on GRCh38 is UNCONFIRMED,
# so the caller downgrades to UNVALIDATED (see module docstring).
ABO_VALIDATION_STATUS = (
    "UNVALIDATED — provisional ABO call. GRCh38 reference-allele identity and "
    "variant coordinates require empirical confirmation, and the caller must be "
    "validated against a known-type control before any clinical use. DO NOT USE "
    "CLINICALLY."
)


def abo_assign_type(sites: dict, defs: dict) -> dict:
    """Assign ABO genotype/serology from per-site observations.

    sites: {rsid: obs|None}, where obs = {"alt_copies": 0|1|2, "gt": str,
            "depth": int, "is_indel": bool}. None = position not observed / hom-ref.
    defs:  loaded abo_allele_defs.json (o_defining_rsid, ab_discriminating_rsids).

    Returns {ABO_type, alleles, confidence, per_site, validation_status}.
    Logic is validated by unit tests on synthetic, coordinate-independent inputs.
    """
    o_rs = defs["o_defining_rsid"]
    ab_rs = list(defs["ab_discriminating_rsids"])

    def _copies(rsid):
        obs = sites.get(rsid)
        return int(obs["alt_copies"]) if obs and obs.get("alt_copies") else 0

    o_copies = _copies(o_rs)                       # 0,1,2 deletion copies
    b_copies = max((_copies(r) for r in ab_rs), default=0)  # B-defining SNP copies

    non_o = 2 - o_copies
    # Clamp B copies to the number of non-O haplotypes available.
    b_copies = min(b_copies, non_o)
    a_copies = non_o - b_copies

    haps = ["O"] * o_copies + ["B"] * b_copies + ["A"] * a_copies
    haps.sort(key=lambda h: {"A": 0, "B": 1, "O": 2}[h])
    allele_names = {"A": "ABO*A1.01", "B": "ABO*B.01", "O": "ABO*O.01.01"}
    alleles = "/".join(allele_names[h] for h in haps)

    present = set(haps)
    if present == {"O"}:
        abo = "O"
    elif "A" in present and "B" in present:
        abo = "AB"
    elif "A" in present:
        abo = "A"
    elif "B" in present:
        abo = "B"
    else:
        abo = "Undetermined"

    # Confidence never exceeds "low" while the caller is UNVALIDATED.
    observed = sum(1 for r in [o_rs, *ab_rs] if sites.get(r) is not None)
    min_depth = min((sites[r]["depth"] for r in [o_rs, *ab_rs]
                     if sites.get(r) and sites[r].get("depth")), default=0)
    if abo == "Undetermined" or min_depth < 10:
        confidence = "low"
    else:
        confidence = "low"  # capped at low: UNVALIDATED
    per_site = []
    for r in [o_rs, *ab_rs]:
        obs = sites.get(r)
        per_site.append({
            "rsid": r,
            "alt_copies": obs["alt_copies"] if obs else 0,
            "gt": obs.get("gt", "0/0") if obs else "0/0",
            "depth": obs.get("depth") if obs else None,
        })
    return {
        "ABO_type": abo,
        "alleles": alleles,
        "confidence": confidence,
        "per_site": per_site,
        "validation_status": ABO_VALIDATION_STATUS,
        "n_sites_observed": observed,
    }


def write_abo_tsv(path: str, sample: str, result: dict) -> None:
    """Write results/<S>/abo/abo_type.tsv (stable OmniGen contract)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "ABO_type", "alleles", "confidence"])
        w.writerow([sample, result.get("ABO_type", "Undetermined"),
                    result.get("alleles", "-"), result.get("confidence", "low")])


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
