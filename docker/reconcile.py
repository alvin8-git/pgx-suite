#!/usr/bin/env python3
"""Variant-tier reconciliation (design A.1).

Star-allele callers disagree on *nomenclature* far more than on *biology*: PyPGx
calls DPYD `c.85T>C (*9A)`, Stargazer `*S10`, Aldy `*9A` — the same variant, three
names. The plain ensemble (pgx-compare.py) compares diplotype strings and marks
that "discordant". This module collapses known synonyms to a canonical token so
the verdict reflects genuine agreement.

Scoped deliberately: only well-established synonyms in `allele_synonyms.json`.
VKORC1 H-haplotypes and NUDT15 *5/*6 are NOT collapsed (they are not synonyms);
those go to array adjudication instead.

Backward-compatible: a gene with no synonym entry canonicalizes to the same
order-insensitive string the ensemble already used, so its verdict is unchanged.

Usage:
    reconcile.py --selftest
    reconcile.py --prove <SAMPLE> [--results results] [--gene DPYD]
"""
import argparse
import glob
import json
import os
import re
import sys

_HERE = os.path.dirname(os.path.realpath(__file__))
_SYN_PATH = os.path.join(_HERE, "allele_synonyms.json")

with open(_SYN_PATH, encoding="utf-8") as _fh:
    _RAW = json.load(_fh)

# Reference-allele synonyms (collapse to "*1") and per-gene token -> canonical.
_REFERENCE = {t.lower() for t in _RAW.get("_reference", [])}
SYNONYMS: dict[str, dict[str, str]] = {}
for _gene, _entry in _RAW.items():
    if _gene.startswith("_"):
        continue
    token_to_canon: dict[str, str] = {}
    for _canon, _aliases in _entry.items():
        token_to_canon[_canon.lower()] = _canon
        for _a in _aliases:
            token_to_canon[_a.lower()] = _canon
    SYNONYMS[_gene] = token_to_canon


def _strip_allele(tok: str) -> str:
    """Reduce one allele token to a comparable base form (lowercased)."""
    t = (tok or "").strip().lower()
    if not t:
        return ""
    t = re.sub(r"\s*\(.*?\)", "", t)   # drop "(*9a)" parentheticals
    t = re.sub(r"x\d+", "", t)          # copy number  *2xn / *106x2 -> *106
    t = t.split("+")[0]                  # *28+*80 -> *28 (base allele)
    t = t.split(";")[0]                  # pypgx "*17;*1;" -> *17
    return t.strip()


def canonicalize_allele(gene: str, tok: str) -> str:
    """Map one allele token to its canonical form for this gene.

    Copy-number / sub-allele stripping is applied ONLY for genes with a curated
    synonym map. For every other gene we keep the token as-is (lowercased) and
    apply only the universal Reference->*1 collapse — so CNV-meaningful calls
    like CYP2D6 `*1x2` or `*36+*10` are NOT flattened into `*1`/`*36`.
    """
    t = (tok or "").strip().lower()
    if not t or t in ("-", "."):
        return t
    if gene in SYNONYMS:
        base = _strip_allele(tok)
        if base in _REFERENCE:
            return "*1"
        return SYNONYMS[gene].get(base, base)
    return "*1" if t in _REFERENCE else t


def canonicalize_diplotype(gene: str, dip: str) -> str:
    """Order-insensitive, synonym-collapsed diplotype.

    For a gene with no synonym entry this matches pgx-compare's
    _normalize_diplotype (sorted halves), so verdicts stay unchanged.
    """
    if not dip or dip == "-":
        return "-"
    if "/" not in dip:
        return canonicalize_allele(gene, dip)
    halves = [canonicalize_allele(gene, p) for p in dip.split("/")]
    return "/".join(sorted(halves))


# ── proof / self-check ────────────────────────────────────────────────────────
def _prove(sample: str, results_dir: str, only_gene: str | None):
    """Show that callers reconcile to the same canonical diplotype on real data."""
    genes = [only_gene] if only_gene else ["DPYD"]
    any_seen = False
    for gene in genes:
        det = glob.glob(os.path.join(results_dir, sample, "genes", gene, "*detail.json"))
        if not det:
            continue
        any_seen = True
        d = json.load(open(det[0]))
        raw, canon = {}, set()
        for t, info in d.get("tools", {}).items():
            dip = info.get("diplotype")
            if dip in (None, "-", ""):
                continue
            c = canonicalize_diplotype(gene, dip)
            raw[t] = (dip, c)
            canon.add(c)
        print(f"[{sample}] {gene}: pipeline verdict = "
              f"{d.get('verdict', {}).get('status', '?')} "
              f"({d.get('verdict', {}).get('consensus_diplotype', '?')})")
        for t, (dip, c) in raw.items():
            print(f"    {t:11} {dip:32} -> {c}")
        verdict = "CONCORDANT after reconciliation" if len(canon) == 1 else \
                  f"still differs: {sorted(canon)}"
        print(f"    => {verdict}\n")
    if not any_seen:
        print(f"(no detail.json found for {sample})")


def _selftest():
    # DPYD: three caller spellings of the same variant collapse to one canonical
    assert canonicalize_diplotype("DPYD", "Reference/c.85T>C (*9A)") == \
           canonicalize_diplotype("DPYD", "*1/*S10") == \
           canonicalize_diplotype("DPYD", "*1/*9A"), "DPYD *9A synonyms must reconcile"
    assert canonicalize_diplotype("DPYD", "*1/*5") == \
           canonicalize_diplotype("DPYD", "*1/*S6") == \
           canonicalize_diplotype("DPYD", "*1/c.1627A>G"), "DPYD *5 synonyms must reconcile"
    # CYP2C19 sub-allele collapse
    assert canonicalize_diplotype("CYP2C19", "*2/*38") == \
           canonicalize_diplotype("CYP2C19", "*2/*2"), "CYP2C19 *38 is a *2 sub-allele"
    # TPMT *1S inflation
    assert canonicalize_diplotype("TPMT", "*1/*1S") == "*1/*1"
    # order-insensitive
    assert canonicalize_diplotype("DPYD", "*9A/*1") == canonicalize_diplotype("DPYD", "*1/*9A")
    # backward-compat: a gene with no entry behaves like plain normalize
    assert canonicalize_diplotype("CYP2D6", "*4/*2") == "*2/*4"
    # we DO NOT over-collapse: VKORC1/NUDT15 are not in the map, stay distinct
    assert canonicalize_diplotype("NUDT15", "*1/*5") != canonicalize_diplotype("NUDT15", "*1/*6")
    assert "VKORC1" not in SYNONYMS, "VKORC1 must be left to adjudication, not collapsed"
    print("reconcile self-test passed")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--prove", metavar="SAMPLE")
    ap.add_argument("--results", default="results")
    ap.add_argument("--gene")
    args = ap.parse_args()
    if args.selftest:
        _selftest()
    elif args.prove:
        _prove(args.prove, args.results, args.gene)
    else:
        ap.error("pass --selftest or --prove <SAMPLE>")


if __name__ == "__main__":
    main()
