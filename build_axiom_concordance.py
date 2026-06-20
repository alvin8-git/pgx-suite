#!/usr/bin/env python3
"""Build an Axiom-array concordance TSV for a validation sample.

Host-side validation utility (NOT part of the container pipeline). Reads the
Axiom microarray truth from Axiom.xlsx and the pipeline's per-gene detail.json,
and emits <results>/<sample>/axiom_concordance.tsv with BOTH array software
columns (Axiom P6 and Axiom P9) where available.

The TSV is then embedded in the HTML report via:
    pgx-report.py ... --axiom <results>/<sample>/axiom_concordance.tsv

Only run this for validation samples — standard clinical samples get no panel.

Usage:
    build_axiom_concordance.py <SAMPLE> [--xlsx Axiom.xlsx] [--results results]
"""
import argparse
import csv
import glob
import json
import os
import re
import sys

# Genes the Axiom panel genotypes, in display order.
ARRAY_GENES = ["ABCG2", "CYP2B6", "CYP2C19", "CYP2C9", "CYP2D6", "CYP3A4",
               "CYP3A5", "CYP4F2", "DPYD", "NAT2", "NUDT15", "RYR1",
               "SLCO1B1", "TPMT", "UGT1A1", "VKORC1"]
# Array values that mean "no call" — row is dropped if the array gave us nothing.
NO_ARRAY = {"", "none", "indeterminate", "wt/wt", "-"}


def core_alleles(dip):
    """Multiset (sorted tuple) of base star alleles, or None if uncallable."""
    if not dip:
        return None
    d = dip.strip()
    if d.lower() in ("discordant", "./.", "-", "indeterminate", "none", ""):
        return None
    out = []
    for h in re.split(r"/", d):
        h = h.strip()
        h = re.sub(r"x\d+", "", h)          # copy number  *2xN -> *2
        h = re.sub(r"[()]", "", h)          # (*1/*10)xN parens
        base = h.split("+")[0].split(";")[0].split("(")[0].strip()
        if base:
            out.append(base)
    return tuple(sorted(out)) if out else None


def array_candidates(val):
    """Array cell may hold comma-separated ambiguous diplotypes; yield each."""
    if not val:
        return []
    return [c for c in (x.strip() for x in val.split(",")) if c]


def match_one(consensus, array_val):
    """Compare consensus diplotype to one array call -> MATCH/PARTIAL/MISMATCH/None."""
    c = core_alleles(consensus)
    if c is None:
        return None
    best = None
    for cand in array_candidates(array_val):
        a = core_alleles(cand)
        if a is None:
            continue
        if set(c) == set(a):
            return "MATCH"
        if set(c) & set(a):
            best = "PARTIAL"
        elif best is None:
            best = "MISMATCH"
    return best


def best_match(consensus, p6, p9):
    rank = {"MATCH": 3, "PARTIAL": 2, "MISMATCH": 1}
    cands = [m for m in (match_one(consensus, p6), match_one(consensus, p9)) if m]
    if not cands:
        return ""  # consensus uncallable (DISCORDANT) — caller applies a fallback
    return max(cands, key=lambda m: rank[m])


# Array reports these genes as nucleotide genotypes; ref allele per gene.
SNP_REF = {"ABCG2": "C", "VKORC1": "G"}


def snp_variant_count(gene, val):
    """# of non-ref nucleotides in an array genotype like 'A/C','C/C','A/A'."""
    ref = SNP_REF[gene]
    return sum(1 for a in re.split(r"/", val) if a.strip() and a.strip() != ref)


def pipeline_variant_count(consensus):
    """# of non-reference alleles in a star-allele consensus, or None if uncallable."""
    if core_alleles(consensus) is None:
        return None
    return sum(1 for a in re.split(r"/", consensus)
               if "ref" not in a.lower() and a.strip() not in ("*1", ""))


def tool_alleles(pertool):
    """Set of base star alleles seen across all tools (for discordant fallback)."""
    out = set()
    for chunk in pertool.split(";"):
        if ":" not in chunk:
            continue
        dip = chunk.split(":", 1)[1]
        c = core_alleles(dip.split(" or ")[0].split("[")[0])
        if c:
            out.update(a for a in c if "ref" not in a.lower() and a not in ("*1",))
    return out


def array_variant_alleles(p6, p9):
    """Non-ref base star alleles the array reported (union over P6/P9 candidates)."""
    out = set()
    for val in (p6, p9):
        for cand in array_candidates(val):
            c = core_alleles(cand)
            if c:
                out.update(a for a in c if "ref" not in a.lower() and a not in ("*1",))
    return out


def adjudicate(gene, consensus, status, p6, p9, match):
    """Credibility direction + plain-English reason when array and NGS disagree.

    The array is NOT ground truth (design A.2). When they differ, say which is
    more credible and why, by mechanism — array false-negative (panel gap),
    array false-positive vs NGS miss (confirm), single-SNP (trust array), or
    NGS-unresolved. MATCH short-circuits to concordant.
    """
    if match == "MATCH":
        return ("concordant", "array and pipeline agree")
    if gene in SNP_REF:
        return ("array_authoritative",
                "single-SNP locus; array genotypes this directly, NGS adds notation noise")
    arr_variant = bool(array_variant_alleles(p6, p9))
    ngs_vc = pipeline_variant_count(consensus)
    # array reference, NGS callers agree on a variant -> array panel gap
    if not arr_variant and status in ("concordant", "majority") and ngs_vc:
        return ("array_FN",
                "array panel gap: NGS callers concordant on a variant the array does not score (NGS authoritative)")
    # array variant, NGS called reference -> array FP or NGS miss
    if arr_variant and ngs_vc in (0, None) and status != "discordant":
        return ("array_FP_or_NGS_FN",
                "array reports a variant NGS did not call; orthogonal confirmation (pileup/Sanger) needed")
    if status == "discordant":
        return ("ngs_unresolved",
                "NGS callers disagree even after reconciliation; manual review")
    return ("review", "nomenclature or array-panel difference; review")


def classify(gene, consensus, p6, p9, pertool):
    """Final agreement verdict, with SNP and discordant-consensus handling."""
    if gene in SNP_REF:
        av = next((snp_variant_count(gene, v) for v in (p6, p9)
                   if v and v.lower() not in NO_ARRAY), None)
        pv = pipeline_variant_count(consensus)
        if av is None or pv is None:       # array/pipeline uncallable -> flag
            return "MISMATCH"
        if av == pv:
            return "MATCH"
        return "PARTIAL" if (av > 0 and pv > 0) else "MISMATCH"

    m = best_match(consensus, p6, p9)
    if m:
        return m
    # consensus is DISCORDANT: fall back to per-tool vs array key alleles
    shared = array_variant_alleles(p6, p9) & tool_alleles(pertool)
    return "PARTIAL" if shared else "MISMATCH"


def load_axiom(xlsx, sample):
    import openpyxl
    wb = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    ws = wb["Sheet1"]
    rows = list(ws.iter_rows(values_only=True))
    hdr = rows[0]
    gi = {g: hdr.index(g) for g in ARRAY_GENES if g in hdr}
    out = {}  # software -> {gene: value}
    for r in rows[1:]:
        if r[0] != sample:
            continue
        sw = (r[2] or "").strip()  # "Axiom P6" / "Axiom P9"
        out[sw] = {g: ("" if r[gi[g]] is None else str(r[gi[g]]).strip())
                   for g in gi}
    return out


def load_pipeline(results, sample):
    out = {}
    for g in ARRAY_GENES:
        det = glob.glob(os.path.join(results, sample, "genes", g, "*detail.json"))
        if not det:
            continue
        d = json.load(open(det[0]))
        v = d.get("verdict", {})
        pertool = "; ".join(
            f"{t}:{info.get('diplotype', '-')}" for t, info in d.get("tools", {}).items()
            if info.get("diplotype") not in (None, "")
        )
        out[g] = {
            "consensus": v.get("consensus_diplotype", ""),
            "status": v.get("status", ""),
            "pertool": pertool,
        }
    return out


def _selftest():
    # star alleles: order-insensitive exact match
    assert classify("CYP2C19", "*17/*2", "*2/*17", "*2/*17", "") == "MATCH"
    # SNP genotype het variant agrees with pipeline het variant call
    assert classify("ABCG2", "*ref/*rs2231142", "A/C", "A/C", "") == "MATCH"
    assert classify("ABCG2", "*ref/*ref", "C/C", "C/C", "") == "MATCH"
    # array normal but pipeline (discordant) found an extra variant -> Review
    assert classify("DPYD", "DISCORDANT", "*1/*1", "*1/*1",
                    "PyPGx:Reference/*9A; Aldy:*1/*9A") == "MISMATCH"
    # discordant consensus but a tool shares the array's key allele -> Partial
    assert classify("UGT1A1", "DISCORDANT", "*1/*28", "*28+60+80+93/*60",
                    "PyPGx:*1/*80; Stargazer:*1/*28") == "PARTIAL"
    # partial overlap on callable consensus
    assert classify("CYP2D6", "*1x2/*36+*10", "Indeterminate", "(*1/*10)xN", "") == "PARTIAL"
    # homozygous-variant SNP genotype, pipeline uncallable -> Review
    assert classify("VKORC1", "DISCORDANT", "A/A", "A/A", "") == "MISMATCH"
    # adjudication directions
    assert adjudicate("DPYD", "Reference/c.85T>C (*9A)", "concordant", "*1/*1", "*1/*1",
                      "MISMATCH")[0] == "array_FN"          # NGS authoritative, array panel gap
    assert adjudicate("VKORC1", "DISCORDANT", "discordant", "A/A", "A/A",
                      "MISMATCH")[0] == "array_authoritative"
    assert adjudicate("RYR1", "DISCORDANT", "discordant", "", "c.6178G>T/WT",
                      "MISMATCH")[0] in ("array_FP_or_NGS_FN", "ngs_unresolved")
    assert adjudicate("CYP2C19", "*1/*17", "concordant", "*1/*17", "*1/*17",
                      "MATCH")[0] == "concordant"
    print("build_axiom_concordance self-test passed")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sample", nargs="?")
    ap.add_argument("--xlsx", default="Axiom.xlsx")
    ap.add_argument("--results", default="results")
    ap.add_argument("--axiom-sample", default=None,
                    help="Sample key to look up in Axiom.xlsx when it differs from "
                         "the results dir name (e.g. ILMN 'H-342' dir vs Axiom 'H_342')")
    ap.add_argument("--selftest", action="store_true", help="run assertion checks and exit")
    args = ap.parse_args()

    if args.selftest:
        _selftest()
        return
    if not args.sample:
        ap.error("sample is required (or pass --selftest)")

    axiom = load_axiom(args.xlsx, args.axiom_sample or args.sample)
    if not axiom:
        sys.exit(f"No Axiom rows for sample {args.axiom_sample or args.sample!r} in {args.xlsx}")
    p6 = axiom.get("Axiom P6", {})
    p9 = axiom.get("Axiom P9", {})
    pipe = load_pipeline(args.results, args.sample)

    out_path = os.path.join(args.results, args.sample, "axiom_concordance.tsv")
    n = 0
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Gene", "Axiom_P6", "Axiom_P9", "Consensus", "Status",
                    "Match", "Direction", "Reason", "PerTool"])
        for g in ARRAY_GENES:
            v6, v9 = p6.get(g, ""), p9.get(g, "")
            # drop rows where neither array software called the gene
            if v6.lower() in NO_ARRAY and v9.lower() in NO_ARRAY:
                continue
            pp = pipe.get(g)
            if pp is None:
                continue  # pipeline didn't assess this gene
            consensus = pp["consensus"]
            match = classify(g, consensus, v6, v9, pp["pertool"])
            direction, reason = adjudicate(g, consensus, pp["status"], v6, v9, match)
            w.writerow([g, v6, v9, consensus, pp["status"], match,
                        direction, reason, pp["pertool"]])
            n += 1
    print(f"wrote {out_path}  ({n} genes)")


if __name__ == "__main__":
    main()
