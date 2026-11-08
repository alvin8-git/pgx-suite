#!/usr/bin/env python3
"""Build a PharmCAT outside-call file from this pipeline's own results.

PharmCAT 3.2.0 does not type HLA-A/HLA-B from a VCF and `pharmcat_pipeline`
has no outside-call option; only the core tool does (`pharmcat.jar -po FILE`,
https://pharmcat.clinpgx.org/using/Outside-Call-Format/). This script writes:

  <OUT>/pharmcat/outside_calls.pharmcat.tsv  the file given to `-po` (gene<TAB>diplotype)
  <OUT>/pharmcat/outside_calls.tsv           record: gene, diplotype, source, passed, reason

Rules (never fabricate a call; an omitted gene is logged with the reason):
  HLA-A, HLA-B : OptiType diplotype from genes/<G>/<G>_*_comparison.tsv, cut to
                 two-field resolution, written as `*15:02/*55:01`.
  CYP2D6       : only when >=2 callers typed it and every caller that typed it
                 reports the same diplotype (alleles sorted, case/space ignored).
"""
import argparse, csv, glob, json, os, re, sys

HLA = ("HLA-A", "HLA-B")
UNTYPED = {"", "-", "NA", "N/A", "NONE", "UNKNOWN"}


def _newest(paths):
    return max(paths, key=os.path.getmtime) if paths else None


def read_comparison(out, gene):
    """Rows of the newest genes/<gene>/<gene>_*_comparison.tsv, or (None, [])."""
    f = _newest(glob.glob(os.path.join(out, "genes", gene, f"{gene}_*_comparison.tsv")))
    if not f:
        return None, []
    with open(f, newline="") as fh:
        return f, list(csv.DictReader(fh, delimiter="\t"))


def typed(row):
    d = (row.get("Diplotype") or "").strip()
    bad = d.upper() in UNTYPED or "NO_CALL" in d.upper() or "/" not in d
    return (row.get("Status") or "").strip().lower() == "ok" and not bad


def norm_dip(d):
    parts = [re.sub(r"\s+", "", p).lower() for p in d.split("/")]
    return "/".join(sorted(parts))


def hla_two_field(allele, gene):
    m = re.match(r"^(?:%s)?\*?(\d+):(\d+)" % re.escape(gene), allele.strip())
    return f"*{m.group(1)}:{m.group(2)}" if m else None


def call_hla(out, gene):
    f, rows = read_comparison(out, gene)
    if not f:
        return None, "", f"no {gene} comparison table"
    ot = [r for r in rows if (r.get("Tool") or "").strip().lower() == "optitype"]
    if not ot or not typed(ot[0]):
        return None, "OptiType", "OptiType did not type this gene"
    alleles = [hla_two_field(a, gene) for a in ot[0]["Diplotype"].split("/")]
    if len(alleles) != 2 or None in alleles:
        return None, "OptiType", f"unparseable diplotype {ot[0]['Diplotype']!r}"
    return "/".join(alleles), "OptiType", ""


def call_cyp2d6(out):
    f, rows = read_comparison(out, "CYP2D6")
    if not f:
        return None, "", "no CYP2D6 comparison table"
    t = [r for r in rows if typed(r)]
    names = ",".join(r["Tool"] for r in t)
    if len(t) < 2:
        return None, names, f"only {len(t)} caller(s) typed CYP2D6 (need >=2)"
    if len({norm_dip(r["Diplotype"]) for r in t}) != 1:
        return None, names, "callers disagree: " + "; ".join(
            f"{r['Tool']}={r['Diplotype']}" for r in t)
    return t[0]["Diplotype"].strip(), names, ""


def build(out):
    pc = os.path.join(out, "pharmcat")
    os.makedirs(pc, exist_ok=True)
    recs = []
    for g in HLA:
        d, src, why = call_hla(out, g)
        recs.append((g, d, src, why))
    d, src, why = call_cyp2d6(out)
    recs.append(("CYP2D6", d, src, why))
    with open(os.path.join(pc, "outside_calls.pharmcat.tsv"), "w") as fh:
        for g, d, _, _ in recs:
            if d:
                fh.write(f"{g}\t{d}\n")
    with open(os.path.join(pc, "outside_calls.tsv"), "w") as fh:
        fh.write("gene\tdiplotype\tsource_tool\tpassed_to_pharmcat\treason_omitted\n")
        for g, d, src, why in recs:
            fh.write(f"{g}\t{d or '-'}\t{src or '-'}\t{'yes' if d else 'no'}\t{why or '-'}\n")
            print(f"[outside-calls] {g}: " + (f"{d} (source {src})" if d else f"omitted - {why}"),
                  file=sys.stderr)
    return recs


def update_provenance(out, recs):
    p = os.path.join(out, "provenance.json")
    if not os.path.exists(p):
        return
    prov = json.load(open(p))
    prov["pharmcat_outside_calls"] = [
        {"gene": g, "diplotype": d, "source_tool": src, "passed": True} if d else
        {"gene": g, "passed": False, "reason": why} for g, d, src, why in recs]
    with open(p + ".tmp", "w") as fh:
        json.dump(prov, fh, indent=2)
    os.replace(p + ".tmp", p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="pipeline results dir (contains genes/, pharmcat/)")
    ap.add_argument("--update-provenance", action="store_true")
    a = ap.parse_args()
    recs = build(a.out)
    if a.update_provenance:
        update_provenance(a.out, recs)
    print(sum(1 for r in recs if r[1]))


if __name__ == "__main__":
    main()
