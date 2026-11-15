#!/usr/bin/env python3
"""Build a PharmCAT outside-call file from this pipeline's own results.

PharmCAT 3.2.0 does not type HLA-A/HLA-B from a VCF and `pharmcat_pipeline`
has no outside-call option; only the core tool does (`pharmcat.jar -po FILE`,
https://pharmcat.clinpgx.org/using/Outside-Call-Format/). This script writes:

  <OUT>/pharmcat/outside_calls.pharmcat.tsv  the file given to `-po` (gene, diplotype, phenotype,
                                             activity score; empty columns allowed)
  <OUT>/pharmcat/outside_calls.tsv           record: gene, call, level, source, passed, reason, per caller

Rules (never fabricate a call; an omitted gene is logged with the reason):
  HLA-A, HLA-B : OptiType diplotype from genes/<G>/<G>_*_comparison.tsv, cut to
                 two-field resolution, written as `*15:02/*55:01`.
  CYP2D6       : by level of agreement among the callers that typed it (>=2 with a
                 comparable call; see reconcile.norm_diplotype for the notation rules).
                 Levels, first match wins:
                   a  all normalised diplotypes equal -> diplotype (+ CPIC activity score
                      when known)                                  "all callers agree"
                   b  diplotypes differ, every caller has a known CPIC activity score and
                      the scores are equal -> activity score only  "activity score agreed n/n"
                   c  otherwise, every caller has the same phenotype class (from the score
                      via the CPIC score table, else the caller's own phenotype) -> phenotype
                      only                                         "phenotype agreed n/n"
                   (b is tried before e, e before c; see the comment in call_cyp2d6)
                   e  Cyrius has a call with a known CPIC activity score and >=1 other caller has
                      the same score -> Cyrius's cleaned diplotype + that score. Dissenting callers
                      (including an explicit Indeterminate) do not block it and are listed in the
                      reason column. Needs a CPIC score for Cyrius's call: no Cyrius call, or one
                      with an unscored allele, never reaches e.
                   d  otherwise omitted, with the reason.
                 The levels live in reconcile.cyp2d6_agreement, shared with the verdict's CYP2D6 phenotype. Folded
                 tandem spellings (`*36+*36+*10` = `*36x2+*10`) and CPIC `x>=3` rows (a lower bound: usable in
                 level c as a phenotype class, never as an equal score) are scored.
                 A caller that reports Indeterminate does not agree with anything: it blocks
                 level c. A caller with neither score nor phenotype is ignored in level c
                 only while >=3 callers remain. Activity values come from docker/cpic_cyp2d6.json
                 (CPIC API); an unknown allele, `rsNNN` tag or unknown copy number gives no
                 score. Per-caller diplotype/score/phenotype are recorded.
                 NB: a phenotype-only call makes PharmCAT annotate every activity score of
                 that phenotype.
"""
import argparse, csv, glob, json, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reconcile  # noqa: E402  (single diplotype normaliser, shared with pgx-compare)

HLA = ("HLA-A", "HLA-B")


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
    return reconcile.is_typed_call(row.get("Status"), row.get("Diplotype"))


def norm_dip(d, gene="CYP2D6"):
    """Comparison key (see reconcile.norm_diplotype); None when not comparable."""
    return reconcile.norm_diplotype(gene, d)


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


_fmt_score = reconcile._fmt_score


def call_cyp2d6(out):
    return _no_lettered(_call_cyp2d6(out))


def _call_cyp2d6(out):
    """-> dict(diplotype, phenotype, score, level, source, why, callers)."""
    res = {"gene": "CYP2D6", "diplotype": None, "phenotype": None, "score": None, "level": "-",
           "source": "", "why": "", "callers": []}
    f, rows = read_comparison(out, "CYP2D6")
    if not f:
        res["why"] = "no CYP2D6 comparison table"
        return res
    t = [r for r in rows if typed(r)]
    keyed = [(r, norm_dip(r["Diplotype"])) for r in t]
    t = [r for r, k in keyed if k is not None]
    skipped = [r["Tool"] for r, k in keyed if k is None]
    res["source"] = ",".join(r["Tool"] for r in t)
    if len(t) < 2:
        res["why"] = (f"only {len(t)} caller(s) with a comparable CYP2D6 call (need >=2)"
                      + (f"; {','.join(skipped)} unphased, excluded" if skipped else ""))
        return res
    shared = reconcile.cyp2d6_agreement([{"tool": r["Tool"], "diplotype": r["Diplotype"], "phenotype": r.get("Phenotype")}
                                         for r in t])
    res.update({k: v for k, v in shared.items() if k != "source" or v})
    return res


def _no_lettered(res):
    """A diplotype with a legacy lettered allele name (`*4N`, `*13C`, `*14A`) is never passed to PharmCAT: it
    cannot import it, and the letter has no cited mapping to a current star number."""
    if res.get("diplotype") and reconcile.LEGACY_LETTERED.search(res["diplotype"]):
        res.update(level="-", source="", why=f"diplotype {res['diplotype']} contains a legacy lettered allele name "
                   "(no cited mapping to a current star allele): no outside call", diplotype=None,
                   phenotype=None, score=None)
    return res


def _passed(r):
    return bool(r["diplotype"] or r["phenotype"] or r["score"] is not None)


def build(out):
    pc = os.path.join(out, "pharmcat")
    os.makedirs(pc, exist_ok=True)
    recs = []
    for g in HLA:
        d, src, why = call_hla(out, g)
        recs.append({"gene": g, "diplotype": d, "phenotype": None, "score": None,
                     "level": "hla" if d else "-", "source": src, "why": why, "callers": []})
    recs.append(call_cyp2d6(out))
    with open(os.path.join(pc, "outside_calls.pharmcat.tsv"), "w") as fh:
        for r in recs:
            if _passed(r):
                cols = [r["gene"], r["diplotype"] or "", r["phenotype"] or "", _fmt_score(r["score"])]
                fh.write("\t".join(cols).rstrip("\t") + "\n")
    with open(os.path.join(pc, "outside_calls.tsv"), "w") as fh:
        fh.write("gene\tdiplotype\tactivity_score\tphenotype\tlevel\tsource\tpassed_to_pharmcat"
                 "\treason_omitted\tper_caller\n")
        for r in recs:
            per = "; ".join(f"{c['tool']}={c['diplotype']}|{_fmt_score(c['score']) or '-'}|"
                            f"{c['phenotype'] or '-'}" for c in r["callers"]) or "-"
            fh.write("\t".join([r["gene"], r["diplotype"] or "-", _fmt_score(r["score"]) or "-",
                                r["phenotype"] or "-", r["level"], r["source"] or "-",
                                "yes" if _passed(r) else "no", r["why"] or "-", per]) + "\n")
            what = " ".join(x for x in (r["diplotype"] or "",
                                        f"phenotype={r['phenotype']}" if r["phenotype"] else "",
                                        f"AS={_fmt_score(r['score'])}" if r["score"] is not None else "") if x)
            print(f"[outside-calls] {r['gene']}: " + (f"{what} (level {r['level']}, {r['source']})" + (f" [{r['why']}]" if r["why"] else "")
                  if _passed(r) else f"omitted - {r['why']}"), file=sys.stderr)
    return recs


def update_provenance(out, recs):
    p = os.path.join(out, "provenance.json")
    if not os.path.exists(p):
        return
    prov = json.load(open(p))
    prov["pharmcat_outside_calls"] = [
        {"gene": r["gene"], "diplotype": r["diplotype"], "phenotype": r["phenotype"],
         "activity_score": r["score"], "level": r["level"], "source_tool": r["source"],
         "per_caller": r["callers"], "passed": True, "note": r["why"] or None} if _passed(r) else
        {"gene": r["gene"], "passed": False, "reason": r["why"], "per_caller": r["callers"]}
        for r in recs]
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
    print(sum(1 for r in recs if _passed(r)))


if __name__ == "__main__":
    main()
