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


def norm_dip(d, gene="CYP2D6"):
    """Comparison key (see reconcile.norm_diplotype); None when not comparable."""
    return reconcile.norm_diplotype(gene, d)


# Preference for which agreeing caller's spelling is cleaned and passed to PharmCAT.
_PREFER = ("Cyrius", "PyPGx", "Stargazer", "StellarPGx", "Aldy")


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


def _fmt_score(x):
    return "" if x is None else f"{x:.2f}".rstrip("0").rstrip(".") if x != int(x) else f"{x:.1f}"


def call_cyp2d6(out):
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
    info = []
    for r in t:
        sc = reconcile.activity_score(r["Diplotype"])
        own = reconcile.phenotype_class(r.get("Phenotype"))
        mapped = reconcile.score_phenotype(sc)
        info.append({"tool": r["Tool"], "diplotype": r["Diplotype"], "key": norm_dip(r["Diplotype"]),
                     "score": sc, "own": own,
                     "class": reconcile.phenotype_class(mapped) if mapped else own})
    res["callers"] = [{"tool": i["tool"], "diplotype": i["diplotype"], "score": i["score"],
                       "phenotype": i["class"]} for i in info]
    n = len(info)
    desc = "; ".join(f"{i['tool']}={i['diplotype']} (score {_fmt_score(i['score']) or 'unknown'}, "
                     f"{i['class'] or 'no phenotype'})" for i in info)
    # a: unanimous diplotype
    if len({i["key"] for i in info}) == 1:
        pick = min(t, key=lambda r: _PREFER.index(r["Tool"]) if r["Tool"] in _PREFER else len(_PREFER))
        res.update(level="a", source="all callers agree",
                   diplotype=reconcile.clean_diplotype("CYP2D6", pick["Diplotype"]),
                   score=reconcile.activity_score(pick["Diplotype"]))
        return res
    # b: every score known and equal
    scores = {i["score"] for i in info}
    if None not in scores and len(scores) == 1:
        res.update(level="b", source=f"activity score agreed {n}/{n}", score=scores.pop())
        return res
    # e: Cyrius (the read-depth copy-number caller) has a call with a CPIC activity score and at least one other
    # caller's call has the same score. Cyrius's cleaned diplotype is passed with that score; callers that
    # disagree, an explicit Indeterminate among them, are listed in `why` and never block this level.
    # Precedence: a, then b (every caller on the score: nothing to dissent, and the passed score-only call is
    # unchanged), then e, then c (e carries a diplotype and a score, c a phenotype only).
    cy = next((i for i in info if i["tool"] == "Cyrius" and i["score"] is not None), None)
    if cy is not None:
        same = [i for i in info if i is not cy and i["score"] == cy["score"]]
        if same:
            others = [i for i in info if i is not cy and i["score"] != cy["score"]]
            note = ("; dissenting: " + "; ".join(
                f"{i['tool']}={i['diplotype']} (score {_fmt_score(i['score']) or 'unknown'}, "
                f"{i['own'] or i['class'] or 'no phenotype'})" for i in others)) if others else ""
            res.update(level="e", source=f"Cyrius + {len(same)} other caller(s) on activity score "
                                         f"{_fmt_score(cy['score'])}",
                       diplotype=reconcile.clean_diplotype("CYP2D6", cy["diplotype"]), score=cy["score"],
                       why=("agreeing on the score: " + ", ".join(i["tool"] for i in same) + note))
            return res
    # c: one phenotype class; an explicit Indeterminate blocks it
    if any("Indeterminate" in (i["own"], i["class"]) for i in info):
        res["why"] = "a caller reports Indeterminate (blocks phenotype agreement): " + desc
        return res
    # Every caller's call must have a CPIC activity score: the phenotype is then CPIC's for that score, never a
    # caller's own label (a caller may print "Normal" over an allele CPIC lists with unknown function), and no
    # caller is left out of the comparison.
    if None in scores:
        res["why"] = "callers disagree and a call has no CPIC activity score: " + desc
        return res
    if len({i["class"] for i in info}) == 1 and info[0]["class"]:
        res.update(level="c", source=f"phenotype agreed {n}/{n}", phenotype=f"{info[0]['class']} Metabolizer")
        return res
    res["why"] = "callers disagree: " + desc
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
