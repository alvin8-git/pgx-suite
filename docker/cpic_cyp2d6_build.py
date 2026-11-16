#!/usr/bin/env python3
"""Refresh docker/cpic_cyp2d6.json from the CPIC API (no allele fact is hard-coded in code).

    python3 docker/cpic_cyp2d6_build.py [--retrieved YYYY-MM-DD] [--out docker/cpic_cyp2d6.json]

Sources (both fetched live):
  alleles          https://api.cpicpgx.org/v1/allele?genesymbol=eq.CYP2D6&select=name,clinicalfunctionalstatus,activityvalue
  score->phenotype https://api.cpicpgx.org/v1/gene_result?genesymbol=eq.CYP2D6&select=result,activityscore
Alleles whose activityvalue is not a plain number (n/a, null) are stored with activity null and never
yield an activity score. A `≥number` value (`*1x≥3`: `≥3.0`) is stored as activity null plus `activity_min`
(a lower bound, never an exact score).
"""
import argparse, datetime, json, os, re, urllib.request

ALLELE_URL = ("https://api.cpicpgx.org/v1/allele?genesymbol=eq.CYP2D6"
              "&select=name,clinicalfunctionalstatus,activityvalue&order=name")
RESULT_URL = ("https://api.cpicpgx.org/v1/gene_result?genesymbol=eq.CYP2D6"
              "&select=result,activityscore&order=activityscore")


def _get(url, name=None):
    """Fetch a CPIC API URL; with PGX_CPIC_SNAPSHOT=<dir> read <dir>/<name>.json instead (an offline replay of a
    saved API response, used to rebuild and to verify the committed tables without network access)."""
    snap = os.environ.get("PGX_CPIC_SNAPSHOT")
    if snap and name:
        with open(os.path.join(snap, name + ".json"), encoding="utf-8") as fh:
            return json.load(fh)
    req = urllib.request.Request(url, headers={"User-Agent": "pgx-suite-cpic-refresh"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


def build(retrieved):
    alleles = {}
    for a in _get(ALLELE_URL, "CYP2D6_allele"):
        v = a.get("activityvalue")
        ok = isinstance(v, str) and re.fullmatch(r"\d+(\.\d+)?", v)
        entry = {"activity": float(v) if ok else None, "function": a.get("clinicalfunctionalstatus")}
        m = re.fullmatch(r"≥(\d+(\.\d+)?)", v) if isinstance(v, str) else None
        if m:
            entry["activity_min"] = float(m.group(1))
        alleles[a["name"]] = entry
    exact, at_least = {}, {}
    for r in _get(RESULT_URL, "CYP2D6_gene_result"):
        s, res = (r.get("activityscore") or "").strip(), r["result"]
        if re.fullmatch(r"\d+(\.\d+)?", s):
            exact[f"{float(s):.2f}"] = res
        elif s.startswith("≥"):
            at_least[f"{float(s[1:]):.2f}"] = res
    return {"gene": "CYP2D6", "retrieved": retrieved,
            "source": {"alleles": ALLELE_URL, "score_phenotype": RESULT_URL},
            "alleles": alleles,
            "score_phenotype": {"exact": exact, "at_least": at_least}}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--retrieved", default=datetime.date.today().isoformat())
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                  "cpic_cyp2d6.json"))
    a = ap.parse_args()
    data = build(a.retrieved)          # fetch first: a failed refresh never truncates the file
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print("wrote", a.out)
