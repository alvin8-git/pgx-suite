#!/usr/bin/env python3
"""Refresh docker/cpic_diplotype.json (CPIC diplotype -> phenotype, star-allele genes) from the CPIC API.

    python3 docker/cpic_diplotype_build.py [--retrieved YYYY-MM-DD] [--out docker/cpic_diplotype.json]

Source (fetched live, per gene): https://api.cpicpgx.org/v1/diplotype?genesymbol=eq.<GENE>
                                 &select=diplotype,generesult  (CPIC diplotype lookup table of each guideline)
Used only to decide whether the several diplotype readings of an unresolved call give one CPIC phenotype or
different ones (pgx-compare.py, reconcile.cpic_diplotype_phenotype). No phenotype is written from memory: the
file is that table, keyed by the order-insensitive diplotype. TPMT/NUDT15 stay in cpic_thiopurine.json, CYP2D6
in cpic_cyp2d6.json (activity score). PGX_CPIC_SNAPSHOT=<dir> replays <dir>/<GENE>_diplotype.json instead of
the network (see cpic_cyp2d6_build.py).
"""
import argparse, datetime, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cpic_cyp2d6_build import _get  # noqa: E402

GENES = ("CYP2B6", "CYP2C19", "CYP2C9", "CYP3A5", "DPYD", "SLCO1B1", "UGT1A1")
URL = "https://api.cpicpgx.org/v1/diplotype?genesymbol=eq.{g}&select=diplotype,generesult&order=diplotype"


AL_URL = "https://api.cpicpgx.org/v1/allele?genesymbol=eq.{g}&select=name,clinicalfunctionalstatus&order=name"
FN_GENES = ("CYP2C19", "DPYD", "TPMT")      # genes with an allele_synonyms.json entry: test_verdict checks each pair's function


def key(dip):
    """Order-insensitive diplotype key (`*3A/*1` == `*1/*3A`)."""
    return "/".join(sorted(h.strip() for h in dip.split("/")))


def build(retrieved):
    out = {"retrieved": retrieved, "source": {g: URL.format(g=g) for g in GENES}}
    for g in GENES:
        table = {}
        for r in _get(URL.format(g=g), g + "_diplotype"):
            k = key(r["diplotype"])
            assert table.get(k, r["generesult"]) == r["generesult"], f"{g} {k} has two CPIC results"
            table[k] = r["generesult"]
        out[g] = table
    out["allele_function"] = {g: {a["name"]: a["clinicalfunctionalstatus"] for a in _get(AL_URL.format(g=g), g + "_allele")}
                              for g in FN_GENES}
    out["source"]["allele_function"] = AL_URL
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--retrieved", default=datetime.date.today().isoformat())
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "cpic_diplotype.json"))
    a = ap.parse_args()
    data = build(a.retrieved)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print("wrote", a.out)
