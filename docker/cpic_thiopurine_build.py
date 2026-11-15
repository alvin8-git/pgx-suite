#!/usr/bin/env python3
"""Refresh docker/cpic_thiopurine.json (TPMT and NUDT15 diplotype -> CPIC phenotype) from the CPIC API.

    python3 docker/cpic_thiopurine_build.py [--retrieved YYYY-MM-DD] [--out docker/cpic_thiopurine.json]

Source (fetched live, per gene): https://api.cpicpgx.org/v1/diplotype?genesymbol=eq.<GENE>
                                 &select=diplotype,generesult  (CPIC thiopurine guideline lookup table)
No phenotype is written from memory: the file is that table, keyed by the order-insensitive diplotype.
PGX_CPIC_SNAPSHOT=<dir> replays <dir>/<GENE>_diplotype.json instead of the network (see cpic_cyp2d6_build.py).
"""
import argparse, datetime, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cpic_cyp2d6_build import _get  # noqa: E402

GENES = ("TPMT", "NUDT15")
URL = "https://api.cpicpgx.org/v1/diplotype?genesymbol=eq.{g}&select=diplotype,generesult&order=diplotype"


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
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--retrieved", default=datetime.date.today().isoformat())
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "cpic_thiopurine.json"))
    a = ap.parse_args()
    data = build(a.retrieved)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print("wrote", a.out)
