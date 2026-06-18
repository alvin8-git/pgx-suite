#!/usr/bin/env python3
"""Tests for variant-tier reconciliation (A.1) + the compute_verdict integration.

Guards two things:
1. reconcile.py canonicalization collapses known synonyms but does NOT over-collapse
   (VKORC1 H-haplotypes / NUDT15 *5 vs *6 stay distinct).
2. compute_verdict in pgx-compare.py counts on the canonical form, so the DPYD
   *9A/*S10/c.85T>C three-way nomenclature split now resolves to concordant.
"""
import importlib.util
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reconcile  # noqa: E402


def _load(modname, filename):
    spec = importlib.util.spec_from_file_location(modname, os.path.join(HERE, filename))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    # 1) reconcile.py's own assertions (canonicalize + no over-collapse)
    reconcile._selftest()

    # 2) synonym data schema
    with open(os.path.join(HERE, "allele_synonyms.json"), encoding="utf-8") as fh:
        raw = json.load(fh)
    assert "_reference" in raw, "reference-synonym list missing"
    assert "VKORC1" not in raw, "VKORC1 must NOT be in the synonym map (adjudication-only)"
    assert "NUDT15" not in raw, "NUDT15 *5/*6 are distinct variants, not synonyms"
    assert "DPYD" in raw and "CYP2C19" in raw

    # 3) compute_verdict integration: the DPYD three-name split must reconcile
    comp = _load("pgxcompare", "pgx-compare.py")
    CR = comp.CallerResult

    def r(dip, pheno="Normal Metabolizer"):
        c = CR(tool="t")
        c.status, c.diplotype, c.phenotype = "ok", dip, pheno
        return c

    dpyd = [r("Reference/c.85T>C (*9A)"), r("*1/*S10"), r("*1/*9A")]
    v = comp.compute_verdict(dpyd, no_call=False, gene="DPYD")
    assert v["status"] == "concordant", f"DPYD should reconcile to concordant, got {v['status']}"
    assert v["n_agree"] == 3 and v["reconciled"] is True

    # without a gene (legacy path) the same three strings stay split
    v_legacy = comp.compute_verdict(dpyd, no_call=False, gene=None)
    assert v_legacy["status"] != "concordant", "legacy path must not silently reconcile"

    # a genuinely discordant gene stays discordant after reconciliation
    disc = [r("*1/*2"), r("*1/*3"), r("*2/*2")]
    assert comp.compute_verdict(disc, no_call=False, gene="CYP2D6")["status"] in ("discordant", "majority")

    print("all reconcile + verdict-integration checks passed")


if __name__ == "__main__":
    try:
        main()
    except AssertionError as e:
        print("FAIL:", e, file=sys.stderr)
        sys.exit(1)
