#!/usr/bin/env python3
"""Tests for the PharmCAT integration (sample-level CPIC reference caller).

parse_pharmcat reads OUT/pharmcat/<base>.report.json (two levels up from the
per-gene output dir) and pulls genes[<gene>].sourceDiplotypes[0].{label,phenotypes}.
PharmCAT is the verdict authority for UGT1A1/CYP2B6/CYP4F2. Pure-Python: PharmCAT
itself (Apptainer SIF) is exercised by the in-image verification, not here.
"""
import importlib.util
import json
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def _load_compare():
    spec = importlib.util.spec_from_file_location("c", os.path.join(HERE, "pgx-compare.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _write_report(out_root, genes):
    pcdir = os.path.join(out_root, "pharmcat")
    os.makedirs(pcdir, exist_ok=True)
    with open(os.path.join(pcdir, "pgx.report.json"), "w") as fh:
        json.dump({"genes": genes}, fh)


def main():
    c = _load_compare()

    # wiring: PharmCAT enabled + authoritative for the three target genes
    support, _, _ = c.load_gene_config()
    for g in ("UGT1A1", "CYP2B6", "CYP4F2"):
        assert support[g].get("pharmcat") is True, f"{g} should enable pharmcat"
        assert c.AUTHORITATIVE[g] == "PharmCAT", f"{g} should be PharmCAT-authoritative"

    with tempfile.TemporaryDirectory() as root:
        _write_report(root, {
            "UGT1A1": {"sourceDiplotypes": [
                {"label": "*1/*80", "phenotypes": ["Indeterminate"]},
                {"label": "*1/*80+*28", "phenotypes": ["Intermediate Metabolizer"]}]},
            "CYP2B6": {"sourceDiplotypes": [{"label": "*1/*4", "phenotypes": ["Rapid Metabolizer"]}]},
        })
        gene_dir = os.path.join(root, "genes", "UGT1A1")
        os.makedirs(gene_dir, exist_ok=True)
        r = c.parse_pharmcat(gene_dir, "UGT1A1", "S1")
        assert r.status == "ok" and r.diplotype == "*1/*80" and r.tool == "PharmCAT"
        assert r.phenotype == "Indeterminate"
        assert "*1/*80+*28" in (r.alternative_diplotypes or "")  # extra candidates kept

        # gene present in report but with no diplotype -> graceful no-call
        _write_report(root, {"UGT1A1": {"sourceDiplotypes": []}})
        assert c.parse_pharmcat(gene_dir, "UGT1A1", "S1").status == "failed"

    # missing report -> graceful failure
    with tempfile.TemporaryDirectory() as root:
        gd = os.path.join(root, "genes", "CYP2B6"); os.makedirs(gd)
        r = c.parse_pharmcat(gd, "CYP2B6", "S1")
        assert r.status == "failed" and "no PharmCAT report" in r.diplotype

    # authoritative verdict: PharmCAT wins over the star-allele vote
    CR = c.CallerResult
    def _r(t, dp, st="ok", phe="-"):
        x = CR(tool=t); x.status, x.diplotype, x.phenotype = st, dp, phe; return x
    rs = [_r("PyPGx", "*1/*28"), _r("Stargazer", "*1/*80"), _r("Aldy", "*1/*28_80"),
          _r("PharmCAT", "*1/*80", phe="Indeterminate")]
    v = c.compute_verdict(rs, no_call=False, gene="UGT1A1")
    assert v["status"] == "concordant" and v["consensus_diplotype"] == "*1/*80"
    assert v.get("authority") == "PharmCAT"

    # ── PharmCAT cross-check (informational; never feeds the verdict) ──────────
    with tempfile.TemporaryDirectory() as root:
        _write_report(root, {
            "CYP2C19": {"sourceDiplotypes": [{"label": "*1/*17", "phenotypes": ["Rapid Metabolizer"]}]},
        })
        gd = os.path.join(root, "genes", "CYP2C19"); os.makedirs(gd)
        # agrees: same normalised phenotype as the verdict
        xc = c.pharmcat_cross_check(gd, "CYP2C19", "S1", {"consensus_phenotype": "rapid_metabolizer"})
        assert xc and xc["agrees"] is True and xc["diplotype"] == "*1/*17", xc
        # differs: verdict phenotype disagrees -> surfaced as a review flag, not an error
        xc2 = c.pharmcat_cross_check(gd, "CYP2C19", "S1", {"consensus_phenotype": "Normal Metabolizer"})
        assert xc2 and xc2["agrees"] is False, xc2
        # authoritative gene -> None (it's the verdict there, not a cross-check)
        _write_report(root, {"UGT1A1": {"sourceDiplotypes": [{"label": "*1/*1", "phenotypes": ["Normal"]}]}})
        gdu = os.path.join(root, "genes", "UGT1A1"); os.makedirs(gdu)
        assert c.pharmcat_cross_check(gdu, "UGT1A1", "S1", {"consensus_phenotype": "Normal"}) is None
        # CYP2D6 -> None (PharmCAT does not call CYP2D6 structurally; Cyrius does)
        gdd = os.path.join(root, "genes", "CYP2D6"); os.makedirs(gdd)
        assert c.pharmcat_cross_check(gdd, "CYP2D6", "S1", {"consensus_phenotype": "x"}) is None
        # cross-check NEVER changes the verdict: a non-authoritative gene stays as voted
        rs2 = [_r("PyPGx", "*1/*17", phe="Rapid Metabolizer"),
               _r("Stargazer", "*1/*17", phe="Rapid Metabolizer")]
        v2 = c.compute_verdict(rs2, no_call=False, gene="CYP2C19")
        assert v2.get("authority") != "PharmCAT", v2   # PharmCAT is not CYP2C19's authority

    print("all PharmCAT integration checks passed")


if __name__ == "__main__":
    try:
        main()
    except AssertionError as e:
        print("FAIL:", e, file=sys.stderr)
        sys.exit(1)
