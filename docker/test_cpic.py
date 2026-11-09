#!/usr/bin/env python3
"""Schema guard for cpic.json + the diplotype_check registry (F7).

cpic.json holds the static CPIC clinical data; pgx-report.py re-attaches the
dynamic diplotype predicates from DIPLOTYPE_CHECKS at load. This test fails if
the data file drifts out of shape or the registry stops covering the genes that
need a direct-diplotype check.
"""
import importlib.util
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# Every gene entry must carry these data keys (diplotype_check is attached in code).
REQUIRED_KEYS = {"desc", "drugs", "high_pheno", "moderate_pheno", "landing_notes"}

# Genes whose interpretation needs a direct diplotype predicate (not phenotype text).
EXPECTED_CHECK_GENES = {
    "DPYD", "VKORC1", "RYR1", "CACNA1S", "MT-RNR1", "CYP4F2",
    "ABCG2", "HLA-A", "HLA-B", "CYP3A4", "IFNL3",
}


def _load_report_module():
    spec = importlib.util.spec_from_file_location("pgxr", os.path.join(HERE, "pgx-report.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    # 1) cpic.json parses and round-trips
    with open(os.path.join(HERE, "cpic.json"), encoding="utf-8") as fh:
        data = json.load(fh)
    assert len(data) >= 31, f"expected >=31 genes in cpic.json, got {len(data)}"

    # 2) every gene has the required data keys, and NO callable leaked into JSON
    for gene, entry in data.items():
        missing = REQUIRED_KEYS - set(entry)
        assert not missing, f"{gene} missing keys in cpic.json: {missing}"
        assert "diplotype_check" not in entry, \
            f"{gene}: diplotype_check must live in code, not cpic.json"

    # 3) the assembled CPIC_DB has the registry attached for exactly the right genes
    m = _load_report_module()
    db = m.CPIC_DB
    assert set(db) == set(data), "CPIC_DB gene set diverged from cpic.json"
    with_check = {g for g, e in db.items() if e.get("diplotype_check")}
    assert with_check == EXPECTED_CHECK_GENES, (
        f"diplotype_check coverage drifted: "
        f"missing={EXPECTED_CHECK_GENES - with_check}, extra={with_check - EXPECTED_CHECK_GENES}"
    )

    # 4) spot-check predicate behaviour (regression anchors)
    assert db["DPYD"]["diplotype_check"]("*2/*4") is True
    assert db["DPYD"]["diplotype_check"]("*1/*1") is False
    assert db["VKORC1"]["diplotype_check"]("rs9923231") is True
    assert db["CYP4F2"]["diplotype_check"]("*3/*1") is True
    assert db["CYP4F2"]["diplotype_check"]("*1/*1") is False
    assert db["CYP2D6"].get("diplotype_check") is None  # phenotype-driven, no predicate

    # 5) NAT2: CPIC phenotype from CPIC allele functions, overriding the callers.
    spec = importlib.util.spec_from_file_location("pgxcompare", os.path.join(HERE, "pgx-compare.py"))
    comp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(comp)
    ph = comp.nat2_cpic_phenotype
    assert ph("*4/*4") == "Rapid Metabolizer"
    for d in ("*4/*5", "*4/*6", "*4/*7", "*6/*4"):
        assert ph(d) == "Intermediate Metabolizer", d
    for d in ("*5/*6", "*6/*6", "*7/*14"):
        assert ph(d) == "Poor Metabolizer", d
    for d in ("*4/*10", "*4/*60", "*4/*999", "*4/*5B", "garbage", "-"):
        assert ph(d) == "Indeterminate", d
    # the strings must land in the categories cpic.json's NAT2 keys expect
    assert m._pheno_cat("Poor Metabolizer") == "poor"
    assert m._pheno_cat("Intermediate Metabolizer") == "intermediate"
    assert m._pheno_cat("Rapid Metabolizer") == "rapid"
    assert m._pheno_cat("Indeterminate") not in db["NAT2"]["moderate_pheno"]
    CR = comp.CallerResult
    def _r(tool, dip, phe="-"):
        c = CR(tool=tool); c.status, c.diplotype, c.phenotype = "ok", dip, phe; return c
    v = comp.compute_verdict([_r("PyPGx", "*4/*6", "Indeterminate"),
                              _r("Stargazer", "*4/*6", "normal_function"),
                              _r("Aldy", "*4/*6")], no_call=False, gene="NAT2")
    assert v["consensus_phenotype"] == "Intermediate Metabolizer"

    print("all cpic.json schema + registry checks passed")


if __name__ == "__main__":
    try:
        main()
    except AssertionError as e:
        print("FAIL:", e, file=sys.stderr)
        sys.exit(1)
