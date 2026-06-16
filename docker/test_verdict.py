#!/usr/bin/env python3
"""Self-check for the single-source verdict (pgx-compare.compute_verdict) and the
report's consumer (pgx-report.verdict_card). No framework.

Run: python3 docker/test_verdict.py
"""
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


cmp = _load("pgx_compare", "pgx-compare.py")
rpt = _load("pgx_report", "pgx-report.py")
CR = cmp.CallerResult


def _r(tool, dip, status="ok", pheno="Normal Metabolizer"):
    return CR(tool=tool, diplotype=dip, phenotype=pheno, status=status)


# ── compute_verdict (the authority) ───────────────────────────────────────────
def test_concordant_is_order_insensitive():
    v = cmp.compute_verdict([_r("PyPGx", "*1/*4"), _r("Aldy", "*4/*1")], no_call=False)
    assert v["status"] == "concordant" and v["n_agree"] == 2, v
    assert v["consensus_diplotype"] == "*1/*4", v


def test_two_vs_two_tie_is_discordant():
    v = cmp.compute_verdict(
        [_r("PyPGx", "*1/*4"), _r("Aldy", "*1/*4"),
         _r("Stargazer", "*1/*2"), _r("StellarPGx", "*1/*2")], no_call=False)
    assert v["status"] == "discordant", v  # the old code silently picked a winner


def test_majority_when_clear_plurality():
    v = cmp.compute_verdict(
        [_r("PyPGx", "*1/*4"), _r("Aldy", "*1/*4"), _r("Stargazer", "*1/*2")], no_call=False)
    assert v["status"] == "majority" and v["n_agree"] == 2 and v["n_called"] == 3, v


def test_low_coverage_is_no_call():
    v = cmp.compute_verdict([_r("PyPGx", "*1/*1")], no_call=True)
    assert v["status"] == "no_call" and v["consensus_diplotype"] == "NO_CALL", v


def test_nothing_called_is_no_data():
    v = cmp.compute_verdict([_r("PyPGx", "-", status="not_run", pheno="-")], no_call=False)
    assert v["status"] == "no_data", v


# ── verdict_card (the report consumer) ────────────────────────────────────────
def test_card_no_call_grey_no_phenotype():
    dip, ph, cc, _, _ = rpt.verdict_card(
        {"status": "no_call", "consensus_diplotype": "NO_CALL", "n_agree": 0, "n_called": 0})
    assert dip == "NO_CALL" and ph == "-" and cc == "card-no-data", (dip, ph, cc)


def test_card_discordant_red_no_phenotype():
    dip, ph, cc, _, _ = rpt.verdict_card({"status": "discordant", "n_agree": 2, "n_called": 4})
    assert dip == "DISCORDANT" and ph == "-" and cc == "card-red", (dip, ph, cc)


def test_card_concordant_green():
    dip, ph, cc, _, _ = rpt.verdict_card(
        {"status": "concordant", "consensus_diplotype": "*1/*4",
         "consensus_phenotype": "Normal Metabolizer", "n_agree": 4, "n_called": 4})
    assert dip == "*1/*4" and cc == "card-green", (dip, cc)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
    print("all verdict self-checks passed")
