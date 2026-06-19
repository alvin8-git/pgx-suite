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
        [_r("PyPGx", "*1/*4", pheno="Normal Metabolizer"),
         _r("Aldy", "*1/*4", pheno="Normal Metabolizer"),
         _r("Stargazer", "*1/*2", pheno="Intermediate Metabolizer"),
         _r("StellarPGx", "*1/*2", pheno="Intermediate Metabolizer")], no_call=False)
    assert v["status"] == "discordant", v  # the old code silently picked a winner


def test_majority_when_clear_plurality():
    # dissenter has a DIFFERENT phenotype, so this stays a string-majority
    v = cmp.compute_verdict(
        [_r("PyPGx", "*1/*4", pheno="Normal Metabolizer"),
         _r("Aldy", "*1/*4", pheno="Normal Metabolizer"),
         _r("Stargazer", "*1/*2", pheno="Intermediate Metabolizer")], no_call=False)
    assert v["status"] == "majority" and v["n_agree"] == 2 and v["n_called"] == 3, v


# ── clinical phenotype tier ───────────────────────────────────────────────────
def test_phenotype_tier_rescues_string_discordance():
    # DPYD-style: 3 distinct diplotype STRINGS (Stargazer *S codes, Aldy rsIDs)
    # but the two phenotype-emitting callers agree; Aldy abstains (Indeterminate).
    v = cmp.compute_verdict(
        [_r("PyPGx", "c.85T>C (*9A)/c.496A>G", pheno="Normal Metabolizer"),
         _r("Stargazer", "*S10/*S20", pheno="normal_metabolizer"),
         _r("Aldy", "*9A/*rs2297595", pheno="Indeterminate")],
        no_call=False, gene="DPYD")
    assert v["status"] == "concordant" and v.get("phenotype_tier") \
        and v.get("authority") == "Phenotype", v


def test_phenotype_tier_not_applied_when_phenotypes_disagree():
    v = cmp.compute_verdict(
        [_r("PyPGx", "*1/*4", pheno="Normal Metabolizer"),
         _r("Aldy", "*1/*5", pheno="Intermediate Metabolizer"),
         _r("Stargazer", "*1/*6", pheno="Poor Metabolizer")],
        no_call=False, gene="CYP2C19")
    assert v["status"] == "discordant", v


def test_phenotype_tier_needs_two_agreeing_callers():
    # only ONE caller emits a usable phenotype; the rest abstain -> no rescue
    v = cmp.compute_verdict(
        [_r("PyPGx", "*1/*4", pheno="Normal Metabolizer"),
         _r("Aldy", "*1/*5", pheno="Indeterminate"),
         _r("Stargazer", "*1/*6", pheno="Indeterminate")],
        no_call=False, gene="CYP2C9")
    assert v["status"] == "discordant", v


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


# Colour follows the VERDICT STATUS, not the raw tool fraction. Regression guard:
# an authority/phenotype-resolved gene is often concordant at only 1/N callers
# (e.g. PharmCAT alone). It must still be GREEN, never red. (Earlier the card
# coloured by n_agree/n_called, so these read as discordant-red.)
def test_card_authority_concordant_is_green_at_one_of_four():
    dip, ph, cc, na, nc = rpt.verdict_card(
        {"status": "concordant", "consensus_diplotype": "*1/*80+*28",
         "consensus_phenotype": "Intermediate Metabolizer", "authority": "PharmCAT",
         "n_agree": 1, "n_called": 4})
    assert cc == "card-green", (cc, na, nc)          # NOT card-red
    assert dip == "*1/*80+*28" and na == 1 and nc == 4   # raw n/N still surfaced


def test_card_phenotype_tier_concordant_is_green():
    _, _, cc, _, _ = rpt.verdict_card(
        {"status": "concordant", "consensus_diplotype": "Reference/c.85T>C (*9A)",
         "consensus_phenotype": "Normal Metabolizer", "authority": "Phenotype",
         "phenotype_tier": True, "n_agree": 2, "n_called": 3})
    assert cc == "card-green", cc


def test_card_majority_is_amber():
    _, _, cc, _, _ = rpt.verdict_card(
        {"status": "majority", "consensus_diplotype": "*1/*2",
         "consensus_phenotype": "Normal Metabolizer", "n_agree": 2, "n_called": 3})
    assert cc == "card-amber", cc


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
    print("all verdict self-checks passed")
