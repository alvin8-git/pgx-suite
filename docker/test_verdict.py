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


def test_cyp2b6_cis_trans_is_phase_ambiguous_not_concordant():
    # PyPGx/Stargazer/Aldy *5/*6 vs StellarPGx/PharmCAT *1/*7: same defining variants (*7 = *5 + *6 in cis).
    # CPIC (api.cpicpgx.org, 2026-11-13): both diplotypes are Intermediate Metabolizer.
    im = "Intermediate Metabolizer"
    rs = [_r("PyPGx", "*5/*6", pheno=im), _r("Stargazer", "*5/*6", pheno=im), _r("Aldy", "*5/*6", pheno=im),
          _r("StellarPGx", "*1/*7", pheno=im), _r("PharmCAT", "*1/*7", pheno=im)]
    v = cmp.compute_verdict(rs, no_call=False, gene="CYP2B6")
    # both readings printed (majority reading first), no phase asserted; all five callers support the variant set
    assert v["phase_ambiguous"] is True and v["consensus_diplotype"] == "*5/*6 or *1/*7", v
    assert v["status"] == "concordant" and v["n_agree"] == 5 and v["n_called"] == 5, v
    assert set(map(tuple, v["alternative_readings"].values())) == {
        ("PyPGx", "Stargazer", "Aldy"), ("StellarPGx", "PharmCAT")}, v
    assert "phase ambiguity" in v["note"] and "*5/*6" in v["note"] and "*1/*7" in v["note"], v
    assert v["consensus_phenotype"] == im, v
    # one caller's own phenotype table differs (seen on a real sample): named, does not block the shared phenotype
    rs1 = [_r("PyPGx", "*5/*6", pheno=im), _r("Stargazer", "*5/*6", pheno="normal_metabolizer"), _r("Aldy", "*5/*6", pheno="Intermediate"),
           _r("StellarPGx", "*1/*7"), _r("PharmCAT", "*1/*7", pheno=im)]
    v = cmp.compute_verdict(rs1, no_call=False, gene="CYP2B6")
    assert v["consensus_phenotype"] == im and v["consensus_diplotype"] == "*5/*6 or *1/*7", v
    # readings with different phenotypes: phenotype not passed
    rs2 = [_r("PyPGx", "*5/*6", pheno="Poor Metabolizer"), _r("PharmCAT", "*1/*7", pheno=im)]
    v = cmp.compute_verdict(rs2, no_call=False, gene="CYP2B6")
    assert v["phase_ambiguous"] and v["consensus_phenotype"].startswith("Indeterminate") and v["status"] == "majority", v
    # a real difference (not a phase swap) is untouched; identical calls stay concordant
    v = cmp.compute_verdict([_r("PyPGx", "*1/*6"), _r("PharmCAT", "*1/*7")], no_call=False, gene="CYP2B6")
    assert v["status"] == "concordant" and "phase_ambiguous" not in v, v
    v = cmp.compute_verdict([_r("PyPGx", "*1/*7"), _r("PharmCAT", "*1/*7")], no_call=False, gene="CYP2B6")
    assert v["status"] == "concordant" and "phase_ambiguous" not in v, v


# ── CYP2D6: Cyrius reports no phenotype; derived from the CPIC activity score ──
def _d2d6(dip, others=None, cy_pheno="-", filt="-"):
    """Cyrius call plus, by default, one corroborating caller on the same diplotype."""
    cy = _r("Cyrius", dip, pheno=cy_pheno)
    cy.alternative_diplotypes = filt
    if others is None:
        others = [("PyPGx", dip, "-")]
    return cmp.compute_verdict([cy] + [_r(t, d, pheno=p) for t, d, p in others], no_call=False, gene="CYP2D6")


def test_cyp2d6_phenotype_derived_from_activity_score():
    for dip, score, ph in [("*2/*4", "1", "Intermediate Metabolizer"), ("*4/*4", "0", "Poor Metabolizer"),
                           ("*2/*2x2", "3", "Ultrarapid Metabolizer"),
                           ("*36+*10/*36+*10", "0.5", "Intermediate Metabolizer"),
                           ("*3/*68+*4", "0", "Poor Metabolizer"), ("*1/*1", "2", "Normal Metabolizer")]:
        v = _d2d6(dip)
        assert v["consensus_phenotype"] == ph, (dip, v)
        assert v["phenotype_source"] == f"CPIC activity score {score} from {dip} (same diplotype: PyPGx)", v
        assert (v["consensus_diplotype"], v["status"], v["n_agree"]) == (dip, "concordant", 2), v


def test_cyp2d6_fail_closed_cases_keep_blank_phenotype():
    for dip in ["*1/*2x3", "*1/*4.013+*4", "*36+*10/*36+*36", "*1/*36+*36+*10", "*1/*2;*1/*4",
                "Cyrius_not_assigned_to_haplotypes"]:
        v = _d2d6(dip)
        assert v["consensus_phenotype"] == "-" and "phenotype_source" not in v, (dip, v)
        assert v["note"].startswith("activity score not available for"), (dip, v)
    v = _d2d6("*1/*2x2", filt="LowQ_high_CN")
    assert v["consensus_phenotype"] == "-" and "phenotype_source" not in v, v
    assert "Cyrius filter LowQ_high_CN: phenotype not derived" in v["note"], v
    v = _d2d6("*1/*4", others=[("PyPGx", "*1/*1", "Normal Metabolizer"), ("Aldy", "*1/*1", "Normal")])
    assert v["consensus_phenotype"] == "-" and (v["n_agree"], v["status"]) == (1, "concordant"), v
    assert "no other caller's diplotype has the same activity score: phenotype not derived" in v["note"], v


def test_cyp2d6_corroboration_is_by_activity_score_not_spelling():
    v = _d2d6("*3/*68+*4", others=[(t, "*3/*4", "-") for t in ("PyPGx", "Stargazer")])
    assert v["consensus_phenotype"] == "Poor Metabolizer", v
    assert v["phenotype_source"] == "CPIC activity score 0 from *3/*68+*4 (same score: PyPGx *3/*4, Stargazer *3/*4)", v
    assert v["n_agree"] == 1 and v["status"] == "concordant", v      # counts untouched
    v = _d2d6("*1/*13", others=[("PyPGx", "*5/*39", "-"), ("StellarPGx", "*5/*39", "-")])
    assert v["consensus_phenotype"] == "Intermediate Metabolizer", v


def test_cyp2d6_lettered_legacy_allele_never_corroborates():
    for other in ("*13C/*39", "*4NX/*4", "*36X/*39", "*13c/*39", "*4n/*4", "*36x/*10"):
        v = _d2d6("*1/*13", others=[("Stargazer", other, "-")])
        assert v["consensus_phenotype"] == "-" and "phenotype_source" not in v, (other, v)
    v = _d2d6("*1/*13", others=[("Stargazer", "*13C/*39", "-"), ("PyPGx", "*5/*39", "-")])
    assert v["consensus_phenotype"] == "Intermediate Metabolizer" and "Stargazer" not in v["phenotype_source"], v
    # an uppercase X copy number is not a letter suffix
    v = _d2d6("*2/*2x2", others=[("Stargazer", "*2/*2X2", "-")])
    assert v["consensus_phenotype"] == "Ultrarapid Metabolizer", v


def test_legacy_lettered_pattern():
    rx = cmp._LEGACY_LETTERED
    for t in ("*13C/*39", "*4N/*4", "*36X+*10", "*13c/*39", "*4n/*4", "*36x/*10", "*4NX/*4"):
        assert rx.search(t), t
    for t in ("*2X2", "*10X2", "*1x2", "*68+*4", "*4.021", "*2x2/*1", "*1/*36+*10"):
        assert not rx.search(t), t


def test_cyp2d6_same_diplotype_rule_survives_lettered_spelling():
    v = _d2d6("*4/*4", others=[("Stargazer", "*4N/*4", "-")])
    assert v["consensus_phenotype"] == "Poor Metabolizer" and v["n_agree"] == 2, v
    assert v["phenotype_source"] == "CPIC activity score 0 from *4/*4 (same diplotype: Stargazer)", v
    v = _d2d6("*1/*36+*10", others=[("Stargazer", "*1/*36X+*10", "-")])
    assert v["consensus_phenotype"] == "Normal Metabolizer" and "same diplotype: Stargazer" in v["phenotype_source"], v
    # both rules apply: the same-diplotype wording wins
    v = _d2d6("*4/*4", others=[("PyPGx", "*4/*4", "-"), ("Aldy", "*4/*4.021", "-")])
    assert "(same diplotype: PyPGx" in v["phenotype_source"], v


def test_cyp2d6_authority_own_phenotype_unchanged():
    v = _d2d6("*2/*4", cy_pheno="Poor Metabolizer")
    assert v["consensus_phenotype"] == "Poor Metabolizer" and "phenotype_source" not in v, v


def test_cyp2d6_caller_class_difference_is_a_note_only():
    v = _d2d6("*2/*4", others=[("Aldy", "*2/*4", "Normal"), ("PyPGx", "*2/*4", "Intermediate Metabolizer"),
                               ("StellarPGx", "*2/*4", "Indeterminate")])
    assert v["consensus_phenotype"] == "Intermediate Metabolizer", v
    assert "Aldy Normal" in v["note"] and "PyPGx" not in v["note"], v
    assert (v["n_agree"], v["n_called"], v["status"]) == (4, 4, "concordant"), v


def test_other_authority_gene_with_blank_phenotype_stays_blank():
    rs = [_r("PharmCAT", "*1/*1", pheno="-"), _r("Stargazer", "*1/*1", pheno="normal_metabolizer")]
    v = cmp.compute_verdict(rs, no_call=False, gene="CYP4F2")
    assert v["consensus_phenotype"] == "-" and "phenotype_source" not in v, v


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
    print("all verdict self-checks passed")
