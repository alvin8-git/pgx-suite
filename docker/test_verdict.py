#!/usr/bin/env python3
"""Self-check for the single-source verdict (pgx-compare.compute_verdict) and the
report's consumer (pgx-report.verdict_card). No framework.

Run: python3 docker/test_verdict.py
"""
import importlib.util
import os
import re

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
    # Stargazer's *S20 and Aldy's *rs2297595 are CPIC c.496A>G through dpyd_notation_map.json: all readings score Normal
    assert v["status"] == "concordant" and v["consensus_phenotype"] == "Normal Metabolizer" and v["authority"] == "Phenotype", v
    # a composite Aldy haplotype (*9A + extra rs ids) is no CPIC allele: unscorable, so (policy 2026-11-16) not counted; the other readings keep Normal
    v = cmp.compute_verdict(
        [_r("PyPGx", "c.85T>C (*9A)/c.496A>G", pheno="Normal Metabolizer"),
         _r("Stargazer", "*S10/*S20", pheno="normal_metabolizer"),
         _r("Aldy", "*9A/*9A+rs1801159+rs2297595", pheno="Indeterminate")],
        no_call=False, gene="DPYD")
    assert v["status"] == "concordant" and v["consensus_phenotype"] == "Normal Metabolizer" and "phenotype_split" not in v, v
    assert "reading *9A/*9A+rs1801159+rs2297595 has no CPIC phenotype and was not counted" in v["note"], v
    # readings CPIC scores alike keep the phenotype tier
    v = cmp.compute_verdict(
        [_r("PyPGx", "c.85T>C (*9A)/*1", pheno="Normal Metabolizer"),
         _r("Stargazer", "*S6/*1", pheno="normal_metabolizer"),
         _r("Aldy", "*1/*1", pheno="Normal Metabolizer")],
        no_call=False, gene="DPYD")
    assert v["status"] == "concordant" and v.get("phenotype_tier") and v["consensus_phenotype"] == "Normal Metabolizer", v


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
    # CPIC scores *1/*7 and *5/*6 alike (Intermediate): CPIC decides, the callers' Indeterminate label goes in the note
    assert v["phase_ambiguous"] and v["consensus_phenotype"] == "Intermediate Metabolizer" and v["status"] == "majority", v
    assert "callers said Indeterminate" in v["note"], v
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


def test_cyp2d6_tandem_and_x3_spellings_are_scored():
    # CPIC keys: `*36x2+*10` 0.25, `*36x2` 0.0, `*1x≥3`/`*2x≥3` ≥3.0 (cpic_cyp2d6.json)
    for dip, ph in [("*1/*36+*36+*10", "Normal Metabolizer"), ("*36+*10/*36+*36", "Intermediate Metabolizer"),
                    ("*1/*4.013+*4", "Intermediate Metabolizer")]:
        v = _d2d6(dip)
        assert v["consensus_phenotype"] == ph and "phenotype_source" in v, (dip, v)
    v = _d2d6("*1/*2x3")                      # >= 4.0 (lower bound): Ultrarapid; only the same diplotype corroborates
    assert v["consensus_phenotype"] == "Ultrarapid Metabolizer", v
    assert v["phenotype_source"] == "CPIC activity score >= 4 from *1/*2x3 (same diplotype: PyPGx)", v
    v = _d2d6("*1/*2x3", others=[("PyPGx", "*1/*2x4", "-")])    # same lower bound is not the same score
    assert v["consensus_phenotype"] == "-", v


def test_cyp2d6_fail_closed_cases_keep_blank_phenotype():
    for dip in ["*1/*2;*1/*4", "Cyrius_not_assigned_to_haplotypes"]:
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
    for other in ("*13C/*39", "*4NX/*4", "*13c/*39", "*4n/*4", "*14A/*1"):
        v = _d2d6("*1/*13", others=[("Stargazer", other, "-")])
        assert v["consensus_phenotype"] == "-" and "phenotype_source" not in v, (other, v)
    v = _d2d6("*1/*13", others=[("Stargazer", "*13C/*39", "-"), ("PyPGx", "*5/*39", "-")])
    assert v["consensus_phenotype"] == "Intermediate Metabolizer" and "Stargazer" not in v["phenotype_source"], v
    # an uppercase X copy number is not a letter suffix
    v = _d2d6("*2/*2x2", others=[("Stargazer", "*2/*2X2", "-")])
    assert v["consensus_phenotype"] == "Ultrarapid Metabolizer", v


def test_legacy_lettered_pattern():
    rx = cmp._LEGACY_LETTERED
    for t in ("*13C/*39", "*4N/*4", "*13c/*39", "*4n/*4", "*4NX/*4", "*14A/*1", "*14b/*1"):
        assert rx.search(t), t
    for t in ("*2X2", "*10X2", "*1x2", "*68+*4", "*4.021", "*2x2/*1", "*1/*36+*10", "*36X+*10", "*36x/*10", "*36X/*39"):
        assert not rx.search(t), t


def test_cyp2d6_same_diplotype_rule_survives_lettered_spelling():
    v = _d2d6("*4/*4", others=[("Stargazer", "*4N/*4", "-")])    # a lettered name is not *4: no corroboration
    assert v["consensus_phenotype"] == "-" and v["n_agree"] == 1, v
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


# ── PharmCAT "Unknown" is no call, not an authority answer ────────────────────
def _pc(dip, pheno="-"):
    return _r("PharmCAT", dip, pheno=pheno)


def test_pharmcat_unknown_does_not_override_agreeing_callers():
    for unk in ("Unknown/Unknown", "unknown / unknown", "No Result"):
        rs = [_pc(unk), _r("PyPGx", "*2/*2"), _r("Stargazer", "*2/*2"), _r("Aldy", "*2/*2")]
        v = cmp.compute_verdict(rs, no_call=False, gene="CYP2B6")
        assert v["consensus_diplotype"] == "*2/*2" and v["status"] == "concordant", (unk, v)
        assert (v["n_agree"], v["n_called"]) == (3, 3), v
        assert "PharmCAT reported no diplotype" in v["note"], v
        assert v.get("authority") != "PharmCAT", v
    # callers that disagree stay disagreeing: the status says so
    rs = [_pc("Unknown/Unknown"), _r("PyPGx", "*2/*2"), _r("Stargazer", "*2/*2"), _r("Aldy", "*2/*2"),
          _r("StellarPGx", "*2/*5", pheno="Intermediate Metabolizer")]
    v = cmp.compute_verdict(rs, no_call=False, gene="CYP2B6")
    assert v["status"] == "majority" and (v["n_agree"], v["n_called"]) == (3, 4), v
    # one uncorroborated caller never becomes a concordant answer: explicit no-call
    v = cmp.compute_verdict([_pc("Unknown/Unknown"), _r("PyPGx", "*1/*3")], no_call=False, gene="CYP4F2")
    assert v["status"] == "no_call" and v["consensus_diplotype"] == "NO_CALL" and "PyPGx *1/*3" in v["note"], v


def test_pharmcat_real_call_stays_authority():
    rs = [_pc("*1/*6"), _r("PyPGx", "*1/*1"), _r("Stargazer", "*1/*1"), _r("Aldy", "*1/*1")]
    v = cmp.compute_verdict(rs, no_call=False, gene="CYP2B6")
    assert v["consensus_diplotype"] == "*1/*6" and v["authority"] == "PharmCAT", v
    assert (v["n_agree"], v["n_called"]) == (1, 4) and "note" in v and "no diplotype" not in v["note"], v
    # on any other gene a no-call abstains too (named in the note, not counted)
    v = cmp.compute_verdict([_pc("Unknown/Unknown"), _r("PyPGx", "*1/*2"), _r("Aldy", "*1/*2")], no_call=False, gene="NUDT15")
    assert v["n_called"] == 2 and "PharmCAT reported no diplotype" in v["note"], v


# ── lettered legacy names are never scored as their base star number ───────────
def test_lettered_legacy_names_have_no_activity_score():
    rc = cmp.reconcile
    for t in ("*14A/*1", "*14B/*1", "*13C/*1", "*4NX/*1", "*36X/*1", "*14a/*1", "*4n/*1", "*13c/*39"):
        base = t.split("/")[0]
        if base == "*36X":      # Aldy's bare tandem marker, not a legacy letter: normalised, see reconcile
            assert rc.clean_diplotype("CYP2D6", t) == "*36/*1", t
            continue
        assert rc.activity_score(t) is None, t
        assert rc.canonicalize_diplotype("CYP2D6", t) != rc.canonicalize_diplotype("CYP2D6", "*1/" + base.rstrip("ABCNXabcnx")), t


def test_copy_number_and_suballele_tokens_still_scored():
    import json
    rc = cmp.reconcile
    table = json.load(open(os.path.join(HERE, "cpic_cyp2d6.json")))["alleles"]
    assert rc.activity_score("*2X2/*1") == rc.activity_score("*2x2/*1") == table["*2x2"]["activity"] + table["*1"]["activity"]
    assert rc.activity_score("*1x2/*1") == table["*1x2"]["activity"] + table["*1"]["activity"]
    assert rc.activity_score("*4.021/*1") == table["*4"]["activity"] + table["*1"]["activity"]
    assert rc.activity_score("*68+*4/*1") == table["*68+*4"]["activity"] + table["*1"]["activity"]
    assert rc.activity_score("*36X+*10/*1") == table["*36+*10"]["activity"] + table["*1"]["activity"]
    # the CPIC table itself says why: *14 and *114 differ in function
    assert table["*14"]["activity"] != table["*114"]["activity"]


# ── phenotype agreement without diplotype agreement: all readings listed, none chosen ──
def test_phenotype_tier_lists_every_reading():
    rs = [_r("PyPGx", "*2/*2"), _r("Stargazer", "*2/*2"), _r("Aldy", "*2/*2"), _r("StellarPGx", "*2/*5")]
    v = cmp.compute_verdict(rs, no_call=False, gene="CYP2B6")
    assert v["consensus_diplotype"] == "*2/*2 or *2/*5", v
    assert v["consensus_phenotype"] == "Normal Metabolizer" and v["status"] == "concordant", v
    assert v["diplotype_resolved"] is False and (v["n_agree"], v["n_called"]) == (4, 4), v
    assert v["diplotype_readings"] == {"*2/*2": ["PyPGx", "Stargazer", "Aldy"], "*2/*5": ["StellarPGx"]}, v
    assert "agree on the phenotype but not the diplotype: *2/*2 (PyPGx, Stargazer, Aldy) or *2/*5 (StellarPGx)" in v["note"], v
    # order is by caller count, then alphabetical, whatever the input order
    v2 = cmp.compute_verdict(list(reversed(rs)), no_call=False, gene="CYP2B6")
    assert v2["consensus_diplotype"] == v["consensus_diplotype"], v2
    # 2 vs 2 tie, same phenotype: alphabetical
    v = cmp.compute_verdict([_r("PyPGx", "*1/*4"), _r("Aldy", "*1/*4"), _r("Stargazer", "*1/*3"),
                             _r("StellarPGx", "*1/*3")], no_call=False, gene="TPMT")
    assert v["consensus_diplotype"] == "*1/*3 or *1/*4" and v["diplotype_resolved"] is False, v


def test_phenotype_disagreement_stays_discordant_and_unanimous_unchanged():
    v = cmp.compute_verdict([_r("PyPGx", "*1/*4", pheno="Normal Metabolizer"), _r("Aldy", "*1/*4", pheno="Normal Metabolizer"),
                             _r("Stargazer", "*1/*2", pheno="Intermediate Metabolizer"),
                             _r("StellarPGx", "*1/*2", pheno="Intermediate Metabolizer")], no_call=False, gene="CYP2B6")
    assert v["status"] == "discordant" and "diplotype_readings" not in v, v
    v = cmp.compute_verdict([_r("PyPGx", "*1/*4"), _r("Aldy", "*4/*1")], no_call=False)
    assert v["consensus_diplotype"] == "*1/*4" and "diplotype_resolved" not in v, v

def test_dropped_pharmcat_no_call_leaks_neither_diplotype_nor_phenotype():
    # HG002 CYP2B6 as stored by a production run: PharmCAT "Unknown/Unknown" with phenotype "No Result"
    rs = [_r("PyPGx", "*2/*2", pheno="Normal Metabolizer"), _r("Stargazer", "*2/*2", pheno="normal_metabolizer"),
          _r("Aldy", "*2/*2", pheno="Normal"), _r("StellarPGx", "*2/*5", pheno="-"),
          _r("PharmCAT", "Unknown/Unknown", pheno="No Result")]
    v = cmp.compute_verdict(rs, no_call=False, gene="CYP2B6")
    assert v["consensus_phenotype"] == "Normal Metabolizer", v
    assert v["consensus_diplotype"] == "*2/*2 or *2/*5" and (v["n_agree"], v["n_called"]) == (3, 4), v
    assert "PharmCAT" not in v["diplotype_readings"].get("*2/*2", []) + v["diplotype_readings"].get("*2/*5", []), v
    # same with the phenotype column the only abstaining text, in the non-tier (vote) path
    rs = [_r("PyPGx", "*1/*2", pheno="Normal Metabolizer"), _r("Aldy", "*1/*2", pheno="Normal Metabolizer"),
          _r("PharmCAT", "Unknown/Unknown", pheno="No Result")]
    v = cmp.compute_verdict(rs, no_call=False, gene="CYP4F2")
    assert v["consensus_phenotype"] == "Normal Metabolizer" and v["n_called"] == 2, v

# ── review defects: abstentions, NAT2 readings, partial no-call, *36X, notes ──
def test_nat2_phenotype_tier_readings_keep_agreed_phenotype():
    rs = [_r("PyPGx", "*4/*5", pheno="Intermediate Metabolizer"), _r("Aldy", "*4/*5", pheno="Intermediate Metabolizer"),
          _r("Stargazer", "*4/*6", pheno="Intermediate Metabolizer")]
    v = cmp.compute_verdict(rs, no_call=False, gene="NAT2")
    assert v["consensus_diplotype"] == "*4/*5 or *4/*6" and v["consensus_phenotype"] == "Intermediate Metabolizer", v
    rs = [_r("PyPGx", "*4/*4", pheno="Rapid Metabolizer"), _r("Aldy", "*4/*4", pheno="Rapid Metabolizer"),
          _r("Stargazer", "*5/*6", pheno="Rapid Metabolizer")]       # readings differ in CPIC phenotype
    v = cmp.compute_verdict(rs, no_call=False, gene="NAT2")
    assert v["consensus_phenotype"] == "-" and v["status"] == cmp.SPLIT_STATUS, v
    assert "different CPIC phenotypes" in v["note"] and "confirm with a targeted assay" in v["note"], v


def test_partial_pharmcat_unknown_is_no_call_for_the_authority():
    rs = [_pc("*1/Unknown", "Indeterminate"), _r("PyPGx", "*1/*28"), _r("Aldy", "*1/*28")]
    v = cmp.compute_verdict(rs, no_call=False, gene="UGT1A1")
    assert v["consensus_diplotype"] == "*1/*28" and (v["n_agree"], v["n_called"]) == (2, 2), v
    assert "partial diplotype" in v["note"], v
    # one other caller only: not corroborated, PharmCAT's row stays
    v = cmp.compute_verdict([_pc("*1/Unknown"), _r("PyPGx", "*1/*28")], no_call=False, gene="UGT1A1")
    assert v["consensus_diplotype"] == "*1/Unknown" and v.get("authority") == "PharmCAT", v


def test_abstaining_caller_is_not_a_dissent_or_a_count():
    rs = [_r("PyPGx", "*1/*2"), _r("Aldy", "*1/*2"), _r("Stargazer", "Indeterminate", pheno="-")]
    v = cmp.compute_verdict(rs, no_call=False, gene="CYP2C9")
    assert v["consensus_diplotype"] == "*1/*2" and (v["n_agree"], v["n_called"]) == (2, 2), v
    assert v["status"] == "concordant" and v.get("diplotype_resolved", True) is True, v
    assert "Stargazer" in v["note"] and v["consensus_diplotype"] != "", v
    # PharmCAT Unknown on a non-authority gene abstains too (unanimous stays unanimous, 2-2 is not a tie)
    rs = [_pc("Unknown/Unknown", "n/a"), _r("PyPGx", "*1/*2"), _r("Aldy", "*1/*2")]
    v = cmp.compute_verdict(rs, no_call=False, gene="CYP2C9")
    assert (v["status"], v["n_agree"], v["n_called"]) == ("concordant", 2, 2), v
    # "-" / empty rows are absent, not abstentions: no note
    v = cmp.compute_verdict([_r("PharmCAT", "-"), _r("PyPGx", "*1/*2"), _r("Aldy", "*1/*2")], no_call=False, gene="CYP2B6")
    assert "not counted" not in v.get("note", ""), v


def test_aldy_tandem_marker_scores_like_call_cyp2d6():
    # Cyrius *1/*10 (no phenotype) and Aldy *1/*36X+*10 share the CPIC score but not the canonical diplotype
    rc = cmp.reconcile
    assert rc.activity_score("*1/*10") == rc.activity_score("*1/*36X+*10") is not None
    v = _d2d6("*1/*10", others=[("Aldy", "*1/*36X+*10", "-")])
    assert v["consensus_phenotype"] != "-" and "same score: Aldy" in v["phenotype_source"], v
    v = _d2d6("*1/*10", others=[("Stargazer", "*1/*13C", "-")])
    assert v["consensus_phenotype"] == "-", v


def test_no_call_is_never_shown_as_a_diplotype():
    # 1. PharmCAT-authority gene, PharmCAT no-call, the others do not corroborate each other
    v = cmp.compute_verdict([_pc("Unknown/Unknown", "n/a"), _r("PyPGx", "*1/*3"),
                             _r("Aldy", "*1/*4", pheno="Intermediate Metabolizer")], no_call=False, gene="CYP4F2")
    assert v["status"] == "no_call" and v["consensus_diplotype"] == "NO_CALL" and v["consensus_phenotype"] == "-", v
    assert "PharmCAT made no call; other callers do not corroborate each other: PyPGx *1/*3; Aldy *1/*4" in v["note"], v
    v = cmp.compute_verdict([_pc("Unknown/Unknown", "n/a")], no_call=False, gene="CYP2B6")
    assert v["status"] == "no_call", v
    # 2. every caller abstains on a non-PharmCAT gene
    v = cmp.compute_verdict([_r("PyPGx", "Indeterminate", pheno="-"), _r("Aldy", "Unknown/Unknown", pheno="-")],
                            no_call=False, gene="CYP2C9")
    assert v["status"] == "no_call" and v["consensus_diplotype"] == "NO_CALL" and "abstained" in v["note"], v
    # 3. the tier helper never returns an empty diplotype
    v = cmp._phenotype_tier_verdict([("x", "Indeterminate")], ["Aldy"], ("Normal Metabolizer", 1, "x"), 1, {})
    assert v["status"] == "no_call" and v["consensus_diplotype"] != "", v

# ── CYP2D6 verdict without Cyrius = the outside-call agreement rule ─────────────
def _vote(rows):
    return cmp.compute_verdict([_r(t, d, pheno=p) for t, d, p in rows], no_call=False, gene="CYP2D6")


def test_cyp2d6_vote_phenotype_is_the_outside_call_phenotype():
    import pgx_outside_calls as oc, tempfile
    cases = [
        # clinical-sample pattern: displayed diplotype *1x2/*36+*10 scores 2.25 (Normal); PyPGx/Aldy label Ultrarapid
        [("PyPGx", "*1x3/*10", "Ultrarapid Metabolizer"), ("Stargazer", "*1x2/*36+*10", "normal_metabolizer"),
         ("Aldy", "*1+*1/*10+*10", "Ultrarapid"), ("StellarPGx", "*36+*10/*1x2", "Normal metaboliser (NM)")],
        # no caller has a score for the unanimous allele: no phenotype, as before for the drug table
        [("PyPGx", "*43/*86", "Indeterminate"), ("Stargazer", "*43/*86", "-"), ("Aldy", "*43/*86", "Indeterminate")],
        # one caller Indeterminate blocks agreement
        [("PyPGx", "*1/*2", "Normal Metabolizer"), ("Aldy", "*1/*2x2", "Indeterminate")],
        # unanimous diplotype
        [("PyPGx", "*4/*10", "Intermediate Metabolizer"), ("Stargazer", "*4/*10", "intermediate_metabolizer")],
        # same score, different diplotypes
        [("PyPGx", "*1/*41", "Normal Metabolizer"), ("Aldy", "*2/*41", "Normal")],
        # lower-bound caller + exact callers of the same class
        [("PyPGx", "*1x3/*10", "Ultrarapid Metabolizer"), ("Aldy", "*2x2/*2x2", "Ultrarapid")],
    ]
    for rows in cases:
        v = _vote(rows)
        with tempfile.TemporaryDirectory() as out:
            d = os.path.join(out, "genes", "CYP2D6"); os.makedirs(d)
            with open(os.path.join(d, "CYP2D6_S_comparison.tsv"), "w") as fh:
                fh.write("Gene\tSample\tBuild\tTool\tDiplotype\tActivityScore\tPhenotype\tStatus\tSVMode\n")
                for t, dip, p in rows:
                    fh.write(f"CYP2D6\tS\tGRCh38\t{t}\t{dip}\t-\t{p}\tok\t-\n")
            o = oc.call_cyp2d6(out)
        # the gene row has a phenotype exactly when the drug table gets a call, and it is the same one
        exp = oc.reconcile.agreement_phenotype(o) if o["level"] != "-" else None
        assert v["consensus_phenotype"] == (exp or "-"), (rows, v, o)
    v = _vote(cases[0])                       # that case: no longer a tie broken by caller order
    assert v["consensus_phenotype"] == "-" and "phenotype not derived" in v["phenotype_note"], v
    assert v["consensus_diplotype"] == "*1x2/*36+*10", v
    v = _vote(cases[3])
    assert v["consensus_phenotype"] == "Intermediate Metabolizer" and v["phenotype_source"].startswith("CPIC activity score 0.25"), v
    # an all-Ultrarapid vote including a lower-bound caller
    assert _vote(cases[5])["consensus_phenotype"] == "Ultrarapid Metabolizer"


# ── TPMT / NUDT15: phenotype from the CPIC diplotype table ────────────────────
def test_thiopurine_phenotype_is_the_cpic_table_row():
    import json
    t = json.load(open(os.path.join(HERE, "cpic_thiopurine.json")))
    # rows quoted in the report (CPIC API /diplotype, thiopurine guideline lookup)
    assert t["TPMT"]["*1/*6"] == "Intermediate Metabolizer" and t["NUDT15"]["*1/*6"] == "Intermediate Metabolizer"
    assert t["NUDT15"]["*3/*6"] == "Poor Metabolizer"
    v = cmp.compute_verdict([_r("PyPGx", "*1/*6", pheno="Indeterminate"), _r("Aldy", "*1/*6", pheno="Indeterminate"),
                             _r("StellarPGx", "*1/*6", pheno="-")], no_call=False, gene="TPMT")   # HG006
    assert v["consensus_phenotype"] == "Intermediate Metabolizer" and v["status"] == "concordant", v
    assert v["phenotype_source"] == "CPIC TPMT diplotype table, row *1/*6", v
    v = cmp.compute_verdict([_r("PyPGx", "*1/*1"), _r("Aldy", "*1/*1S", pheno="Normal")], no_call=False, gene="TPMT")
    assert v["consensus_phenotype"] == "Normal Metabolizer", v      # *1S is *1
    v = cmp.compute_verdict([_r("PyPGx", "*6/*1", pheno="Indeterminate"), _r("Aldy", "*1/*6", pheno="Indeterminate")],
                            no_call=False, gene="NUDT15")
    assert v["consensus_phenotype"] == "Intermediate Metabolizer", v
    # a caller label that contradicts CPIC is replaced and named
    v = cmp.compute_verdict([_r("PyPGx", "*1/*6", pheno="Normal Metabolizer"), _r("Aldy", "*1/*6", pheno="Normal")],
                            no_call=False, gene="NUDT15")
    assert v["consensus_phenotype"] == "Intermediate Metabolizer" and "differs from the CPIC phenotype" in v["note"], v
    # a discordant verdict shows no single caller's reading as the CPIC row
    v = cmp.compute_verdict([_r("PyPGx", "*1/*3", pheno="Intermediate Metabolizer"), _r("Stargazer", "*2/*6", pheno="-"),
                             _r("Aldy", "*5/*6", pheno="Indeterminate"), _r("StellarPGx", "*1/*1", pheno="-")],
                            no_call=False, gene="NUDT15")
    assert "phenotype_source" not in v, v
    assert v["status"] == "discordant" and v["consensus_phenotype"] == "-", v    # HG005 shape: PyPGx label dropped
    v = cmp.compute_verdict([_r("PyPGx", "*1/*5", pheno="Intermediate Metabolizer"), _r("Aldy", "*5/*6", pheno="-")],
                            no_call=False, gene="NUDT15")
    assert v["consensus_phenotype"] == "-" and "caller labels: PyPGx Intermediate Metabolizer" in v["note"], v


def test_nudt15_allele_split_withholds_the_phenotype():
    rs = [_r("PyPGx", "*1/*1", pheno="Normal Metabolizer"), _r("Stargazer", "*1/*6", pheno="-", status="not_run"),
          _r("Aldy", "*1/*6", pheno="Indeterminate"), _r("StellarPGx", "*1/*1", pheno="-")]
    v = cmp.compute_verdict(rs, no_call=False, gene="NUDT15", informational={"PharmCAT": "*1/*6"})   # HG006
    assert v["consensus_phenotype"] == "-" and "phenotype_source" not in v, v
    assert "callers disagree on NUDT15 *6: Aldy, PharmCAT report it; PyPGx, StellarPGx do not" in v["note"], v
    v = cmp.compute_verdict(rs, no_call=False, gene="NUDT15")          # without the PharmCAT row: still split
    assert v["consensus_phenotype"] == "-" and "callers disagree on NUDT15 *6: Aldy report it; PyPGx, StellarPGx do not" in v["note"], v
    # callers that agree on every allele keep the CPIC lookup
    v = cmp.compute_verdict([_r("PyPGx", "*1/*6"), _r("Aldy", "*6/*1")], no_call=False, gene="NUDT15",
                            informational={"PharmCAT": "*1/*6"})
    assert v["consensus_phenotype"] == "Intermediate Metabolizer", v


def test_cyp2d6_majority_by_score_level_m():
    # clinical-sample pattern: three callers *1/*41 (1.25), Aldy has an unscored allele; HG002_12X pattern: two Indeterminate
    v = _vote([("PyPGx", "*1/*41", "Normal Metabolizer"), ("Stargazer", "*1/*41", "normal_metabolizer"),
               ("Aldy", "*1+*2/*119", "Ultrarapid"), ("StellarPGx", "*1/*41", "Normal metaboliser (NM)")])
    assert v["consensus_phenotype"] == "Normal Metabolizer" and "majority by activity score 1.25 (3/4 callers)" in v["phenotype_source"], v
    assert "dissenting: Aldy=*1+*2/*119" in v["note"], v
    v = _vote([("PyPGx", "*2/*4", "Intermediate Metabolizer"), ("Stargazer", "*2/*4", "intermediate_metabolizer"),
               ("Aldy", "*4K/*160", "-"), ("StellarPGx", "*39/*4", "-")])
    assert v["consensus_phenotype"] == "Intermediate Metabolizer", v    # *39/*4 also scores 1.0: three callers share it
    v = _vote([("PyPGx", "*2/*4", "x"), ("Stargazer", "*2/*4", "x"), ("Aldy", "*1/*1", "Indeterminate")])
    assert v["consensus_phenotype"] == "-", v      # only two callers share a score: below MAJORITY_MIN
    v = _vote([("PyPGx", "*1x3/*10", "-"), ("Stargazer", "*1/*41", "-"), ("Aldy", "*1/*41", "-")])
    assert v["consensus_phenotype"] == "-", v      # two share a score
    import pgx_outside_calls as oc
    assert oc.reconcile.MAJORITY_MIN == 3


def test_report_key_findings_keep_authority_genes():
    src = open(os.path.join(HERE, "pgx-report.py")).read()
    assert 'if n_agree < 2 and not auth_tool:' in src and '"authority_tool"' in src
    assert 'concordance_warn = n_agree < n_called\n' in src      # an authority gene keeps the disagreement warning
    assert "decided_by" not in src


# ── authority genes: the verdict says who decided ─────────────────────────────
def test_authority_verdict_names_who_decided():
    v = cmp.compute_verdict([_r("VCF-Check", "*1/*2", pheno="Normal"), _r("PyPGx", "*1/*1"), _r("Aldy", "*1/*1")],
                            no_call=False, gene="VKORC1")
    assert v["status"] == "concordant" and v["n_agree"] == 1 and v["authority"] == "VCF-Check", v
    assert "decided_by" not in v, v                       # the one key is `authority`
    v = cmp.compute_verdict([_pc("*1/*28"), _r("PyPGx", "*1/*1")], no_call=False, gene="UGT1A1")
    assert v.get("authority") == "PharmCAT", v
    v = cmp.compute_verdict([_r("PyPGx", "*1/*1"), _r("Aldy", "*1/*1")], no_call=False, gene="CYP2C9")
    assert "authority" not in v, v                        # a vote is not an authority decision
    snake = open(os.path.join(HERE, "Snakefile")).read()
    assert 'v.get("authority") in (None, "", "Phenotype")' in snake and "DecidedBy" in snake


def test_majority_phenotype_comes_only_from_the_majority_callers():
    v = cmp.compute_verdict([_r("PyPGx", "*1/*2", pheno="-"), _r("Aldy", "*1/*2", pheno="-"),
                             _r("Stargazer", "*1/*3", pheno="Poor Metabolizer")], no_call=False, gene="CYP2C9")
    assert v["status"] == "majority" and v["consensus_diplotype"] == "*1/*2" and v["consensus_phenotype"] == "-", v
    v = cmp.compute_verdict([_r("PyPGx", "*1/*2", pheno="Intermediate Metabolizer"), _r("Aldy", "*1/*2", pheno="-"),
                             _r("Stargazer", "*1/*3", pheno="Poor Metabolizer")], no_call=False, gene="CYP2C9")
    assert v["consensus_phenotype"] == "Intermediate Metabolizer", v


def test_nudt15_zygosity_and_spelling():
    rs = [_r("PyPGx", "*3/*3", pheno="-"), _r("Aldy", "*1/*3", pheno="-")]
    v = cmp.compute_verdict(rs, no_call=False, gene="NUDT15")
    assert v["consensus_phenotype"] == "-" and "callers disagree on NUDT15 *3 zygosity: PyPGx reports 2 copies; Aldy reports 1 copy" in v["note"], v
    # a sub-allele suffix and Aldy's X marker are spelling, not a split
    v = cmp.compute_verdict([_r("PyPGx", "*3.001/*1"), _r("Aldy", "*1/*3"), _r("Stargazer", "*3X/*1")],
                            no_call=False, gene="NUDT15")
    assert v["consensus_phenotype"] == "Intermediate Metabolizer", v


def test_cyp2d6_lower_bound_cyrius_agrees_with_outside_call():
    import pgx_outside_calls as oc, tempfile
    rows = [("Cyrius", "*1x3/*10", "-"), ("PyPGx", "*1x3/*10", "-"), ("Aldy", "*1/*10", "Normal")]
    v = cmp.compute_verdict([_r(t, d, pheno=p) for t, d, p in rows], no_call=False, gene="CYP2D6")
    with tempfile.TemporaryDirectory() as out:
        d = os.path.join(out, "genes", "CYP2D6"); os.makedirs(d)
        with open(os.path.join(d, "CYP2D6_S_comparison.tsv"), "w") as fh:
            fh.write("Gene\tSample\tBuild\tTool\tDiplotype\tActivityScore\tPhenotype\tStatus\tSVMode\n")
            for t, dip, p in rows:
                fh.write(f"CYP2D6\tS\tGRCh38\t{t}\t{dip}\t-\t{p}\tok\t-\n")
        o = oc.call_cyp2d6(out)
    assert v["consensus_phenotype"] == "Ultrarapid Metabolizer", v
    assert o["level"] == "e" and o["diplotype"] == "*1x3/*10" and o["score"] is None, o   # never an exact score
    assert oc.reconcile.agreement_phenotype(o) == v["consensus_phenotype"], (o, v)


def test_level_m_is_blocked_by_an_explicit_indeterminate():
    v = _vote([("PyPGx", "*1/*41", "Normal Metabolizer"), ("Stargazer", "*1/*41", "normal_metabolizer"),
               ("StellarPGx", "*1/*41", "Normal metaboliser (NM)"), ("Aldy", "*1/*1", "Indeterminate")])
    assert v["consensus_phenotype"] == "-", v        # Aldy scores 2.0 and says Indeterminate
    v = _vote([("PyPGx", "*2/*4", "Intermediate Metabolizer"), ("Stargazer", "*2/*4", "intermediate_metabolizer"),
               ("Aldy", "*4K/*160", "Indeterminate"), ("StellarPGx", "*39/*4", "Indeterminate")])
    assert v["consensus_phenotype"] == "-", v        # HG002_12X shape: StellarPGx's scoreable Indeterminate blocks m


# ── owner decision 2026-11-15: readings with different CPIC phenotypes -> "-" + status label ──
def _tier(gene, readings, label="Poor Metabolizer"):
    """Callers whose own phenotype labels agree (so the phenotype tier fires) but whose diplotypes are `readings`."""
    tools = ["PyPGx", "Stargazer", "Aldy", "StellarPGx"]
    return cmp.compute_verdict([_r(t, readings[i % len(readings)], pheno=label) for i, t in enumerate(tools)],
                               no_call=False, gene=gene)


def test_split_cpic_phenotype_has_no_phenotype_and_its_own_status():
    for gene, rd, a, b in (("CYP3A5", ["*3/*3", "*1/*1"], "Poor Metabolizer", "Normal Metabolizer"),
                           ("CYP3A5", ["*1/*3", "*1/*1"], "Intermediate Metabolizer", "Normal Metabolizer"),
                           ("CYP2C19", ["*1/*2", "*2/*36"], "Intermediate Metabolizer", "Poor Metabolizer")):
        v = _tier(gene, rd)
        assert v["consensus_phenotype"] == "-" and v["status"] == "concordant, diplotype unresolved", v
        assert v["diplotype_resolved"] is False and not v.get("phenotype_tier"), v
        assert f"{rd[0]} {a}" in v["note"] and f"{rd[1]} {b}" in v["note"] and "confirm with a targeted assay" in v["note"], v
        assert v["phenotype_by_reading"] == {rd[0]: a, rd[1]: b}, v


def test_same_cpic_phenotype_readings_keep_phenotype_and_status():
    for gene, rd, ph in (("CYP2C19", ["*1/*2", "*1/*2B"], "Intermediate Metabolizer"),
                         ("CYP3A5", ["*3/*3", "*3/*3B"], "Poor Metabolizer"),
                         ("CYP2B6", ["*2/*2", "*2/*5"], "Normal Metabolizer")):
        v = _tier(gene, rd, label=ph)
        assert v["consensus_phenotype"] == ph and v["status"] == "concordant" and v.get("phenotype_tier"), v
        assert "phenotype_split" not in v, v


import contextlib


@contextlib.contextmanager
def _stub(obj, name, value):
    """Replace obj.name for the block and always restore the original."""
    orig = getattr(obj, name)
    setattr(obj, name, value)
    try:
        yield
    finally:
        setattr(obj, name, orig)


def _dpyd(readings_by_tool, label="Intermediate Metabolizer"):
    return cmp.compute_verdict([_r(t, d, pheno=label) for t, d in readings_by_tool.items()], no_call=False, gene="DPYD")


_STAR2A = "c.1905+1G>A (*2A)/c.85T>C (*9A)"
_ALDY2A = "*2A+rs1801160+rs45589337+rs1801265/*rs2297595"     # *2A (no function) + normal-function tags


def test_unscorable_reading_set_aside_only_after_the_component_check():
    # owner policy 2026-11-16: a composite Aldy reading whose components are all known and add no deviating allele
    v = _dpyd({"PyPGx": _STAR2A, "Stargazer": "*S10/*S3", "Aldy": _ALDY2A, "PharmCAT": _STAR2A})
    assert v["consensus_phenotype"] == "Intermediate Metabolizer" and v["status"] == "concordant", v
    assert "phenotype_split" not in v and v["unscorable_readings_ignored"] == [_ALDY2A], v
    assert f"reading {_ALDY2A} has no CPIC phenotype and was not counted" in v["note"], v
    assert "every component has a CPIC allele function" in v["note"], v
    assert v["note"].count("CPIC phenotype of the other readings is Intermediate Metabolizer") == 1, v
    assert "_ALDY2A" not in v["phenotype_source"] and _ALDY2A not in v["phenotype_source"], v
    # a differing caller label is replaced, named first, then the clause
    v = _dpyd({"PyPGx": _STAR2A, "Stargazer": "*S10/*S3", "Aldy": _ALDY2A, "PharmCAT": _STAR2A}, label="Normal Metabolizer")
    assert v["consensus_phenotype"] == "Intermediate Metabolizer", v
    assert v["note"].index("callers said Normal Metabolizer") < v["note"].index("was not counted"), v


def _split(support, label="Normal Metabolizer", gene="DPYD"):
    """_split_by_cpic_phenotype on a hand-built "A or B" verdict; support = {reading: [callers]}. (compute_verdict's own
    comparison key strips `+rsNNN` tags, so a lone-tag composite never reaches this rule through it.)"""
    return cmp._split_by_cpic_phenotype({"consensus_diplotype": " or ".join(support), "consensus_phenotype": label,
                                         "status": "concordant", "phenotype_tier": True, "authority": "Phenotype",
                                         "diplotype_readings": support}, gene)


def test_unscorable_reading_with_an_extra_deviating_component_is_a_split():
    # Aldy *1/*1+rs67376798 (c.2846A>T, decreased function) beside *1/*1: evidence the others lack
    assert cmp.reconcile.reading_function_profile("DPYD", "*1/*1+rs67376798")[0], "c.2846A>T must be a deviating component"
    v = _split({"*1/*1": ["PyPGx", "Stargazer", "PharmCAT"], "*1/*1+rs67376798": ["Aldy"]})
    assert v["consensus_phenotype"] == "-" and v["status"] == cmp.SPLIT_STATUS and "unscorable_readings_ignored" not in v, v
    assert "*1/*1+rs67376798 no CPIC phenotype" in v["note"], v
    # the same composite with only normal-function tags is set aside
    v = _split({"*1/*1": ["PyPGx", "Stargazer", "PharmCAT"], "*1/*1+rs1801159": ["Aldy"]})
    assert v["consensus_phenotype"] == "Normal Metabolizer" and v["status"] == "concordant", v
    # a component whose function CPIC's table does not give is unknown: split
    assert cmp.reconcile.reading_function_profile("DPYD", "*1/*1+rs0000001")[1] == ["*rs0000001"]
    v = _split({"*1/*1": ["PyPGx", "Stargazer", "PharmCAT"], "*1/*1+rs0000001": ["Aldy"]})
    assert v["consensus_phenotype"] == "-" and v["status"] == cmp.SPLIT_STATUS, v
    # one copy in the kept reading does not cover two in the unscorable one
    v = _split({"*1/*2A": ["PyPGx", "Stargazer", "PharmCAT"], "*2A+rs1801159/*2A+rs2297595": ["Aldy"]},
               label="Intermediate Metabolizer")
    assert v["consensus_phenotype"] == "-" and v["status"] == cmp.SPLIT_STATUS, v


def test_extra_variant_tag_votes_only_when_its_cpic_function_is_normal():
    R = cmp.reconcile
    four = lambda aldy: [_r("PyPGx", "*1/*1"), _r("Stargazer", "*1/*1"), _r("PharmCAT", "*1/*1"), _r("Aldy", aldy)]
    # Aldy *1/*1+rs67376798 (DPYD c.2846A>T, Decreased function) no longer votes as *1/*1
    assert R.canonicalize_diplotype("DPYD", "*1/*1+rs67376798") != R.canonicalize_diplotype("DPYD", "*1/*1")
    assert [x[1:] for x in R.canonical_tags("DPYD", "*1/*1+rs67376798")] == [("rs67376798", "Decreased function")]
    v = cmp.compute_verdict(four("*1/*1+rs67376798"), no_call=False, gene="DPYD")
    assert v["consensus_phenotype"] == "-" and v["status"] != "concordant", v          # was a 4/4 concordant Normal
    assert "Aldy *1/*1+rs67376798: tag rs67376798 is Decreased function in CPIC's allele table" in v["note"], v
    v = cmp.compute_verdict([_r("PyPGx", "*1/*1"), _r("Stargazer", "*1/*1"), _r("PharmCAT", "*1/*1"),
                             _r("Aldy", "*1/*1+rs67376798", pheno="Intermediate Metabolizer")], no_call=False, gene="DPYD")
    assert v["consensus_phenotype"] != "Normal Metabolizer" or v["status"] != "concordant", v
    # all-Normal tags: unchanged, still one reading, no note
    v = cmp.compute_verdict(four("*1/*1+rs1801159"), no_call=False, gene="DPYD")
    assert v["status"] == "concordant" and (v["n_agree"], v["n_called"]) == (4, 4) and "tag" not in v.get("note", ""), v
    # a tag CPIC's table does not give: its own reading
    v = cmp.compute_verdict(four("*1/*1+rs0000001"), no_call=False, gene="DPYD")
    assert v["consensus_phenotype"] == "-" and "tag rs0000001 is unknown" in v["note"], v
    # a gene with a synonym map but no function table keeps every tag (nothing can be shown Normal)
    with _stub(R, "_allele_functions", lambda g: None):
        assert R.canonicalize_diplotype("DPYD", "*1/*1+rs1801159") != R.canonicalize_diplotype("DPYD", "*1/*1")
        assert [x[1:] for x in R.canonical_tags("DPYD", "*1/*1+rs1801159")] == [("rs1801159", "unknown")]
    # two callers on the tagged reading are a reading of their own, not a vote for *1/*1
    v = cmp.compute_verdict([_r("PyPGx", "*1/*1"), _r("PharmCAT", "*1/*1"),
                             _r("Aldy", "*1/*1+rs67376798"), _r("Stargazer", "*1/*1+rs67376798")], no_call=False, gene="DPYD")
    assert v["consensus_phenotype"] == "-" and v["status"] != "concordant", v


def test_tag_vote_review_fixes():
    R = cmp.reconcile
    c = lambda d, g="DPYD": R.canonicalize_diplotype(g, d)
    # (1) an HGVS tag keeps its own `+` and resolves (c.1905+1G>A is *2A, No function)
    assert R._split_tags("*1+c.1905+1G>A") == ("*1", ["c.1905+1g>a"]), R._split_tags("*1+c.1905+1G>A")
    assert c("*1/*1+c.1905+1G>A") == c("*1/*1+rs3918290") != c("*1/*1")
    assert [x[2] for x in R.canonical_tags("DPYD", "*1/*1+c.1905+1G>A")] == ["No function"]
    assert R._strip_allele("c.1905+1G>A") == "c.1905+1g>a"
    # (4) only the blocking tags are in the key; a Normal tag beside them is dropped
    assert c("*1/*1+rs1801159+rs67376798") == c("*1/*1+rs67376798")
    # (5) one variant in three notations is one key
    assert c("*1/*1+rs67376798") == c("*1/*1+*rs67376798") == c("*1/*1+c.2846A>T") != c("*1/*1")
    # (3) the "not counted" note only names a caller whose key differs from the consensus reading
    same = cmp.compute_verdict([_r(t, "*1/*1+rs67376798") for t in ("PyPGx", "Stargazer", "PharmCAT", "Aldy")], no_call=False, gene="DPYD")
    assert "not counted as agreeing" not in same["note"] and "keeps a tag that is not Normal function" in same["note"], same
    alone = cmp.compute_verdict([_r("PyPGx", "*1/*1"), _r("Stargazer", "*1/*1"), _r("PharmCAT", "*1/*1"),
                                 _r("Aldy", "*1/*1+rs67376798", pheno="Intermediate Metabolizer")], no_call=False, gene="DPYD")
    assert alone["status"] == "majority" and alone["note"].count("not counted as agreeing") == 1 and "Aldy *1/*1+rs67376798" in alone["note"], alone
    assert "PyPGx *1/*1:" not in alone["note"], alone
    # (6) TPMT: a Normal-function tag does not stop the CPIC lookup; a tag the table does not give does (and is named)
    assert R.thiopurine_phenotype("TPMT", "*1/*3A+*1") == R.thiopurine_phenotype("TPMT", "*1/*3A") is not None
    assert R.thiopurine_phenotype("TPMT", "*1/*3A+rs1142345") is None
    v = cmp.compute_verdict([_r("PyPGx", "*1/*3A+rs1142345"), _r("Stargazer", "*1/*3A")], no_call=False, gene="TPMT")
    assert "rs1142345 is unknown in CPIC's allele table" in v["note"], v
    dock = open(os.path.join(HERE, "..", "Dockerfile"), encoding="utf-8").read()
    for f in ("cpic_diplotype.json", "dpyd_notation_map.json", "allele_synonyms.json"):
        assert f"COPY docker/{f}" in dock, f      # canonicalize_allele needs these at run time; a missing one fails loud
    # (7) the tag list comes from the keys the vote used: an unplaceable call and a VCF-check call have none
    assert R.canonical_tags("DPYD", "*1(het)+*2A(het) [unphased]") == [] and R.canonicalize_diplotype("DPYD", "*1(het)+*2A(het) [unphased]") is R.UNRESOLVED
    assert R.canonical_tags("DPYD", "*2A(het)") == [] and R.canonical_tags("DPYD", "-") == []
    # (8) the cleaned component is returned by _component_function
    assert R._component_function("DPYD", "rs0000001")[0] == "*rs0000001"


def test_unscorable_reading_backed_by_more_callers_is_a_split():
    base = {"consensus_diplotype": "*1/*1 or *1/*1+rs1801159", "consensus_phenotype": "Normal Metabolizer",
            "status": "concordant", "phenotype_tier": True, "authority": "Phenotype"}
    v = cmp._split_by_cpic_phenotype(dict(base, diplotype_readings={"*1/*1": ["PyPGx"], "*1/*1+rs1801159": ["Stargazer", "Aldy"]}), "DPYD")
    assert v["consensus_phenotype"] == "-" and v["status"] == cmp.SPLIT_STATUS, v
    v = cmp._split_by_cpic_phenotype(dict(base, diplotype_readings={"*1/*1": ["PyPGx", "Stargazer"], "*1/*1+rs1801159": ["Aldy", "PharmCAT"]}), "DPYD")
    assert v["consensus_phenotype"] == "Normal Metabolizer" and v["status"] == "concordant", v      # equal support is not more
    v = cmp._split_by_cpic_phenotype(dict(base), "DPYD")        # support unknown: cannot be shown, so a split
    assert v["consensus_phenotype"] == "-" and v["status"] == cmp.SPLIT_STATUS, v


def test_ambiguous_cpic_row_is_a_split_never_ignored():
    assert cmp.reconcile.cpic_diplotype_lookup("DPYD", "*1/*1") == ("scored", "Normal Metabolizer")
    assert cmp.reconcile.cpic_diplotype_lookup("DPYD", "*1/*1+rs1801159")[0] == "no_row"
    def lookup(g, d):
        return ("ambiguous", None) if d == "*1/*2" else ("scored", "Normal Metabolizer")
    with _stub(cmp.reconcile, "cpic_diplotype_lookup", lookup):
        v = cmp.compute_verdict([_r("PyPGx", "*1/*1"), _r("Stargazer", "*1/*1"), _r("Aldy", "*1/*2"), _r("PharmCAT", "*1/*1")],
                                no_call=False, gene="CYP2C19")
    assert v["consensus_phenotype"] == "-" and v["status"] == cmp.SPLIT_STATUS, v
    assert v["phenotype_by_reading"]["*1/*2"] is None, v
    assert cmp._reading_phenotype("CYP2C19", "*1/*2")[2] == "scored"


def test_nat2_indeterminate_reading_is_a_split_never_ignored():
    with _stub(cmp, "nat2_cpic_phenotype", lambda d: "Indeterminate" if d == "*4/*99" else "Normal Metabolizer"):
        assert cmp._reading_phenotype("NAT2", "*4/*99") == (True, None, "indeterminate")
        v = cmp.compute_verdict([_r("PyPGx", "*4/*4"), _r("Stargazer", "*4/*4"), _r("Aldy", "*4/*99"), _r("PharmCAT", "*4/*4")],
                                no_call=False, gene="NAT2")
    assert v["consensus_phenotype"] == "-" and v["status"] == cmp.SPLIT_STATUS, v


def test_lettered_unscorable_allele_still_counts_toward_the_split():
    # CYP2C19 *1A is deliberately unscorable (owner decision 2026-11-16) and has no CPIC allele function
    assert cmp.reconcile.reading_function_profile("CYP2C19", "*1A/*2")[1] == ["*1A"]
    v = _tier("CYP2C19", ["*1/*2", "*1A/*2"], label="Intermediate Metabolizer")
    assert v["consensus_phenotype"] == "-" and v["status"] == cmp.SPLIT_STATUS and "*1A/*2 no CPIC phenotype" in v["note"], v
    # unknown *99 is no allele CPIC defines: split
    v = _tier("CYP2C19", ["*1/*2", "*1/*99"], label="Intermediate Metabolizer")
    assert v["consensus_phenotype"] == "-" and v["status"] == cmp.SPLIT_STATUS, v


def test_unscorable_policy_constant_and_split_rules():
    orig = cmp.UNSCORABLE_READING_IS_A_SPLIT
    assert orig is False
    try:
        cmp.UNSCORABLE_READING_IS_A_SPLIT = True            # the 2026-11-15 rule: any unscorable reading splits
        v = _dpyd({"PyPGx": _STAR2A, "Stargazer": "*S10/*S3", "Aldy": _ALDY2A, "PharmCAT": _STAR2A})
        assert v["consensus_phenotype"] == "-" and v["status"] == cmp.SPLIT_STATUS, v
    finally:
        cmp.UNSCORABLE_READING_IS_A_SPLIT = orig
    # scorable readings that disagree: still a split
    v = _tier("CYP2C19", ["*1/*2", "*2/*36", "*1/*99"], label="Intermediate Metabolizer")
    assert v["consensus_phenotype"] == "-" and v["status"] == cmp.SPLIT_STATUS, v
    # the only scorable reading has no usable phenotype: still a split
    with _stub(cmp.reconcile, "cpic_diplotype_lookup", lambda g, d: ("scored", "Indeterminate") if d == "*1/*2" else ("no_row", None)):
        v = _tier("CYP2C19", ["*1/*2", "*1/*99"], label="Intermediate Metabolizer")
    assert v["consensus_phenotype"] == "-" and v["status"] == cmp.SPLIT_STATUS, v
    # a gene with no CPIC diplotype table is unchanged
    v = _tier("CYP4F2", ["*1/*3", "*3/*3"], label="Intermediate Metabolizer")
    assert v["consensus_phenotype"] == "Intermediate Metabolizer" and v["status"] == "concordant", v


def test_report_sentence_for_a_set_aside_reading():
    import importlib.util
    spec = importlib.util.spec_from_file_location("pgx_report_t", os.path.join(HERE, "pgx-report.py"))
    rep = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rep)
    v = _dpyd({"PyPGx": _STAR2A, "Stargazer": "*S10/*S3", "Aldy": _ALDY2A, "PharmCAT": _STAR2A})
    tools = {t: {"status": "ok", "diplotype": d, "phenotype": "Intermediate Metabolizer"}
             for t, d in {"PyPGx": _STAR2A, "Stargazer": "*S10/*S3", "Aldy": _ALDY2A, "PharmCAT": _STAR2A}.items()}
    _, res = rep._pipeline_narrative("DPYD", v, tools)
    assert "set aside" in res and _ALDY2A in res and "were not all in agreement" in res, res
    assert "every reading agreed" not in res and "all callers agreed" not in res, res


def test_phase_ambiguous_cyp2b6_follows_cpic_phenotypes():
    # *5/*6 and *1/*7 hold the same variants; CPIC scores both Intermediate -> phenotype kept, status unchanged
    rs = [_pc("*5/*6", "Intermediate Metabolizer"), _r("PyPGx", "*1/*7", pheno="Intermediate Metabolizer"),
          _r("Stargazer", "*5/*6", pheno="Intermediate Metabolizer"), _r("Aldy", "*5/*6", pheno="Intermediate Metabolizer")]
    v = cmp.compute_verdict(rs, no_call=False, gene="CYP2B6")
    assert v.get("phase_ambiguous") and v["consensus_phenotype"] == "Intermediate Metabolizer", v
    assert v["status"] == "concordant", v
    rc = cmp.reconcile
    assert rc.cpic_diplotype_phenotype("CYP2B6", "*1/*6") == rc.cpic_diplotype_phenotype("CYP2B6", "*4/*9") \
        == "Intermediate Metabolizer"
    # a phase-ambiguous pair CPIC scores differently gets "-" and the status label (stub the table lookup)
    with _stub(cmp.reconcile, "cpic_diplotype_lookup",
               lambda g, d: ("scored", "Poor Metabolizer" if d == "*1/*7" else "Intermediate Metabolizer")):
        v = cmp.compute_verdict(rs, no_call=False, gene="CYP2B6")
    assert v["consensus_phenotype"] == "-" and v["status"] == cmp.SPLIT_STATUS and "phase_ambiguous" in v, v


def test_cpic_diplotype_table_lookup_is_cpic_data():
    rc = cmp.reconcile
    assert rc.cpic_diplotype_phenotype("CYP3A5", "*3/*3") == "Poor Metabolizer"
    assert rc.cpic_diplotype_phenotype("CYP2C19", "*36/*2") == "Poor Metabolizer"       # order-insensitive
    assert rc.cpic_diplotype_phenotype("CYP2D6", "*1/*1") is None                       # CYP2D6 has its own rule
    assert rc.cpic_diplotype_phenotype("CYP4F2", "*1/*3") is None                       # no CPIC table


def test_dpyd_synonym_readings_are_one_reading():
    # PyPGx HGVS+star name, Stargazer S-name and Aldy star name of the same variant are one reading
    rc = cmp.reconcile
    assert len({rc.canonicalize_diplotype("DPYD", d) for d in
                ("*1/c.1905+1G>A (*2A)", "*1/*2A", "Reference/c.1905+1G>A", "*S1/*1")}) == 1
    rs = [_r("PyPGx", "*1/c.1905+1G>A (*2A)", pheno="Intermediate Metabolizer"),
          _r("Stargazer", "*1/*S1", pheno="Intermediate Metabolizer"), _r("Aldy", "*1/*2A", pheno="Intermediate Metabolizer"),
          _r("PharmCAT", "*1/*2A", pheno="Intermediate Metabolizer")]
    v = cmp.compute_verdict(rs, no_call=False, gene="DPYD")
    assert " or " not in v["consensus_diplotype"] and v["status"] == "concordant" and (v["n_agree"], v["n_called"]) == (4, 4), v


def test_every_synonym_pair_has_the_same_cpic_function():
    """A synonym folds two names into one allele: where CPIC defines both, their function must be identical."""
    import json
    fn = json.load(open(os.path.join(HERE, "cpic_diplotype.json")))["allele_function"]
    syn = json.load(open(os.path.join(HERE, "allele_synonyms.json")))
    n = 0
    for gene, ent in syn.items():
        if gene.startswith("_"):
            continue
        table = {k.lower(): v for k, v in fn.get(gene, {}).items()}
        for canon, aliases in ent.items():
            group = [canon] + aliases
            # CPIC names carry the HGVS and star name together (`c.85T>C (*9A)`): index each part
            known = {}
            for name, f in fn.get(gene, {}).items():
                for part in [name] + re.findall(r"\(([^)]*)\)", name) + [re.sub(r"\s*\(.*\)", "", name)]:
                    known[part.strip().lower()] = f
            fs = {known[t.lower()] for t in group if t.lower() in known}
            assert len(fs) <= 1, (gene, group, fs)
            n += len(fs)
    assert n >= 4, n          # CPIC defines DPYD *9A, *5, *2A, *13 and c.2846A>T, CYP2C19 *2: the check is not vacuous
    assert syn["CYP2C19"]["*2"] == [] and fn["CYP2C19"]["*38"] != fn["CYP2C19"]["*2"]


import reconcile as R


def test_dpyd_notation_map_matches_its_sources():
    """dpyd_notation_map.json is rebuilt from Stargazer's shipped tables + the embedded CPIC allele definitions and must
    equal the committed file; CPIC API comparison under PGX_CPIC_VERIFY=1."""
    import json, importlib.util
    spec = importlib.util.spec_from_file_location("nmb", os.path.join(HERE, "dpyd_notation_map_build.py"))
    b = importlib.util.module_from_spec(spec); spec.loader.exec_module(b)
    m = json.load(open(os.path.join(HERE, "dpyd_notation_map.json")))
    cp = m["cpic_alleles"]
    sg = b.stargazer_codes()
    assert m["source"]["stargazer"]["star_table_sha256"] == b.sha256(os.path.join(b.SG_DIR, "star_table.tsv"))
    assert m["source"]["stargazer"]["snp_table_sha256"] == b.sha256(os.path.join(b.SG_DIR, "snp_table.tsv"))
    d = m["DPYD"]
    for code, e in d["stargazer"].items():
        assert e["cpic"] in cp, code
        want = sorted(f"{p}:{c}" for p, c in sg[code.upper().replace("*S", "*S")]["core"])
        assert sorted(e["stargazer_core"]) == want, code
        if e["method"] == "position+alt":
            assert [(x["pos"], x["change"]) for x in cp[e["cpic"]]] == sorted(sg[code.upper()]["core"]), code
        else:                                   # indel: one variant, same dbSNP id on both sides
            assert len(cp[e["cpic"]]) == 1 and [cp[e["cpic"]][0]["rsid"]] == e["rsid"], code
        for rs in e["rsid"]:
            assert rs in {x["rsid"] for x in cp[e["cpic"]]}, code
    assert set(d["stargazer"]) | set(d["unscorable_stargazer"]) == {c.lower() for c in sg}
    for tag, e in d["aldy"].items():
        v = cp[e["cpic"]]
        assert len(v) == 1 and v[0]["rsid"] == e["rsid"] and tag == "*" + e["rsid"], tag
    # no rs id may name two CPIC single-variant alleles (would be ambiguous: left out of the map)
    rs_count = {}
    for n, v in cp.items():
        if len(v) == 1:
            rs_count[v[0]["rsid"]] = rs_count.get(v[0]["rsid"], 0) + 1
    assert all(("*" + r) in d["aldy"] for r, c in rs_count.items() if c == 1)
    assert not any(("*" + r) in d["aldy"] for r, c in rs_count.items() if c > 1)
    for n in d["cpic_star_names"].values():
        assert n["cpic"] in cp
    # every CPIC name the map points at exists in CPIC's diplotype table, and the committed file equals a rebuild
    allnames = {h.strip() for k in json.load(open(os.path.join(HERE, "cpic_diplotype.json")))["DPYD"] for h in k.split("/")}
    for sec in ("stargazer", "aldy", "cpic_star_names"):
        for k, e in d[sec].items():
            assert e["cpic"] in allnames, (k, e)
    if os.environ.get("PGX_CPIC_VERIFY"):
        fresh = b.build(m["retrieved"])
        assert fresh["DPYD"] == d and fresh["cpic_alleles"] == cp, "dpyd_notation_map.json is out of date with CPIC"
    # a Stargazer code that agrees with the hand-kept synonym map must agree on the allele
    syn = json.load(open(os.path.join(HERE, "allele_synonyms.json")))["DPYD"]
    for canon, aliases in syn.items():
        for a in aliases:
            if a in d["stargazer"]:
                assert any(x.lower() in (canon, *aliases) or canon in R._cpic_allele_key("DPYD", d["stargazer"][a]["cpic"]).lower()
                           for x in [a]) or R._cpic_allele_key("DPYD", d["stargazer"][a]["cpic"]) == canon, a


def test_dpyd_notation_readings_get_the_cpic_phenotype():
    """Stargazer *S codes and Aldy *rs tags score through the map; composite Aldy names and unmapped codes stay unscorable."""
    f = R.cpic_diplotype_phenotype
    assert f("DPYD", "*1/*S20") == f("DPYD", "Reference/c.496A>G") == "Normal Metabolizer"
    assert f("DPYD", "*1/*rs2297595") == "Normal Metabolizer"
    assert f("DPYD", "*S10/*S66") == f("DPYD", "c.2279C>T/c.85T>C (*9A)") == "Intermediate Metabolizer"
    assert f("DPYD", "*S10/*S3") == f("DPYD", "c.1905+1G>A (*2A)/c.85T>C (*9A)")
    assert f("DPYD", "*1/*4") == f("DPYD", "Reference/c.1601G>A (*4)") is not None
    assert f("DPYD", "*9A/*9A+rs1801159+rs2297595") is None      # composite haplotype: not a CPIC allele
    assert f("DPYD", "*1/*S8") is None                            # Stargazer *S8 (3 variants) is no CPIC allele
    assert R._cpic_allele_key("CYP2C19", "*S20") == "*s20"        # the map is DPYD only


def test_legacy_letter_fold_table_matches_pharmvar_snapshot():
    """The fold table is derived from the dated snapshot of PharmVar's rendered Legacy Label column; with
    PGX_PHARMVAR_PAGES=<dir of <GENE>.txt page texts> the snapshot is also compared with freshly rendered pages."""
    import json, pharmvar_legacy_build as PB
    snap = json.load(open(os.path.join(HERE, "pharmvar_legacy_snapshot.json")))
    tab = json.load(open(os.path.join(HERE, "pharmvar_legacy_suballeles.json")))
    der = PB.derive(snap, tab["checked"])
    assert der == tab, "pharmvar_legacy_suballeles.json is not what pharmvar_legacy_build.derive gives for the snapshot"
    assert {d["pharmvar_version"] for d in snap["genes"].values()} == {"6.2.29"}
    assert tab["fold"]["CYP3A5"]["*3B"]["pharmvar_suballele"] == "CYP3A5*3.002" and tab["fold"]["CYP3A5"]["*3B"]["core_allele"] == "*3"
    assert tab["fold"]["CYP2C19"]["*2B"]["core_allele"] == "*2"
    assert "*1A" not in tab["fold"]["CYP2C19"] and tab["redesignated"]["CYP2C19"]["*1A"]["core_allele"] == "*38"
    for g, labels in tab["fold"].items():          # every folded label sits under the core of the same number
        for lab, e in labels.items():
            assert lab.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ") == e["core_allele"], (g, lab)
    assert not {"TPMT", "UGT1A1", "DPYD"} & set(tab["fold"]) and not {"TPMT", "UGT1A1", "DPYD"} & set(snap["genes"])
    pages = os.environ.get("PGX_PHARMVAR_PAGES")
    if pages:
        for g, d in snap["genes"].items():
            assert PB.parse(g, open(os.path.join(pages, g + ".txt"), encoding="utf-8").read()) == d, g


def test_legacy_letter_fold_refuses_unlisted_labels():
    f = R.fold_legacy_letter
    assert f("CYP3A5", "*3B") == "*3" and f("CYP2C19", "*2B") == "*2" and f("CYP2C9", "*3A") == "*3"
    assert f("CYP2C19", "*1A") == "*1A"              # re-designated (now *38.001): not folded to *1
    assert f("CYP2C19", "*2Z") == "*2Z" and f("CYP3A5", "*2B") == "*2B"   # letter not in PharmVar's column for that gene
    assert f("TPMT", "*3A") == "*3A" and f("UGT1A1", "*28A") == "*28A"    # genes PharmVar does not maintain
    assert f("CYP2C19", "*2") == "*2" and f("DPYD", "*9A") == "*9A"
    assert R.cpic_diplotype_phenotype("CYP3A5", "*3/*3B") == R.cpic_diplotype_phenotype("CYP3A5", "*3/*3") == "Poor Metabolizer"
    assert R.cpic_diplotype_phenotype("CYP2C19", "*1/*1A") is None      # ambiguous: stays unscorable
    assert R.cpic_diplotype_phenotype("CYP2C19", "*1/*2Z") is None


def test_split_rule_review_fixes():
    sp = cmp._split_by_cpic_phenotype
    base = {"consensus_diplotype": "*2/*2 or *2/*3", "consensus_phenotype": "Intermediate Metabolizer",
            "status": "concordant", "phenotype_tier": True, "authority": "Phenotype"}
    # all readings one CPIC phenotype that differs from the callers' label: CPIC's phenotype is reported, label in the note
    v = sp(dict(base), "CYP2C19")
    assert v["consensus_phenotype"] == "Poor Metabolizer" and v["status"] == "concordant" and "callers said Intermediate" in v["note"], v
    # same for a phase-ambiguous verdict (no phenotype_tier flag)
    v = sp({"consensus_diplotype": "*2/*2 or *2/*3", "consensus_phenotype": "Normal Metabolizer", "status": "majority",
            "phase_ambiguous": True}, "CYP2C19")
    assert v["consensus_phenotype"] == "Poor Metabolizer" and v["status"] == "majority", v
    # an " or " verdict with neither flag is still checked: a CPIC split withdraws the phenotype, level kept
    v = sp({"consensus_diplotype": "*1/*2 or *2/*2", "consensus_phenotype": "Poor Metabolizer", "status": "majority",
            "authority": "Phenotype", "phenotype_tier": True}, "CYP2C19")
    assert v["status"] == cmp.SPLIT_STATUS_MAJORITY == "majority, diplotype unresolved" and v["consensus_phenotype"] == "-", v
    assert "authority" not in v and v["phenotype_split"] and "phenotype_tier" not in v, v
    v = sp({"consensus_diplotype": "*1/*2 or *2/*2", "consensus_phenotype": "Poor Metabolizer", "status": "concordant"}, "CYP2C19")
    assert v["status"] == cmp.SPLIT_STATUS, v
    # all readings Indeterminate: no phenotype, like a split
    with _stub(cmp, "_reading_phenotype", lambda g, d: (True, "Indeterminate", "scored")):
        v = sp(dict(base), "CYP2C19")
        assert v["consensus_phenotype"] == "-" and v["status"] == cmp.SPLIT_STATUS and "no usable phenotype" in v["note"], v


def test_axiom_adjudicate_split_status():
    import importlib.util
    import os
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "build_axiom_concordance.py")
    sp = importlib.util.spec_from_file_location("bac", root)
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    for st in ("concordant", "majority"):
        base = m.adjudicate("CYP2C19", "*1/*2", st, "*1/*1", "*1/*1", "MISMATCH")
        for suffix in (", diplotype unresolved",):
            assert m.adjudicate("CYP2C19", "*1/*2", st + suffix, "*1/*1", "*1/*1", "MISMATCH") == base, st


def test_reconcile_review_fixes():
    R = cmp.reconcile
    assert R._strip_allele("rs1801265+rs2297595") == "rs1801265" and R._strip_allele("*28+*80") == "*28"
    assert R._strip_allele("c.1905+1G>A") == "c.1905+1g>a"
    assert R.fold_legacy_letter("CYP3A5", "*3b") == "*3" and R.fold_legacy_letter("CYP3A5", "*3B") == "*3"
    # a None (ambiguous) row in the first tables does not stop the folded lookup
    exact, coll = R._dip_tables()["CYP3A5"]
    k = R._dip_key("CYP3A5", "*3/*3B", False)
    saved = exact.get(k, "absent")
    exact[k] = None
    try:
        assert R.cpic_diplotype_phenotype("CYP3A5", "*3/*3B") == "Poor Metabolizer"
    finally:
        if saved == "absent":
            exact.pop(k)
        else:
            exact[k] = saved
    # a missing lookup file fails loud at use, naming the file, not at import
    old = R._LEGACY_PATH, R._legacy_fold
    R._LEGACY_PATH, R._legacy_fold = "/nonexistent/pharmvar_legacy_suballeles.json", None
    try:
        try:
            R.fold_legacy_letter("CYP3A5", "*3B")
            raise AssertionError("missing file did not raise")
        except RuntimeError as e:
            assert "pharmvar_legacy_suballeles.json" in str(e)
    finally:
        R._LEGACY_PATH, R._legacy_fold = old


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
    print("all verdict self-checks passed")
