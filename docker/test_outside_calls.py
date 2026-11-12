#!/usr/bin/env python3
"""Tests for pgx_outside_calls.py (PharmCAT outside-call builder). Pure Python."""
import importlib.util, json, os, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("oc", os.path.join(HERE, "pgx_outside_calls.py"))
oc = importlib.util.module_from_spec(spec); spec.loader.exec_module(oc)

HDR = "Gene\tSample\tBuild\tTool\tDiplotype\tActivityScore\tPhenotype\tStatus\tSVMode\n"


def table(out, gene, rows, base="S_30X"):
    """rows: (tool, diplotype, status[, phenotype])"""
    d = os.path.join(out, "genes", gene); os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, f"{gene}_{base}_comparison.tsv"), "w") as fh:
        fh.write(HDR)
        for tool, dip, st, *ph in rows:
            fh.write(f"{gene}\t{base}\tGRCh38\t{tool}\t{dip}\t-\t{ph[0] if ph else '-'}\t{st}\t-\n")


def passed(out):
    """gene -> diplotype column (first two columns of the PharmCAT file)."""
    p = os.path.join(out, "pharmcat", "outside_calls.pharmcat.tsv")
    return {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p)} if os.path.exists(p) else {}


def _cols(out, gene="CYP2D6"):
    """(diplotype, phenotype, score) strings of the PharmCAT line, or None when omitted."""
    p = os.path.join(out, "pharmcat", "outside_calls.pharmcat.tsv")
    for l in open(p):
        c = l.rstrip("\n").split("\t")
        if c[0] == gene:
            return tuple((c + ["", "", ""])[1:4])
    return None


def _decide(rows):
    """-> (pharmcat columns or None, record row dict) for CYP2D6 rows (tool, dip[, pheno])."""
    with tempfile.TemporaryDirectory() as out:
        table(out, "CYP2D6", [(r[0], r[1], "ok", *r[2:]) for r in rows])
        recs = oc.build(out)
        return _cols(out), recs[-1]


def _cyp2d6(rows):
    cols, rec = _decide(rows)
    return (cols[0] if cols else None), rec["why"] or rec["level"]


def test_notation_only_differences():
    # Aldy tandem X marker; tandem order; sub-allele suffix; legacy letter sub-allele
    d, _ = _cyp2d6([("PyPGx", "*2/*36+*10"), ("Stargazer", "*2/*36+*10"), ("Aldy", "*2/*36X+*10"),
                    ("StellarPGx", "*36+*10/*2"), ("Cyrius", "*10+*36/*2")])
    assert d == "*10+*36/*2", d          # Cyrius spelling is the one passed (cleaned)
    d, _ = _cyp2d6([("PyPGx", "*2/*36+*10"), ("Aldy", "*2/*36X+*10")])
    assert d == "*2/*36+*10", d          # Aldy's X marker never reaches PharmCAT
    d, _ = _cyp2d6([("PyPGx", "*3/*4"), ("Aldy", "*3/*4.021")])
    assert d == "*3/*4", d
    d, _ = _cyp2d6([("PyPGx", "*6/*35"), ("Aldy", "*6B/*35")])
    assert d == "*6/*35", d
    # rs tags are a real difference (never stripped) and give no score; with only two callers
    # and no phenotype for Aldy nothing can be agreed
    d, why = _cyp2d6([("PyPGx", "*2/*41"), ("Aldy", "*2+rs61737946+rs61745683/*41")])
    assert d is None and "no CPIC activity score" in why, why
    # copy number is a real difference: Normal (2.0) vs Ultrarapid (3.0) -> omitted
    d, why = _cyp2d6([("PyPGx", "*2/*2"), ("Aldy", "*2/*2x2")])
    assert d is None and "disagree" in why, why


def test_activity_score_and_phenotype():
    import reconcile as rc
    cases = {"*2/*4": 1.0, "*2/*36+*10": 1.25, "*36+*10/*2": 1.25, "*10+*36/*2": 1.25,
             "*2/*36X+*10": 1.25, "*3/*4.021": 0.0, "*3/*68+*4": 0.0, "*10/*10": 0.5,
             "*36+*10/*49": 0.75, "*2/*2x2": 3.0, "*4/*4x2": 0.0, "*1/*1": 2.0}
    for dip, sc in cases.items():
        assert rc.activity_score(dip) == sc, (dip, rc.activity_score(dip))
    # unknown: allele not in the CPIC table (*79), rs tag, unknown copy number, no `/`
    for dip in ("*2/*79", "*2+rs61737946/*41", "*2x\u22653/*1", "*2", "*2/*999"):
        assert rc.activity_score(dip) is None, dip
    assert [rc.score_phenotype(x) for x in (0.0, 0.25, 1.0, 1.25, 2.25, 2.5, 3.0, 3.25, 4.5, 6.5, None)] == [
        "Poor Metabolizer", "Intermediate Metabolizer", "Intermediate Metabolizer", "Normal Metabolizer",
        "Normal Metabolizer", "Ultrarapid Metabolizer", "Ultrarapid Metabolizer", "Ultrarapid Metabolizer",
        "Ultrarapid Metabolizer", "Ultrarapid Metabolizer", None]
    assert rc.phenotype_class("Intermediate Metabolizer") == "Intermediate"
    assert rc.phenotype_class("Indeterminate") == "Indeterminate" and rc.phenotype_class("-") is None


def test_giab_levels():
    five = ("PyPGx", "Stargazer", "Aldy", "StellarPGx", "Cyrius")
    # HG002: five identical -> (a), score 1.0
    cols, rec = _decide([(t, "*2/*4") for t in five])
    assert rec["level"] == "a" and cols == ("*2/*4", "", "1.0") and rec["source"] == "all callers agree", (cols, rec)
    # *2/*36+*10 x4, Aldy `*36X`: (a) with score 1.25
    cols, rec = _decide([(t, "*2/*36+*10") for t in ("PyPGx", "Stargazer", "StellarPGx", "Cyrius")]
                        + [("Aldy", "*2/*36X+*10")])
    assert rec["level"] == "a" and cols == ("*2/*36+*10", "", "1.25"), (cols, rec)
    # HG001: Cyrius `*68+*4` -> diplotypes differ, all scores 0.0 -> (b)
    cols, rec = _decide([("PyPGx", "*3/*4"), ("Stargazer", "*3/*4"), ("Aldy", "*3/*4.021"),
                         ("StellarPGx", "*3/*4"), ("Cyrius", "*3/*68+*4")])
    assert rec["level"] == "b" and cols == ("", "", "0.0") and rec["source"] == "activity score agreed 5/5", (cols, rec)
    # HG005: scores 0.5 and 0.75 differ; Cyrius (0.75) is matched by Aldy and StellarPGx -> (e), the
    # Cyrius diplotype with 0.75 (before the Cyrius rule this was (c), phenotype only)
    cols, rec = _decide([("PyPGx", "*10/*36+*10"), ("Stargazer", "*10/*36+*10"), ("Aldy", "*36X+*10/*49"),
                         ("StellarPGx", "*49/*36+*10"), ("Cyrius", "*49/*36+*10")])
    assert rec["level"] == "e" and cols == ("*49/*36+*10", "", "0.75") \
        and rec["source"] == "Cyrius + 2 other caller(s) on activity score 0.75" \
        and "PyPGx" not in rec["why"].split("dissenting")[0] and "dissenting: PyPGx=*10/*36+*10" in rec["why"] \
        and "Stargazer=*10/*36+*10" in rec["why"], (cols, rec)
    # level b keeps precedence over e: five scores equal -> score only (the HG001 control above)
    # level c still applies when Cyrius gave no call, or no one shares its score
    cols, rec = _decide([("PyPGx", "*10/*36+*10"), ("Stargazer", "*10/*36+*10"), ("Aldy", "*36X+*10/*49"),
                         ("StellarPGx", "*49/*36+*10")])
    assert rec["level"] == "c" and cols == ("", "Intermediate Metabolizer", ""), (cols, rec)
    cols, rec = _decide([("PyPGx", "*1/*10"), ("Aldy", "*1/*10"), ("Cyrius", "*1/*1", "Normal Metabolizer"),
                         ("StellarPGx", "*1/*2")])
    assert rec["level"] == "e", (cols, rec)      # Cyrius 2.0 == StellarPGx 2.0; PyPGx/Aldy 1.25 dissent
    cols, rec = _decide([("PyPGx", "*1/*10"), ("Aldy", "*1/*10"), ("Cyrius", "*1/*4")])
    assert cols is None, (cols, rec)             # Cyrius 1.0 shared by nobody
    # HG006: 0.5 vs 0.5 (CPIC *36+*10 = 0.25) -> (b)
    cols, rec = _decide([(t, "*36+*10/*36+*10") for t in ("PyPGx", "Stargazer", "StellarPGx", "Cyrius")]
                        + [("Aldy", "*10/*36X+*10")])
    assert rec["level"] == "b" and cols == ("", "", "0.5"), (cols, rec)
    # HG007: Aldy `*10/*10` = 0.5 == `*10/*36+*10` = 0.5 -> (b)
    cols, rec = _decide([(t, "*10/*36+*10") for t in ("PyPGx", "Stargazer", "StellarPGx", "Cyrius")]
                        + [("Aldy", "*10/*10")])
    assert rec["level"] == "b" and cols == ("", "", "0.5"), (cols, rec)
    # HG004: *79 is not in the CPIC table (score unknown); Normal (2.0) vs Ultrarapid (3.0) -> (d)
    # (d) still holds without Cyrius; with Cyrius's 3.0 shared by Stargazer/StellarPGx it is (e)
    cols, rec = _decide([("PyPGx", "*2/*2"), ("Stargazer", "*2/*2x2"), ("Aldy", "*2/*79"),
                         ("StellarPGx", "*2/*2x2")])
    assert cols is None and rec["level"] == "-" and "disagree" in rec["why"], (cols, rec)
    cols, rec = _decide([("PyPGx", "*2/*2"), ("Stargazer", "*2/*2x2"), ("Aldy", "*2/*79"),
                         ("StellarPGx", "*2/*2x2"), ("Cyrius", "*2/*2x2")])
    assert rec["level"] == "e" and cols == ("*2/*2x2", "", "3.0") and "PyPGx=*2/*2" in rec["why"] \
        and "Aldy=*2/*79" in rec["why"], (cols, rec)
    # four agree on *2/*41; the fifth has +rsNNN tags and reports Indeterminate: not (a) (rs tags
    # differ), not (b) (its score is unknown), and (c) is blocked by the explicit Indeterminate
    # -> omitted
    rows = [(t, "*2/*41") for t in ("PyPGx", "Stargazer", "StellarPGx", "Cyrius")]
    rows.append(("Aldy", "*2+rs61737946+rs61745683/*41", "Indeterminate"))
    cols, rec = _decide([r for r in rows if r[0] != "Cyrius"])
    assert cols is None and "Indeterminate" in rec["why"], (cols, rec)
    # with Cyrius present the three callers on *2/*41 are enough for (e); the Indeterminate caller does not
    # block it but is listed
    cols, rec = _decide(rows)
    assert rec["level"] == "e" and cols[0] == "*2/*41" and "Aldy=" in rec["why"] and "Indeterminate" in rec["why"], (cols, rec)
    # same, but the fifth caller has no phenotype at all: its call has no CPIC activity score, and no caller
    # is left out of the comparison -> omitted (review 2026-11-12: a caller is never silently ignored)
    rows[-1] = ("Aldy", "*2+rs61737946+rs61745683/*41")
    cols, rec = _decide([r for r in rows if r[0] != "Cyrius"])
    assert cols is None and "no CPIC activity score" in rec["why"], (cols, rec)
    # an allele CPIC does not score (*79): the caller's own "Normal" label is never used -> omitted
    cols, rec = _decide([("PyPGx", "*2/*2"), ("Stargazer", "*2/*79", "Normal Metabolizer"),
                         ("Aldy", "*1/*2")])
    assert cols is None and "no CPIC activity score" in rec["why"], (cols, rec)
    # copy number written with a capital X scores like x2
    assert oc.reconcile.activity_score("*2X2/*1") == oc.reconcile.activity_score("*2x2/*1") == 3.0
    for bad in ("Not Normal", "Normal to Ultrarapid", "non-poor"):
        assert oc.reconcile.phenotype_class(bad) is None, bad
    assert oc.reconcile.phenotype_class("unknown") == "Indeterminate"
    assert oc.reconcile.phenotype_class("Poor metaboliser (PM)") == "Poor" and oc.reconcile.phenotype_class("normal_metabolizer") == "Normal"


def test_cyrius_activity_score_level():
    """Option B (2026-11-13): Cyrius's call passes when its CPIC score is shared by >=1 other caller."""
    five = ("PyPGx", "Stargazer", "Aldy", "StellarPGx", "Cyrius")
    # PyPGx *1/*41 (1.25) dissents; Stargazer/Aldy *9/*41 and Cyrius *41/*9 are 0.5 -> (e); before: blocked
    cols, rec = _decide([("PyPGx", "*1/*41"), ("Stargazer", "*9/*41"), ("Aldy", "*9/*41"),
                         ("StellarPGx", "*41/*9"), ("Cyrius", "*41/*9")])
    assert oc.reconcile.activity_score("*1/*41") == 1.25 and oc.reconcile.activity_score("*9/*41") == 0.5
    assert rec["level"] == "e" and cols == ("*41/*9", "", "0.5"), (cols, rec)
    assert "dissenting: PyPGx=*1/*41" in rec["why"], rec
    assert rec["source"] == "Cyrius + 3 other caller(s) on activity score 0.5", rec
    # Cyrius no call -> unchanged from before (b/c/d logic only)
    cols, rec = _decide([("PyPGx", "*1/*41"), ("Stargazer", "*9/*41"), ("Aldy", "*9/*41"),
                         ("StellarPGx", "*41/*9")])
    assert rec["level"] != "e", rec
    # Cyrius call without a CPIC score (*134 unscored) -> never e
    cols, rec = _decide([("PyPGx", "*1/*9"), ("Stargazer", "*1/*9"), ("Aldy", "*1/*9"),
                         ("StellarPGx", "*1/*9"), ("Cyrius", "*134/*9")])
    assert oc.reconcile.activity_score("*134/*9") is None
    assert rec["level"] != "e" and cols is None, (cols, rec)
    cols, rec = _decide([("PyPGx", "*1/*10"), ("Stargazer", "*1/*10"), ("Aldy", "*1/*10"),
                         ("StellarPGx", "*1/*36+*10"), ("Cyrius", "*1/*36+*36+*10")])
    assert oc.reconcile.activity_score("*1/*36+*36+*10") is None
    assert rec["level"] != "e", (cols, rec)
    # provenance carries the audit note
    with tempfile.TemporaryDirectory() as out:
        table(out, "CYP2D6", [(t, d, "ok") for t, d in (("PyPGx", "*2/*2"), ("Aldy", "*2/*79"),
                              ("Stargazer", "*2/*2x2"), ("Cyrius", "*2/*2x2"))])
        recs = oc.build(out)
        json.dump({"sample": "S"}, open(os.path.join(out, "provenance.json"), "w"))
        oc.update_provenance(out, recs)
        c = {r["gene"]: r for r in json.load(open(os.path.join(out, "provenance.json")))["pharmcat_outside_calls"]}["CYP2D6"]
        assert c["passed"] and c["level"] == "e" and c["activity_score"] == 3.0 and "dissenting" in c["note"], c


def test_record_and_pharmcat_columns():
    with tempfile.TemporaryDirectory() as out:
        table(out, "CYP2D6", [(t, d, "ok") for t, d in
                              (("PyPGx", "*3/*4"), ("Aldy", "*3/*4.021"), ("Cyrius", "*3/*68+*4"))])
        recs = oc.build(out)
        line = open(os.path.join(out, "pharmcat", "outside_calls.pharmcat.tsv")).read()
        assert line == "CYP2D6\t\t\t0.0\n", repr(line)      # empty diplotype/phenotype columns
        rec = open(os.path.join(out, "pharmcat", "outside_calls.tsv")).read().splitlines()
        assert rec[0].split("\t")[:6] == ["gene", "diplotype", "activity_score", "phenotype", "level", "source"]
        row = rec[-1].split("\t")
        assert row[0] == "CYP2D6" and row[4] == "b" and "Cyrius=*3/*68+*4|0.0|Poor" in row[-1], row
        json.dump({"sample": "S"}, open(os.path.join(out, "provenance.json"), "w"))
        oc.update_provenance(out, recs)
        g = {r["gene"]: r for r in json.load(open(os.path.join(out, "provenance.json")))["pharmcat_outside_calls"]}
        c = g["CYP2D6"]
        assert c["passed"] and c["level"] == "b" and c["activity_score"] == 0.0 and c["diplotype"] is None
        assert [x["tool"] for x in c["per_caller"]] == ["PyPGx", "Aldy", "Cyrius"]


def main():
    with tempfile.TemporaryDirectory() as out:
        # HLA four-field input is cut to two fields; basename is not assumed
        table(out, "HLA-B", [("OptiType", "HLA-B*15:02:01/HLA-B*55:01", "ok")])
        table(out, "HLA-A", [("OptiType", "-", "failed")])
        # CYP2D6: unanimity ignores order/case/space; one dissenter blocks it
        table(out, "CYP2D6", [("PyPGx", "*4/*2", "ok"), ("Stargazer", "*2 / *4", "ok"),
                              ("Aldy", "-", "failed")])
        oc.build(out)
        p = passed(out)
        assert p == {"HLA-B": "*15:02/*55:01", "CYP2D6": "*4/*2"}, p
        rec = open(os.path.join(out, "pharmcat", "outside_calls.tsv")).read()
        assert "HLA-A\t-\t" in rec and "no\t" in rec and "OptiType" in rec

        table(out, "CYP2D6", [("PyPGx", "*1/*2", "ok"), ("Aldy", "*1/*4", "ok")])
        oc.build(out)
        assert "CYP2D6" not in passed(out), "disagreement must omit CYP2D6"

        table(out, "CYP2D6", [("PyPGx", "*1/*2", "ok"), ("Aldy", "*1/*2", "failed")])
        oc.build(out)
        assert "CYP2D6" not in passed(out), "single typing caller must omit CYP2D6"

        # missing inputs: nothing fabricated, nothing raised
        with tempfile.TemporaryDirectory() as empty:
            oc.build(empty)
            assert passed(empty) == {}

        # provenance gets the record
        json.dump({"sample": "S"}, open(os.path.join(out, "provenance.json"), "w"))
        oc.update_provenance(out, oc.build(out))
        prov = json.load(open(os.path.join(out, "provenance.json")))
        g = {r["gene"]: r for r in prov["pharmcat_outside_calls"]}
        assert g["HLA-B"]["passed"] and g["HLA-B"]["source_tool"] == "OptiType"
        assert not g["CYP2D6"]["passed"] and "reason" in g["CYP2D6"]
    test_notation_only_differences()
    test_activity_score_and_phenotype()
    test_giab_levels()
    test_cyrius_activity_score_level()
    test_record_and_pharmcat_columns()
    print("outside-call tests: PASS")


if __name__ == "__main__":
    main()
