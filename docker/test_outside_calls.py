#!/usr/bin/env python3
"""Tests for pgx_outside_calls.py (PharmCAT outside-call builder). Pure Python."""
import importlib.util, json, os, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("oc", os.path.join(HERE, "pgx_outside_calls.py"))
oc = importlib.util.module_from_spec(spec); spec.loader.exec_module(oc)

HDR = "Gene\tSample\tBuild\tTool\tDiplotype\tActivityScore\tPhenotype\tStatus\tSVMode\n"


def table(out, gene, rows, base="S_30X"):
    d = os.path.join(out, "genes", gene); os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, f"{gene}_{base}_comparison.tsv"), "w") as fh:
        fh.write(HDR)
        for tool, dip, st in rows:
            fh.write(f"{gene}\t{base}\tGRCh38\t{tool}\t{dip}\t-\t-\t{st}\t-\n")


def passed(out):
    p = os.path.join(out, "pharmcat", "outside_calls.pharmcat.tsv")
    return dict(l.rstrip("\n").split("\t") for l in open(p)) if os.path.exists(p) else {}


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
    print("outside-call tests: PASS")


if __name__ == "__main__":
    main()
