#!/usr/bin/env python3
"""Self-check for the coverage gate + failed-tool marking in pgx-compare.py.

No framework: run directly ->  python3 docker/test_coverage_gate.py
Exercises the real pgx-compare.py CLI against synthesized caller outputs.
"""
import csv
import io
import json
import os
import shutil
import subprocess
import tempfile
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
CMP = os.path.join(HERE, "pgx-compare.py")


def _pypgx_zip(d, dip="*1/*1"):
    """Reproduce PyPGx's results.zip -> data.tsv end-state (its wild-type call)."""
    os.makedirs(os.path.join(d, "pypgx"), exist_ok=True)
    buf = io.StringIO()
    w = csv.writer(buf, delimiter="\t")
    w.writerow(["Name", "Genotype", "Phenotype", "ActivityScore",
                "Haplotype1", "Haplotype2", "VariantData", "CNV", "AlternativePhase"])
    w.writerow(["S", dip, "Normal Metabolizer", "2.0", "*1", "*1", "", "Normal", ""])
    with zipfile.ZipFile(os.path.join(d, "pypgx", "results.zip"), "w") as zf:
        zf.writestr("p/data.tsv", buf.getvalue())


def _run(d, *extra):
    subprocess.run(
        ["python3", CMP, "--gene", "CYP2D6", "--sample", "S", "--output-dir", d, *extra],
        check=True, capture_output=True, text=True,
    )
    with open(os.path.join(d, "CYP2D6_S_comparison.tsv")) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    with open(os.path.join(d, "CYP2D6_S_detail.json")) as fh:
        detail = json.load(fh)
    return {r["Tool"]: r for r in rows}, detail


def test_low_coverage_forces_no_call():
    d = tempfile.mkdtemp()
    try:
        _pypgx_zip(d, "*1/*1")  # caller confidently calls wild-type from no data
        by, detail = _run(d, "--region-depth", "3.0", "--min-depth", "10")
        assert by["PyPGx"]["Diplotype"] == "NO_CALL", by["PyPGx"]
        assert by["PyPGx"]["Status"] == "no_call", by["PyPGx"]
        assert detail["coverage"]["verdict"] == "NO_CALL", detail["coverage"]
    finally:
        shutil.rmtree(d)


def test_adequate_coverage_preserves_call():
    d = tempfile.mkdtemp()
    try:
        _pypgx_zip(d, "*1/*4")
        by, detail = _run(d, "--region-depth", "45.0", "--min-depth", "10")
        assert by["PyPGx"]["Diplotype"] == "*1/*4", by["PyPGx"]
        assert detail["coverage"]["verdict"] == "OK", detail["coverage"]
    finally:
        shutil.rmtree(d)


def test_failed_tool_distinguished_from_not_run():
    d = tempfile.mkdtemp()  # no tool outputs at all
    try:
        by, _ = _run(d, "--failed-tools", "Aldy")
        assert by["Aldy"]["Status"] == "failed", by["Aldy"]
        assert by["Stargazer"]["Status"] == "not_run", by["Stargazer"]
    finally:
        shutil.rmtree(d)


if __name__ == "__main__":
    test_low_coverage_forces_no_call()
    test_adequate_coverage_preserves_call()
    test_failed_tool_distinguished_from_not_run()
    print("all coverage-gate self-checks passed")
