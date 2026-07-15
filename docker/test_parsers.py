#!/usr/bin/env python3
"""Fixture-based unit tests for the pgx-compare.py output parsers (T1).

The parsers translate nine different tool output formats and are the most
format-fragile code in the repo: a tool version bump that changes a column name
or marker silently turns a real call into status=failed (or a wrong call). Each
test writes a representative fixture in the exact layout the parser globs for,
runs the real parser, and asserts the parsed CallerResult.

Fixtures are synthetic but format-accurate, derived from each parser's contract.

NOT covered here: the 3 vcf_check parsers (CACNA1S / G6PD / UGT1A1) shell out to
bcftools on a gene VCF, so they need bcftools + a VCF fixture — an in-container
test (bcftools is not on the host). Tracked as a follow-up.

Run: python3 docker/test_parsers.py
"""
import csv
import importlib.util
import io
import json
import os
import shutil
import tempfile
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


cmp = _load("pgx_compare", "pgx-compare.py")


def test_pypgx():
    d = tempfile.mkdtemp()
    try:
        os.makedirs(os.path.join(d, "pypgx"))
        buf = io.StringIO()
        w = csv.writer(buf, delimiter="\t")
        w.writerow(["Name", "Genotype", "Phenotype", "ActivityScore",
                    "Haplotype1", "Haplotype2", "CNV", "VariantData", "AlternativePhase"])
        w.writerow(["S", "*1/*4", "Intermediate Metabolizer", "1.0", "*1", "*4", "Normal", "", ""])
        with zipfile.ZipFile(os.path.join(d, "pypgx", "results.zip"), "w") as z:
            z.writestr("CYP2D6-pipeline/data.tsv", buf.getvalue())
        r = cmp.parse_pypgx(d, "CYP2D6", "S")
        assert r.status == "ok" and r.diplotype == "*1/*4", vars(r)
        assert r.activity_score == "1.0" and r.phenotype == "Intermediate Metabolizer", vars(r)
    finally:
        shutil.rmtree(d)


def test_stargazer():
    d = tempfile.mkdtemp()
    try:
        os.makedirs(os.path.join(d, "stargazer"))
        with open(os.path.join(d, "stargazer", "genotype-calls.tsv"), "w") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["hap1_main", "hap2_main", "dip_score", "phenotype"])
            w.writerow(["*1", "*4", "1.0", "intermediate metabolizer"])
        r = cmp.parse_stargazer(d, "CYP2D6", "S")
        assert r.status == "ok" and r.diplotype == "*1/*4" and r.activity_score == "1.0", vars(r)


    finally:
        shutil.rmtree(d)


def test_stargazer_negative_dip_score_is_not_run():
    # Stargazer emits a negative dip_score when it cannot make a confident call.
    d = tempfile.mkdtemp()
    try:
        os.makedirs(os.path.join(d, "stargazer"))
        with open(os.path.join(d, "stargazer", "genotype-calls.tsv"), "w") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["hap1_main", "hap2_main", "dip_score", "phenotype"])
            w.writerow(["*1", "*4", "-1.0", "x"])
        r = cmp.parse_stargazer(d, "CYP2D6", "S")
        assert r.status == "not_run", vars(r)
    finally:
        shutil.rmtree(d)


def test_aldy():
    d = tempfile.mkdtemp()
    try:
        os.makedirs(os.path.join(d, "aldy"))
        lines = [
            "#Sample\tGene\tSolutionID\tMajor\tMinor\tCopy\tAllele\tLocation\tType\tCoverage\tEffect\tdbSNP",
            "#Solution 1: cpic=intermediate_metabolizer;cpic_score=1.0",
            "S\tCYP2D6\t1\t*1/*4\t*1/*4\t1\t\t\t\t\t\t",
        ]
        with open(os.path.join(d, "aldy", "CYP2D6.aldy"), "w") as fh:
            fh.write("\n".join(lines) + "\n")
        r = cmp.parse_aldy(d, "CYP2D6", "S")
        assert r.status == "ok" and r.diplotype == "*1/*4", vars(r)
        assert r.activity_score == "1.0" and r.phenotype == "Intermediate Metabolizer", vars(r)
    finally:
        shutil.rmtree(d)


def test_stellarpgx():
    d = tempfile.mkdtemp()
    try:
        ad = os.path.join(d, "stellarpgx", "cyp2d6", "alleles")
        os.makedirs(ad)
        txt = ("Initially computed CN = 2\n"
               "Result:\n*1/*4\n"
               "Activity score:\n1.0\n"
               "Metaboliser status:\nNormal metaboliser\n")
        with open(os.path.join(ad, "res.alleles"), "w") as fh:
            fh.write(txt)
        r = cmp.parse_stellarpgx(d, "CYP2D6", "S")
        assert r.status == "ok" and r.diplotype == "*1/*4", vars(r)
        assert r.activity_score == "1.0" and r.phenotype == "Normal metaboliser", vars(r)
        assert r.copy_number == "2", vars(r)
    finally:
        shutil.rmtree(d)


def test_optitype_hla_a_b_and_c():
    d = tempfile.mkdtemp()
    try:
        os.makedirs(os.path.join(d, "optitype"))
        with open(os.path.join(d, "optitype", "S_result.tsv"), "w") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["", "A1", "A2", "B1", "B2", "C1", "C2", "Reads", "Objective"])
            w.writerow(["0", "A*01:01", "A*02:01", "B*57:01", "B*07:02", "C*01:02", "C*02:02", "1234", "56.7"])
        ra = cmp.parse_optitype(d, "HLA-A", "S")
        assert ra.status == "ok" and ra.diplotype == "HLA-A*01:01/HLA-A*02:01", vars(ra)
        rb = cmp.parse_optitype(d, "HLA-B", "S")
        assert rb.status == "ok" and rb.diplotype == "HLA-B*57:01/HLA-B*07:02", vars(rb)
        # HLA-C: OptiType always emitted C1/C2, but the parser used to discard them,
        # which left the downstream HLA-C*06:02 psoriasis association unscreenable.
        rc = cmp.parse_optitype(d, "HLA-C", "S")
        assert rc.status == "ok" and rc.diplotype == "HLA-C*01:02/HLA-C*02:02", vars(rc)
        assert rc.haplotype1 == "HLA-C*01:02" and rc.haplotype2 == "HLA-C*02:02", vars(rc)
    finally:
        shutil.rmtree(d)


def test_optitype_hla_c_homozygous_psoriasis_allele():
    # Regression: a C*06:02 homozygote must read as a real call, not "absent".
    d = tempfile.mkdtemp()
    try:
        os.makedirs(os.path.join(d, "optitype"))
        with open(os.path.join(d, "optitype", "S_result.tsv"), "w") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["", "A1", "A2", "B1", "B2", "C1", "C2", "Reads", "Objective"])
            w.writerow(["0", "A*30:01", "A*30:01", "B*13:01", "B*38:01", "C*06:02", "C*06:02", "361.0", "347.99"])
        rc = cmp.parse_optitype(d, "HLA-C", "S")
        assert rc.status == "ok" and rc.diplotype == "HLA-C*06:02/HLA-C*06:02", vars(rc)
    finally:
        shutil.rmtree(d)


def test_mutserve():
    d = tempfile.mkdtemp()
    try:
        os.makedirs(os.path.join(d, "mt-rnr1"))
        payload = {
            "diplotype": "m.1555A>G",
            "phenotype": "Increased risk of aminoglycoside-induced hearing loss",
            "variants": [{"label": "m.1555A>G", "af": 0.99, "type": "homoplasmic"}],
        }
        with open(os.path.join(d, "mt-rnr1", "S_mtrna1_result.json"), "w") as fh:
            json.dump(payload, fh)
        r = cmp.parse_mutserve(d, "MT-RNR1", "S")
        assert r.status == "ok" and r.diplotype == "m.1555A>G", vars(r)
    finally:
        shutil.rmtree(d)


def test_missing_output_is_not_run():
    # An empty output dir (caller produced nothing) must read as not_run, NOT a call.
    d = tempfile.mkdtemp()
    try:
        for fn in (cmp.parse_pypgx, cmp.parse_stargazer, cmp.parse_aldy,
                   cmp.parse_stellarpgx, cmp.parse_optitype, cmp.parse_mutserve):
            r = fn(d, "CYP2D6", "S")
            assert r.status == "not_run", (fn.__name__, vars(r))
    finally:
        shutil.rmtree(d)


if __name__ == "__main__":
    for _name, _fn in sorted(globals().items()):
        if _name.startswith("test_") and callable(_fn):
            _fn()
    print("all parser self-checks passed")
