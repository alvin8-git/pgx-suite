#!/usr/bin/env python3
"""Tests for the Cyrius CYP2D6 caller integration (design C).

parse_cyrius reads Cyrius's Sample/Genotype/Filter TSV. This checks genotype
extraction (incl. CNV xN and hybrids), no-call handling, graceful missing
output, and that CYP2D6 is wired to cyrius in genes.tsv. Pure-Python: Cyrius
itself is not run here (it needs a BAM); that path is covered by the in-image
smoke run.
"""
import importlib.util
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def _load_compare():
    spec = importlib.util.spec_from_file_location("c", os.path.join(HERE, "pgx-compare.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _write_cyrius_tsv(gene_dir, sample, genotype, filt="PASS"):
    cdir = os.path.join(gene_dir, "cyrius")
    os.makedirs(cdir, exist_ok=True)
    with open(os.path.join(cdir, f"{sample}.tsv"), "w") as fh:
        fh.write("Sample\tGenotype\tFilter\n")
        fh.write(f"{sample}\t{genotype}\t{filt}\n")


def main():
    c = _load_compare()

    # CYP2D6 is wired to cyrius
    support, _, _ = c.load_gene_config()
    assert support["CYP2D6"].get("cyrius") is True, "CYP2D6 must enable cyrius in genes.tsv"
    assert not any(support[g].get("cyrius") for g in support if g != "CYP2D6"), \
        "only CYP2D6 should use cyrius"

    with tempfile.TemporaryDirectory() as d:
        # straightforward diplotype
        _write_cyrius_tsv(d, "S1", "*1/*4")
        r = c.parse_cyrius(d, "CYP2D6", "S1")
        assert r.status == "ok" and r.diplotype == "*1/*4" and r.tool == "Cyrius"
        assert r.haplotype1 == "*1" and r.haplotype2 == "*4"

    with tempfile.TemporaryDirectory() as d:
        # CNV / hybrid genotype with a non-PASS filter preserved as alt note
        _write_cyrius_tsv(d, "S2", "*68+*4/*10xN", filt="LowQ")
        r = c.parse_cyrius(d, "CYP2D6", "S2")
        assert r.status == "ok" and r.diplotype == "*68+*4/*10xN"
        assert r.alternative_diplotypes == "LowQ"

    with tempfile.TemporaryDirectory() as d:
        # explicit no-call
        _write_cyrius_tsv(d, "S3", "None", filt="LowDepth")
        r = c.parse_cyrius(d, "CYP2D6", "S3")
        assert r.status == "failed" and "no-call" in r.diplotype

    with tempfile.TemporaryDirectory() as d:
        # missing output -> graceful failure, not a crash
        r = c.parse_cyrius(d, "CYP2D6", "S4")
        assert r.status == "failed" and "no Cyrius output" in r.diplotype

    print("all Cyrius integration checks passed")


if __name__ == "__main__":
    try:
        main()
    except AssertionError as e:
        print("FAIL:", e, file=sys.stderr)
        sys.exit(1)
