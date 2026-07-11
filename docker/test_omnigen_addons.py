#!/usr/bin/env python3
"""Fixture-based unit tests for omnigen_addons.py — the deterministic logic of
the OmniGen add-on features (HLA class II, ABO). These run WITHOUT a BAM,
container, or full pipeline run. (mtDNA + Y haplogroups moved to OmniGen.)

Run: python3 docker/test_omnigen_addons.py
"""
import importlib.util
import os
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


oa = _load("omnigen_addons", "omnigen_addons.py")


# ── 1. HLA class II — arcasHLA parse + disease-association logic ─────────────
def test_arcashla_parse_het_and_hom():
    geno = {"A": ["A*01:01:01", "A*02:01:01"],
            "DQB1": ["DQB1*06:02:01", "DQB1*03:01:01"],
            "DQA1": ["DQA1*01:02:01"]}  # homozygous single entry
    h = oa.parse_arcashla_genotype(geno, "HLA-DQB1")
    assert h == ("HLA-DQB1*06:02:01", "HLA-DQB1*03:01:01"), h
    hom = oa.parse_arcashla_genotype(geno, "HLA-DQA1")
    assert hom == ("HLA-DQA1*01:02:01", "HLA-DQA1*01:02:01"), hom
    assert oa.parse_arcashla_genotype(geno, "HLA-DRB1") is None


def test_hla2_celiac_dq25():
    calls = {"HLA-DQA1": "HLA-DQA1*05:01/HLA-DQA1*01:02",
             "HLA-DQB1": "HLA-DQB1*02:01/HLA-DQB1*03:01"}
    f = oa.hla2_disease_findings(calls)
    conds = [x["condition"] for x in f]
    assert any("Celiac" in c for c in conds), f


def test_hla2_narcolepsy():
    calls = {"HLA-DQB1": "HLA-DQB1*06:02/HLA-DQB1*03:01"}
    f = oa.hla2_disease_findings(calls)
    assert any("Narcolepsy" in x["condition"] for x in f), f


def test_hla2_no_risk():
    calls = {"HLA-DQA1": "HLA-DQA1*01:02/HLA-DQA1*01:02",
             "HLA-DQB1": "HLA-DQB1*05:01/HLA-DQB1*05:01"}
    assert oa.hla2_disease_findings(calls) == [], oa.hla2_disease_findings(calls)


# ── 2. ABO — assignment logic (coordinate-independent synthetic sites) ───────
_DEFS = {"o_defining_rsid": "rs8176719",
         "ab_discriminating_rsids": ["rs8176746", "rs8176747"]}


def _obs(alt_copies, depth=40, is_indel=False):
    return {"alt_copies": alt_copies, "gt": {0: "0/0", 1: "0/1", 2: "1/1"}[alt_copies],
            "depth": depth, "is_indel": is_indel}


def test_abo_OO_homozygous_deletion():
    sites = {"rs8176719": _obs(2, is_indel=True), "rs8176746": None, "rs8176747": None}
    r = oa.abo_assign_type(sites, _DEFS)
    assert r["ABO_type"] == "O", r
    assert r["alleles"] == "ABO*O.01.01/ABO*O.01.01", r
    assert r["confidence"] == "low" and r["validation_status"].startswith("UNVALIDATED"), r


def test_abo_AO_het_deletion_no_B_snp():
    sites = {"rs8176719": _obs(1, is_indel=True), "rs8176746": None, "rs8176747": None}
    r = oa.abo_assign_type(sites, _DEFS)
    assert r["ABO_type"] == "A", r  # one O + one A
    assert set(r["alleles"].split("/")) == {"ABO*A1.01", "ABO*O.01.01"}, r


def test_abo_AA_reference_no_variants():
    # Homozygous-reference at all 3 sites -> two non-O, no B SNP -> A/A.
    sites = {"rs8176719": None, "rs8176746": None, "rs8176747": None}
    r = oa.abo_assign_type(sites, _DEFS)
    assert r["ABO_type"] == "A" and r["alleles"] == "ABO*A1.01/ABO*A1.01", r


def test_abo_AB_het_B_snps_no_deletion():
    sites = {"rs8176719": None, "rs8176746": _obs(1), "rs8176747": _obs(1)}
    r = oa.abo_assign_type(sites, _DEFS)
    assert r["ABO_type"] == "AB", r
    assert set(r["alleles"].split("/")) == {"ABO*A1.01", "ABO*B.01"}, r


def test_abo_BB_homozygous_B():
    sites = {"rs8176719": None, "rs8176746": _obs(2), "rs8176747": _obs(2)}
    r = oa.abo_assign_type(sites, _DEFS)
    assert r["ABO_type"] == "B" and r["alleles"] == "ABO*B.01/ABO*B.01", r


def test_abo_BO_het_deletion_plus_B_snp():
    sites = {"rs8176719": _obs(1, is_indel=True), "rs8176746": _obs(1), "rs8176747": _obs(1)}
    r = oa.abo_assign_type(sites, _DEFS)
    assert r["ABO_type"] == "B", r  # one O + one B
    assert set(r["alleles"].split("/")) == {"ABO*B.01", "ABO*O.01.01"}, r


def test_abo_tsv_round_trip():
    d = tempfile.mkdtemp()
    p = os.path.join(d, "abo", "abo_type.tsv")
    oa.write_abo_tsv(p, "HG002", {"ABO_type": "A", "alleles": "ABO*A1.01/ABO*O.01.01",
                                  "confidence": "low"})
    row = oa.read_abo_tsv(p)
    assert row["ABO_type"] == "A" and row["confidence"] == "low", row


def test_abo_defs_file_loads_and_is_flagged_unvalidated():
    defs = oa.load_abo_defs()
    assert defs["o_defining_rsid"] == "rs8176719", defs
    assert "UNVALIDATED" in defs["validation_status"], defs
    assert len(defs["ab_discriminating_rsids"]) == 2, defs


if __name__ == "__main__":
    for _name, _fn in sorted(globals().items()):
        if _name.startswith("test_") and callable(_fn):
            _fn()
    print("all omnigen_addons self-checks passed")
