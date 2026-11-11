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


# ── 2. ABO — read-level assignment on synthetic pileups ─────────────────────
_DEFS = oa.load_abo_defs()
_B = {"rs7853989": ("G", "C"), "rs8176743": ("C", "T"),
      "rs8176746": ("G", "T"), "rs8176747": ("C", "G")}


def _pile(non_o_frac, b_frac, depth=40, cover_b=True):
    """Synthetic pileup: non_o_frac of reads carry the C insertion at c.261;
    b_frac of reads carry the B base at the four A/B SNPs."""
    sites = {}
    n_ins = round(depth * non_o_frac)
    sites["rs8176719"] = {"depth": depth, "ref": "T", "ref_count": depth,
                          "alt_counts": {}, "ins_count": n_ins,
                          "ins_seqs": {"C": n_ins} if n_ins else {}}
    for rs, (ref, b) in _B.items():
        if not cover_b:
            sites[rs] = {"depth": 0, "ref": ref, "ref_count": 0, "alt_counts": {}}
            continue
        nb = round(depth * b_frac)
        sites[rs] = {"depth": depth, "ref": ref, "ref_count": depth - nb,
                     "alt_counts": {b: nb} if nb else {}}
    return {"sites": sites}


def test_abo_hom_ref_is_OO_not_AA():
    r = oa.abo_assign_type(_pile(0.0, 0.0), _DEFS)
    assert r["ABO_type"] == "O" and r["alleles"] == "ABO*O/ABO*O", r
    assert r["evidence"]["o_site_depth"] == 40, r


def test_abo_AO_het_insertion_no_B():
    r = oa.abo_assign_type(_pile(0.5, 0.0), _DEFS)
    assert r["ABO_type"] == "A" and r["alleles"] == "ABO*A/ABO*O", r


def test_abo_AA_hom_insertion_no_B():
    r = oa.abo_assign_type(_pile(1.0, 0.0), _DEFS)
    assert r["ABO_type"] == "A" and r["alleles"] == "ABO*A/ABO*A", r


def test_abo_AB_het_B_snps_hom_insertion():
    r = oa.abo_assign_type(_pile(1.0, 0.5), _DEFS)
    assert r["ABO_type"] == "AB" and r["alleles"] == "ABO*A/ABO*B", r


def test_abo_BB_hom_B():
    r = oa.abo_assign_type(_pile(1.0, 1.0), _DEFS)
    assert r["ABO_type"] == "B" and r["alleles"] == "ABO*B/ABO*B", r


def test_abo_BO_het_insertion_plus_het_B():
    r = oa.abo_assign_type(_pile(0.5, 0.5), _DEFS)
    assert r["ABO_type"] == "B" and r["alleles"] == "ABO*B/ABO*O", r


def test_abo_no_coverage_is_not_typed_never_default():
    assert oa.abo_assign_type(None, _DEFS)["ABO_type"] == "Not typed"
    assert oa.abo_assign_type({"sites": {}}, _DEFS)["ABO_type"] == "Not typed"
    r = oa.abo_assign_type(_pile(0.0, 0.0, depth=4), _DEFS)
    assert r["ABO_type"] == "Not typed" and r["confidence"] == "none", r


def test_abo_non_o_without_B_site_coverage_is_not_typed():
    r = oa.abo_assign_type(_pile(1.0, 0.0, cover_b=False), _DEFS)
    assert r["ABO_type"] == "Not typed", r
    # but O/O needs only the O site
    assert oa.abo_assign_type(_pile(0.0, 0.0, cover_b=False), _DEFS)["ABO_type"] == "O"


def test_abo_pileup_parser_counts_insertions():
    import abo_pileup as ap
    p = ap.parse_pileup_bases("..,+1c..+1C,^F.$,-2AC*", "T")
    assert p["ref"] == 8 and p["ins"] == ["C", "C"] and p["dele"] == ["AC"] and p["star"] == 1, p


def test_abo_tsv_round_trip():
    d = tempfile.mkdtemp()
    p = os.path.join(d, "abo", "abo_type.tsv")
    oa.write_abo_tsv(p, "S1", {"ABO_type": "A", "alleles": "ABO*A/ABO*O",
                               "confidence": "high", "flags": [], "n_sites_observed": 5,
                               "evidence": {"o_site_depth": 33}})
    row = oa.read_abo_tsv(p)
    assert row["ABO_type"] == "A" and row["confidence"] == "high" and row["o_site_depth"] == "33", row


def test_abo_defs_file_coordinates_and_reference_is_O():
    defs = oa.load_abo_defs()
    pos = {v["rsid"]: v["pos_grch38"] for v in defs["variants"]}
    assert pos["rs8176719"] == 133257521 and pos["rs8176746"] == 133255935, pos
    assert pos["rs8176747"] == 133255928 and pos["rs7853989"] == 133256205, pos
    assert "NON-O allele is the one-base" in defs["reference_allele_note"], defs
    assert len(defs["ab_discriminating_rsids"]) == 4, defs


if __name__ == "__main__":
    for _name, _fn in sorted(globals().items()):
        if _name.startswith("test_") and callable(_fn):
            _fn()
    print("all omnigen_addons self-checks passed")
