"""Self-checks for hla2_consensus.py (synthetic inputs only; no BAM, no container)."""
import os
import tempfile

import hla2_consensus as hc


def test_fields_and_rare():
    assert hc.fields("HLA-DRB1*10:01:01:02", 2) == "DRB1*10:01"
    assert hc.fields("DRB1*16:76", 1) == "DRB1*16"
    assert hc.is_rare_name("DRB1*11:284") and not hc.is_rare_name("DRB1*11:01:01")


def test_agree_2field_order_independent():
    r = hc.consensus_gene("DRB1", ["DRB1*10:01:01", "DRB1*04:02:01"],
                          {"alleles": [["DRB1*04:02:01"], ["DRB1*10:01:03"]],
                           "abundance": [30, 30], "quality": [20, 20], "secondary": []})
    assert r["agreement"] == "agree_2field", r
    assert sorted(r["consensus"]) == ["DRB1*04:02", "DRB1*10:01"], r


def test_agree_1field_collapses_rare_name():
    r = hc.consensus_gene("DRB1", ["DRB1*11:284", "DRB1*13:01:01"],
                          {"alleles": [["DRB1*11:01:01"], ["DRB1*13:01:01"]],
                           "abundance": [30, 30], "quality": [20, 20], "secondary": []})
    assert r["agreement"] == "agree_1field", r
    assert "DRB1*11" in r["consensus"] and "DRB1*13:01" in r["consensus"], r
    assert r["collapsed_rare"] == ["DRB1*11:284"], r


def test_disagree_never_reports_a_guess():
    r = hc.consensus_gene("DRB1", ["DRB1*10:01:01", "DRB1*13:01:01"],
                          {"alleles": [["DRB1*15:01:01"], ["DRB1*16:02:01"]],
                           "abundance": [30, 30], "quality": [20, 20], "secondary": []})
    assert r["agreement"] == "disagree" and r["consensus"] == ["DRB1*?", "DRB1*?"], r
    assert r["consensus_resolution"] == "none", r


def test_partial_disagree_keeps_shared_allele():
    r = hc.consensus_gene("DRB1", ["DRB1*10:01:01", "DRB1*13:01:01"],
                          {"alleles": [["DRB1*10:01:01"], ["DRB1*16:02:01"]],
                           "abundance": [30, 30], "quality": [20, 20], "secondary": []})
    assert r["agreement"] == "disagree" and r["consensus"] == ["DRB1*10:01", "DRB1*?"], r
    assert r["consensus_resolution"] == "partial", r


def test_tied_t1k_options_match_group():
    r = hc.consensus_gene("DRB1", ["DRB1*04:02:01", "DRB1*07:01:01"],
                          {"alleles": [["DRB1*07:01:01"], ["DRB1*04:01:01", "DRB1*04:04:01"]],
                           "abundance": [30, 18], "quality": [21, 11], "secondary": []})
    assert r["agreement"] == "agree_1field" and "DRB1*04" in r["consensus"], r


def test_single_typer_and_not_typed_are_explicit():
    r = hc.consensus_gene("DQB1", ["DQB1*02:189", "DQB1*05:01:01"], None)
    assert r["agreement"] == "single_typer" and r["consensus"][0] == "DQB1*02", r
    assert any("no independent confirmation" in f for f in r["flags"]), r
    r = hc.consensus_gene("DQB1", None, None)
    assert r["agreement"] == "not_typed" and r["consensus"] == [], r


def test_low_support_flags():
    ev = {"reads_aligned": 12, "pairs": [{"alleles": ["a", "b"], "explained_pct": 50.0}],
          "abundance": [], "zygosity": None}
    r = hc.consensus_gene("DRB1", ["DRB1*01:01:01", "DRB1*03:01:01"],
                          {"alleles": [["DRB1*01:01:01"], ["DRB1*03:01:01"]],
                           "abundance": [3.0, 30.0], "quality": [0, 20], "secondary": []}, ev)
    kinds = " ".join(r["flags"])
    assert "arcasHLA 12 reads" in kinds and "explains 50%" in kinds and "T1K DRB1*01:01:01" in kinds, r


def test_parsers_roundtrip():
    d = tempfile.mkdtemp()
    p = os.path.join(d, "t.tsv")
    open(p, "w").write("HLA-DRB1\t2\tHLA-DRB1*10:01:01\t30.1\t21\tHLA-DRB1*04:04:01\t11.0\t6\tHLA-DRB1*04:02:01;12.4;7\n"
                       "HLA-DQA1\t1\tHLA-DQA1*01:05:01\t48.1\t30\t.\t0\t-1\t\n")
    t = hc.parse_t1k_tsv(p)
    assert t["DRB1"]["alleles"] == [["DRB1*10:01:01"], ["DRB1*04:04:01"]], t
    assert t["DRB1"]["secondary"][0]["allele"] == "DRB1*04:02:01", t
    assert t["DQA1"]["alleles"][0] == t["DQA1"]["alleles"][1], t   # '.' = homozygous
    v = os.path.join(d, "v.log")
    open(v, "w").write(
        "[genotype] Genotyping HLA-DRB1\n[genotype] 159 reads aligned to HLA-DRB1 in 52 classes\n"
        "[genotype] Top alleles by abundance:\n\t\tallele   abundance\n\t\tDRB1*04:02:01              36.12%\n\n"
        "[genotype] Pairs by % explained reads:\n\t\tallele pair   explained\n"
        "\t\tDRB1*10:01:16, DRB1*04:02:01        84.28%\n\n[genotype] Checking zygosity\n"
        "[genotype] Likely heterozygous: minor/major nonshared count 0.79\n")
    ev = hc.parse_arcas_verbose(v)["DRB1"]
    assert ev["reads_aligned"] == 159 and ev["pairs"][0]["explained_pct"] == 84.28, ev
    assert ev["abundance"][0]["allele"] == "DRB1*04:02:01" and "heterozygous" in ev["zygosity"], ev


if __name__ == "__main__":
    for _n, _f in sorted(globals().items()):
        if _n.startswith("test_") and callable(_f):
            _f()
    print("all hla2_consensus self-checks passed")
