#!/usr/bin/env python3
"""Tests for reconcile.norm_diplotype and the nAgree/nCalled count in compute_verdict.

Run: python3 docker/test_normalise.py
"""
import importlib.util, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reconcile as rc

spec = importlib.util.spec_from_file_location("pgx_compare", os.path.join(HERE, "pgx-compare.py"))
cmp = importlib.util.module_from_spec(spec); spec.loader.exec_module(cmp)


def eq(a, b, gene="CYP2D6"):
    ka, kb = rc.norm_diplotype(gene, a), rc.norm_diplotype(gene, b)
    assert ka is not None and ka == kb, (a, b, ka, kb)


def ne(a, b, gene="CYP2D6"):
    assert rc.norm_diplotype(gene, a) != rc.norm_diplotype(gene, b), (a, b)


def r(tool, dip, pheno="-"):
    c = cmp.CallerResult(tool=tool)
    c.status, c.diplotype, c.phenotype = "ok", dip, pheno
    return c


def test_rules():
    eq("*2/*36X+*10", "*2/*36+*10")            # Aldy tandem marker
    eq("*3/*4.021", "*3/*4")                   # sub-allele suffix
    eq("*6B/*35", "*6/*35")                    # CYP2D6 legacy letter
    eq("*49/*36+*10", "*36+*10/*49")           # tandem + haplotype order
    eq("*2 / *4", "*4/*2")
    ne("*2/*41", "*2+rs61737946+rs61745683/*41")   # rs tags are real differences
    ne("*2/*2", "*2/*2x2")                     # copy number is real
    ne("*3/*4", "*3/*68+*4")
    eq("*3A/*1", "*1/*3A", gene="TPMT")
    ne("*3A/*1", "*3C/*1", gene="TPMT")        # letters kept outside CYP2D6
    # VCF-Check: single listed allele only
    eq("*6(het) [unphased]", "*1/*6", gene="UGT1A1")
    eq("*6(hom) [unphased]", "*6/*6", gene="UGT1A1")
    assert rc.norm_diplotype("UGT1A1", "*6(het)+*60(het) [unphased]") is rc.UNRESOLVED
    assert rc.norm_diplotype("CYP2D6", "-") == "-"


def test_nagree_authority_gene():
    # CYP2D6: Cyrius is the authority; all five callers say *2/*4 -> 5/5, not 1/5
    v = cmp.compute_verdict([r(t, "*2/*4") for t in
                             ("PyPGx", "Stargazer", "Aldy", "StellarPGx", "Cyrius")],
                            no_call=False, gene="CYP2D6")
    assert (v["n_agree"], v["n_called"]) == (5, 5), v
    # spelling-only differences still agree
    v = cmp.compute_verdict([r("PyPGx", "*2/*36+*10"), r("Aldy", "*2/*36X+*10"),
                             r("Cyrius", "*36+*10/*2")], no_call=False, gene="CYP2D6")
    assert (v["n_agree"], v["n_called"]) == (3, 3), v
    # real dispute: Cyrius authority, two callers differ
    v = cmp.compute_verdict([r("PyPGx", "*2/*2"), r("Stargazer", "*2/*2x2"),
                             r("Cyrius", "*2/*2x2")], no_call=False, gene="CYP2D6")
    assert (v["n_agree"], v["n_called"]) == (2, 3), v


def test_nagree_ugt1a1_unphased():
    # PharmCAT authority; four star callers on *1/*6; VCF-Check lists two het alleles
    rows = [r(t, "*1/*6") for t in ("PyPGx", "Stargazer", "Aldy", "PharmCAT")]
    rows.append(r("VCF-Check", "*6(het)+*60(het) [unphased]"))
    v = cmp.compute_verdict(rows, no_call=False, gene="UGT1A1")
    assert (v["n_agree"], v["n_called"]) == (4, 4), v
    assert "VCF-Check" in v.get("note", ""), v
    # single het allele maps to *1/*6 and counts as agreeing
    rows[-1] = r("VCF-Check", "*6(het) [unphased]")
    v = cmp.compute_verdict(rows, no_call=False, gene="UGT1A1")
    assert (v["n_agree"], v["n_called"]) == (5, 5), v


def test_nagree_vote_gene():
    v = cmp.compute_verdict([r("PyPGx", "*1/*4"), r("Aldy", "*4/*1"), r("Stargazer", "*1/*1")],
                            no_call=False, gene="CYP2C9")
    assert (v["n_agree"], v["n_called"], v["status"]) == (2, 3, "majority"), v


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
    print("normalise + nAgree tests: PASS")
