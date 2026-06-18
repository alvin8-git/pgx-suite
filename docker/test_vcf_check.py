#!/usr/bin/env python3
"""Schema + registration guard for the RYR1 VCF-Check variant table (design D).

The RYR1 malignant-hyperthermia variant list is curated as data
(vcf_check_variants.json), not code. This test fails if the table drifts out of
shape, loses coordinates, or the parser stops being wired into pgx-compare.
Pure-Python: no bcftools, so the parser's VCF query path is checked only for
graceful failure (missing VCF), not live calling.
"""
import importlib.util
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RYR1_REGION = (38430690, 38590564)  # chr19, must match genes.tsv


def _load_compare():
    spec = importlib.util.spec_from_file_location("pgxcompare", os.path.join(HERE, "pgx-compare.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    with open(os.path.join(HERE, "vcf_check_variants.json"), encoding="utf-8") as fh:
        data = json.load(fh)
    ryr1 = data["RYR1"]
    variants = ryr1["variants"]
    assert ryr1["chrom"] == "chr19"
    assert len(variants) >= 90, f"expected >=90 RYR1 MH variants, got {len(variants)}"

    seen = {}
    for v in variants:
        for k in ("pos", "ref", "alt", "hgvs_c", "status"):
            assert v.get(k) not in (None, ""), f"{v.get('hgvs_c')} missing {k}"
        assert v["status"] in ("Pathogenic", "Likely Pathogenic"), \
            f"{v['hgvs_c']} unexpected status {v['status']}"
        assert RYR1_REGION[0] <= v["pos"] <= RYR1_REGION[1], \
            f"{v['hgvs_c']} pos {v['pos']} outside RYR1 region"
        assert set(v["ref"]) <= set("ACGT") and set(v["alt"]) <= set("ACGT"), \
            f"{v['hgvs_c']} non-ACGT alleles"
        seen[v["hgvs_c"]] = v

    # regression anchor: the validation variant the star-allele callers all missed
    assert "c.1021G>A" in seen, "c.1021G>A (P_2000 validation variant) missing from table"
    assert seen["c.1021G>A"]["pos"] == 38448712 and seen["c.1021G>A"]["alt"] == "A"

    # parser is registered and degrades gracefully without a VCF
    comp = _load_compare()
    assert "RYR1" in comp._VCF_CHECK_PARSERS, "parse_ryr1_vcf not registered"
    res = comp.parse_ryr1_vcf("/nonexistent", "RYR1", "x")
    assert res.status == "failed" and "VCF not found" in res.diplotype

    print(f"all vcf_check (RYR1) checks passed — {len(variants)} MH variants")


if __name__ == "__main__":
    try:
        main()
    except AssertionError as e:
        print("FAIL:", e, file=sys.stderr)
        sys.exit(1)
