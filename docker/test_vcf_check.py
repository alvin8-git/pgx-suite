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

    # single-SNP genes (path1): VKORC1 / IFNL3 defs + parser registration
    for g, rs, pos in (("VKORC1", "rs9923231", 31096368), ("IFNL3", "rs12979860", 39248147)):
        snp = data[g]["snp"]
        assert data[g]["mode"] == "single_snp"
        assert snp["rsid"] == rs and snp["pos"] == pos
        assert snp["report_ref"] and snp["report_alt"] and "phenotypes" in snp
        assert comp._VCF_CHECK_PARSERS.get(g) is comp.parse_single_snp_vcf
        r = comp.parse_single_snp_vcf("/nonexistent", g, "x")
        assert r.status == "failed" and "VCF not found" in r.diplotype

    # VCF-authoritative verdict policy: the VCF-Check call wins over star-allele
    # callers using incompatible nomenclature.
    assert all(comp.AUTHORITATIVE[g] == "VCF-Check"
               for g in ("RYR1", "CACNA1S", "G6PD", "VKORC1", "IFNL3"))
    CR = comp.CallerResult
    def _r(tool, dip, st="ok", phe="-"):
        c = CR(tool=tool); c.status, c.diplotype, c.phenotype = st, dip, phe; return c
    rs = [_r("PyPGx", "Reference/Reference"), _r("Stargazer", "*1/*1"),
          _r("Aldy", "*H1/*H1"), _r("VCF-Check", "G/A", phe="Intermediate warfarin sensitivity (G/A)")]
    v = comp.compute_verdict(rs, no_call=False, gene="VKORC1")
    assert v["status"] == "concordant" and v["consensus_diplotype"] == "G/A"
    assert v.get("authority") == "VCF-Check"

    # ── G6PD: full CPIC table, sex-aware, never "Normal" on a partial list ────
    g6 = data["G6PD"]
    assert len(g6["sites"]) == 171 and len(g6["alleles"]) >= 180
    gaohe = (154546061, "T", ["C"], [1, 1])           # rs137852340, c.95A>G
    full = {s["pos"]: 30 for s in g6["sites"]}
    c = comp.g6pd_call([gaohe], g6, full, "MALE")
    assert "G6PD Deficient" in c["phenotype"] and "hemizygous" in c["phenotype"], c
    assert c["diplotype"].startswith("Gaohe(hemizygous)"), c
    c = comp.g6pd_call([(154546061, "T", ["C"], [0, 1])], g6, full, "FEMALE")
    assert c["phenotype"].startswith("G6PD Variable"), c
    c = comp.g6pd_call([gaohe], g6, full, "FEMALE")
    assert c["phenotype"].startswith("G6PD Deficient") and "two deficient" in c["phenotype"], c
    # reference: wording is "not found among checked sites", never "Normal"
    partial = dict(full); partial[154546061] = 0
    c = comp.g6pd_call([], g6, partial, "MALE")
    assert "among the sites checked" in c["phenotype"] and "170 of 171" in c["phenotype"], c
    assert "normal" not in c["phenotype"].lower(), c
    # an untyped site must not be read as reference: Gaohe site has no depth, no call
    assert comp.g6pd_call([], g6, {}, "MALE")["diplotype"] == "not typed"
    # G202A alone (no 376G) is CPIC allele Asahi (class III) -> Deficient
    c = comp.g6pd_call([(154536002, "C", ["T"], [1, 1])], g6, full, "MALE")
    assert [x["name"] for x in c["alleles"]] == ["Asahi"] and "Deficient" in c["phenotype"], c
    # A- (202A+376G) resolves to the specific allele, not the nested A allele
    c = comp.g6pd_call([(154536002, "C", ["T"], [1]), (154535277, "T", ["C"], [1])], g6, full, "MALE")
    assert [x["name"] for x in c["alleles"]] == ["A- 202A_376G"], c
    # nested allele on leftover copies: 1/1 202A + 0/1 376G (female) = A- + Asahi = two deficient
    c = comp.g6pd_call([(154536002, "C", ["T"], [1, 1]), (154535277, "T", ["C"], [0, 1])], g6, full, "FEMALE")
    assert sorted(x["name"] for x in c["alleles"]) == ["A- 202A_376G", "Asahi"], c
    assert c["phenotype"].startswith("G6PD Deficient") and "two deficient" in c["phenotype"], c
    # a 0/1 on chrX in a male is flagged, never promoted to a hemizygous Deficient call
    c = comp.g6pd_call([(154546061, "T", ["C"], [0, 1])], g6, full, "MALE")
    assert c["phenotype"].startswith("Indeterminate") and c["male_het"] == ["154546061C"], c
    # REF padded with trailing context (bcftools) still matches the CPIC deletion (Brighton, class I)
    c = comp.g6pd_call([(154532055, "CTCTT", ["CT"], [1], 50.0, 20)], g6, full, "MALE")
    assert [x["name"] for x in c["alleles"]] == ["Brighton"] and "Deficient" in c["phenotype"], c
    # a record at a CPIC site that cannot be matched is never ignored: no negative over it
    c = comp.g6pd_call([(154532055, "CTCT", ["CA"], [1], 50.0, 20)], g6, full, "MALE")
    assert c["phenotype"].startswith("Indeterminate") and "No G6PD deficiency" not in c["phenotype"], c
    # a low-quality record is never a Deficient call, and never a clean negative
    c = comp.g6pd_call([(154546061, "T", ["C"], [1, 1], 4.0, 1)], g6, full, "MALE")
    assert c["phenotype"].startswith("Indeterminate") and "low-quality" in c["phenotype"], c
    # hemizygous deficiency + a male 0/1 artifact elsewhere: Deficient stays the leading phrase
    c = comp.g6pd_call([gaohe, (154536002, "C", ["T"], [0, 1])], g6, full, "MALE")
    assert c["phenotype"].startswith("G6PD Deficient") and "QC:" in c["phenotype"], c
    # another caller's "Variable" is dissent; the discordant phenotype never reads as a negative
    neg = _r("VCF-Check", "Reference (x)", phe="No G6PD deficiency variant found among the sites checked")
    v = comp.compute_verdict([_r("PharmCAT", "B/Viangchan", phe="Variable"), neg], False, "G6PD")
    assert v["status"] == "discordant" and v["consensus_phenotype"].startswith("Indeterminate"), v
    # VCF-Check "G6PD not typed ..." must not be turned into a false discordance
    nt = _r("VCF-Check", "not typed", phe="G6PD not typed (key site(s) below 10x: Mediterranean 3x)")
    v = comp.compute_verdict([_r("PyPGx", "MALE/Gaohe", phe="G6PD Deficient"), nt], False, "G6PD")
    assert v["status"] != "discordant", v
    # authority must not hide a deficiency reported by other callers
    ref = _r("VCF-Check", "Reference (x)", phe="No G6PD deficiency variant found among the sites checked")
    rs = [_r("PyPGx", "MALE/Gaohe", phe="G6PD Deficient"), ref]
    assert comp.compute_verdict(rs, False, "G6PD")["status"] == "discordant"

    # ── MT-RNR1 table: CPIC increased-risk alleles at the right rCRS positions ──
    mt = {v["label"]: v for v in data["MT-RNR1"]["variants"]}
    for lab, pos in (("m.1555A>G", 1555), ("m.1494C>T", 1494), ("m.1095T>C", 1095)):
        assert mt[lab]["pos"] == pos and mt[lab]["cpic_class"] == "increased risk", lab

    # ── ERROR is never NEGATIVE: depth unknown / key site low -> not typed ──
    g6 = data["G6PD"]
    full = {s_["pos"]: 50 for s_ in g6["sites"]}
    assert comp.g6pd_call([], g6, full, "MALE")["diplotype"].startswith("Reference")
    c = comp.g6pd_call([], g6, None, "MALE")          # depth could not be checked
    assert c["diplotype"] == "not typed" and "could not be checked" in c["phenotype"], c
    assert "No G6PD deficiency" not in c["phenotype"]
    for pos, lab in ((154534419, "c.563C>T"), (154533122, "c.871G>A"), (154534495, "c.487G>A"),
                     (154532269, "c.1376G>T"), (154532257, "c.1388G>A"), (154532390, "c.1360C>T"),
                     (154532969, "c.1024C>T"), (154536002, "A-"), (154535277, "A-")):
        low = dict(full); low[pos] = 3
        c = comp.g6pd_call([], g6, low, "MALE")
        assert c["diplotype"] == "not typed" and lab in c["phenotype"] and "3x" in c["phenotype"], (pos, c)
    # sex-aware key-site depth (>= 5 reads per X copy): male 6x typed, male 4x not, female 8x not
    low = dict(full); low[154532257] = 6
    assert comp.g6pd_call([], g6, low, "MALE")["diplotype"].startswith("Reference")
    low[154532257] = 4
    assert comp.g6pd_call([], g6, low, "MALE")["diplotype"] == "not typed"
    low[154532257] = 8
    assert comp.g6pd_call([], g6, low, "FEMALE")["diplotype"] == "not typed"
    assert comp.g6pd_call([], g6, low, None)["diplotype"] == "not typed"
    # a non-key site below min_cov does not block a negative
    low = dict(full); low[154546061] = 0              # Gaohe is not a key site
    assert comp.g6pd_call([], g6, low, "MALE")["diplotype"].startswith("Reference")
    # a positive stays positive with depth unknown or key sites low
    gaohe_rec = [(154546061, "T", ["C"], [1])]
    assert "Deficient" in comp.g6pd_call(gaohe_rec, g6, None, "MALE")["phenotype"]
    low = dict(full); low[154534419] = 0
    assert "Deficient" in comp.g6pd_call(gaohe_rec, g6, low, "MALE")["phenotype"]

    # ── pgx-mt.sh: an unparseable mutserve data row aborts, never skipped ──
    import re, subprocess, tempfile
    sh = open(os.path.join(HERE, "pgx-mt.sh")).read()
    code = re.search(r"python3 - <<'PYEOF'[^\n]*\n(.*?)\nPYEOF\n", sh, re.S).group(1)
    with tempfile.TemporaryDirectory() as td:
        pos = next(iter(data["MT-RNR1"]["variants"]))["pos"]
        dep = os.path.join(td, "d.tsv"); open(dep, "w").write(f"chrM\t{pos}\t500\n")
        def run(row):
            raw = os.path.join(td, "raw.txt")
            open(raw, "w").write("ID\tFilter\tPos\tRef\tVariant\tVariantLevel\tCoverage\n" + row + "\n")
            env = dict(os.environ, SAMPLE="S", MUTSERVE_RAW=raw, DEPTH=dep,
                       VARTABLE=os.path.join(HERE, "vcf_check_variants.json"))
            return subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True)
        assert run(f"S\tPASS\t{pos}\tA\tG\t0.99\t500").returncode == 0
        bad = run(f"S\tPASS\tnotanumber\tA\tG\t0.99\t500")
        assert bad.returncode != 0 and "could not be parsed" in bad.stderr, bad
        assert run("S\tPASS\t1555").returncode != 0      # short row (IndexError)
        # a filtered call at an increased-risk site is never read as reference (all sites at 500x)
        open(dep, "w").write("".join(f"chrM\t{v['pos']}\t500\n" for v in data["MT-RNR1"]["variants"]))
        assert '"diplotype": "Reference"' in run("S\tPASS\t3000\tA\tG\t0.99\t500").stdout
        f = run("S\tSTRAND_BIAS\t1555\tA\tG\t0.99\t2000")
        assert f.returncode == 0 and '"diplotype": "not typed"' in f.stdout and "STRAND_BIAS" in f.stdout, f
        # m.807 has two table alleles (A>G, A>C): both must be reportable
        open(dep, "w").write("".join(f"chrM\t{v['pos']}\t500\n" for v in data["MT-RNR1"]["variants"]))
        ok = run(f"S\tPASS\t807\tA\tG\t0.99\t500")
        assert ok.returncode == 0 and "m.807A>G" in ok.stdout, ok

    print(f"all vcf_check checks passed — RYR1 {len(variants)} MH variants + VKORC1/IFNL3 single-SNP")


def test_vcf_check():
    main()


if __name__ == "__main__":
    try:
        main()
    except AssertionError as e:
        print("FAIL:", e, file=sys.stderr)
        sys.exit(1)
