#!/usr/bin/env bash
# /opt/pgx/test.sh — Phase 1 smoke tests (no user BAM/reference needed)
# Run: docker run --privileged ... pgx-suite:latest bash /opt/pgx/test.sh
set -euo pipefail

PASS=0
FAIL=0

run_test() {
    local name="$1"
    shift
    printf "  %-50s" "$name"
    if "$@" > /tmp/_pgx_test_out 2>&1; then
        echo "PASS"
        ((PASS++)) || true
    else
        echo "FAIL"
        head -5 /tmp/_pgx_test_out | sed 's/^/    /'
        ((FAIL++)) || true
    fi
}

echo "============================================================"
echo " PGx Suite Docker Container — Smoke Tests"
echo " Reference build: GRCh38"
echo "============================================================"
echo ""

# ── PyPGx ─────────────────────────────────────────────────────────────────
echo "[1/8] PyPGx"
run_test "pypgx --version"             pypgx --version
run_test "pypgx import"                python3 -c "import pypgx; print('OK')"
run_test "pypgx run-ngs-pipeline help" pypgx run-ngs-pipeline --help
# Note: full pipeline test requires a BAM file (provided in Phase 2)
echo ""

# ── Stargazer ─────────────────────────────────────────────────────────────
echo "[2/8] Stargazer"
run_test "stargazer --version"   stargazer --version
run_test "stargazer wrapper exists" test -x /usr/local/bin/stargazer

# Functional test: VCF-only mode. The example/ dir is excluded from the lean
# GRCh38-only image (.dockerignore: stargazer-grc38-2.0.3/example/), so this
# SKIPs unless that data is mounted in. The Stargazer binary itself is covered
# by --version above and by the real per-sample runs in pgx-run.sh.
_STARGAZER_EX="/opt/stargazer/example/getrm-cyp2d6-vdr.joint.filtered.vcf"
if [[ -f "$_STARGAZER_EX" ]]; then
    run_test "stargazer CYP2D6 VCF-only (grc38)" \
        bash -c "cd /opt/stargazer/example && stargazer \
            -o /tmp/stargazer_test_out \
            -d wgs \
            -t cyp2d6 \
            -a grc38 \
            -i getrm-cyp2d6-vdr.joint.filtered.vcf"
else
    echo "  [SKIP] stargazer VCF-only — example data not in lean image"
    echo "         (stargazer-grc38-2.0.3/example/ excluded by .dockerignore)"
fi
echo ""

# ── Aldy ──────────────────────────────────────────────────────────────────
echo "[3/8] Aldy"
run_test "aldy version"    python3 -c "import aldy; print('aldy', aldy.__version__)"
run_test "aldy import"     python3 -c "import aldy; print('OK')"
run_test "aldy built-in test suite (78 tests)" \
    aldy test
echo ""

# ── StellarPGx ────────────────────────────────────────────────────────────
echo "[4/8] StellarPGx"
run_test "nextflow --version"    nextflow -version
run_test "apptainer --version"   apptainer --version

if [[ -d "/pgx/stellarpgx" && -f "/pgx/containers/stellarpgx-dev.sif" ]]; then
    run_test "StellarPGx test pipeline (CYP2D6/hg38)" \
        bash -c "cd /pgx/stellarpgx && nextflow run main.nf \
            -profile standard,test \
            --out_dir /tmp/stellarpgx_test_out"
else
    echo "  [SKIP] StellarPGx volumes not mounted"
    echo "         Requires: -v \$(pwd)/StellarPGx:/pgx/stellarpgx"
    echo "                   -v \$(pwd)/StellarPGx/containers:/pgx/containers"
fi
echo ""

# ── Cyrius (CYP2D6 authority) ───────────────────────────────────────────────
echo "[5/8] Cyrius"
run_test "cyrius star_caller present" test -f /opt/cyrius/star_caller.py
run_test "cyrius star_caller --help"  python3 /opt/cyrius/star_caller.py --help
echo ""

# ── PharmCAT (UGT1A1/CYP2B6/CYP4F2 authority) ───────────────────────────────
echo "[6/8] PharmCAT"
run_test "pharmcat positions VCF baked" test -f /opt/pgx/pharmcat_positions.vcf.bgz
if [[ -f "/pgx/containers/pharmcat.sif" ]]; then
    run_test "pharmcat_pipeline --version (SIF)" \
        apptainer exec /pgx/containers/pharmcat.sif pharmcat_pipeline --version
else
    echo "  [SKIP] pharmcat.sif not mounted — pull it (see README 'Pulling the PharmCAT SIF')"
fi
echo ""

# ── HLA / mtDNA callers ─────────────────────────────────────────────────────
echo "[7/8] OptiType + mutserve"
run_test "mutserve JAR valid"   java -jar /usr/local/bin/mutserve.jar --version
if [[ -f "/pgx/containers/optitype.sif" ]]; then
    run_test "optitype --help (SIF)" \
        apptainer exec /pgx/containers/optitype.sif OptiTypePipeline.py --help
else
    echo "  [SKIP] optitype.sif not mounted — pull it (see README 'Installation')"
fi
echo ""

# ── Reconciliation + verdict authority (VCF-Check, synonyms, phenotype tier) ─
echo "[8/8] Reconciliation & verdict authority"
run_test "reconcile.py self-test"        python3 /opt/pgx/reconcile.py --selftest
run_test "allele_synonyms.json loads"    python3 -c "import json; json.load(open('/opt/pgx/allele_synonyms.json'))"
run_test "vcf_check_variants.json loads"  python3 -c "import json; json.load(open('/opt/pgx/vcf_check_variants.json'))"
run_test "pgx-compare.py imports"         python3 -c "import importlib.util as u; s=u.spec_from_file_location('c','/opt/pgx/pgx-compare.py'); m=u.module_from_spec(s); s.loader.exec_module(m); assert m.AUTHORITATIVE['CYP2D6']=='Cyrius'"
run_test "pgx_altcheck.py self-test"      python3 /opt/pgx/pgx_altcheck.py --selftest
# ── OmniGen add-ons (HLA class II, ABO) ──────────────────────────────────────
run_test "omnigen_addons imports"         python3 -c "import importlib.util as u; s=u.spec_from_file_location('o','/opt/pgx/omnigen_addons.py'); m=u.module_from_spec(s); s.loader.exec_module(m); assert m.abo_assign_type({'rs8176719':None,'rs8176746':None,'rs8176747':None}, m.load_abo_defs())['ABO_type']=='A'"
run_test "abo_allele_defs.json UNVALIDATED" python3 -c "import json; d=json.load(open('/opt/pgx/abo_allele_defs.json')); assert 'UNVALIDATED' in d['validation_status']"
run_test "ABO+arcasHLA wired in compare"  python3 -c "import importlib.util as u; s=u.spec_from_file_location('c','/opt/pgx/pgx-compare.py'); m=u.module_from_spec(s); s.loader.exec_module(m); assert m.AUTHORITATIVE['ABO']=='VCF-Check'; assert m.GENE_SUPPORT['HLA-DQB1'].get('arcashla')"
run_test "hla2 script executable"         sh -c "test -x /opt/pgx/pgx-hla2.sh"
echo ""

# ── Summary ───────────────────────────────────────────────────────────────
echo "============================================================"
echo " Results: ${PASS} PASSED  |  ${FAIL} FAILED"
echo "============================================================"
[[ "${FAIL}" -eq 0 ]]
