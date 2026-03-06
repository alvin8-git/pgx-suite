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
echo "[1/4] PyPGx"
run_test "pypgx --version"             pypgx --version
run_test "pypgx import"                python3 -c "import pypgx; print('OK')"
run_test "pypgx run-ngs-pipeline help" pypgx run-ngs-pipeline --help
# Note: full pipeline test requires a BAM file (provided in Phase 2)
echo ""

# ── Stargazer ─────────────────────────────────────────────────────────────
echo "[2/4] Stargazer"
run_test "stargazer --version"   stargazer --version
run_test "stargazer wrapper exists" test -x /usr/local/bin/stargazer

# Functional test: VCF-only mode (no BAM/GDF needed; bundled 70-sample VCF)
# Note: -c and -g (CNV analysis) require grc38-aligned depth data; VCF-only mode
# works with the hg19-aligned example VCF since variant coordinates are remapped.
run_test "stargazer CYP2D6 VCF-only (grc38)" \
    bash -c "cd /opt/stargazer/example && stargazer \
        -o /tmp/stargazer_test_out \
        -d wgs \
        -t cyp2d6 \
        -a grc38 \
        -i getrm-cyp2d6-vdr.joint.filtered.vcf"
echo ""

# ── Aldy ──────────────────────────────────────────────────────────────────
echo "[3/4] Aldy"
run_test "aldy version"    python3 -c "import aldy; print('aldy', aldy.__version__)"
run_test "aldy import"     python3 -c "import aldy; print('OK')"
run_test "aldy built-in test suite (78 tests)" \
    aldy test
echo ""

# ── StellarPGx ────────────────────────────────────────────────────────────
echo "[4/4] StellarPGx"
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

# ── Summary ───────────────────────────────────────────────────────────────
echo "============================================================"
echo " Results: ${PASS} PASSED  |  ${FAIL} FAILED"
echo "============================================================"
[[ "${FAIL}" -eq 0 ]]
