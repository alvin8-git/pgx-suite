#!/usr/bin/env bash
# Smoke test: run the Snakemake pipeline on a small GIAB HG002 slice (4 genes)
# and check the reconciled verdict and diplotype of each gene.
#
#   tests/smoke/run.sh [REF_DIR]
#
# REF_DIR  directory holding hg38.fa + hg38.fa.fai (default: <repo>/GRCh38)
# env:     IMAGE (default pgx-suite:latest), SMOKE_OUT (default <repo>/results/smoke),
#          CORES (default 4)
# Needs the same resources as run_pgx_suite.sh: the built image, StellarPGx/,
# StellarPGx/containers/ (stellarpgx-dev.sif, pharmcat.sif) and pypgx/pypgx-bundle/.
# Exit 0 = every expected value matched; 1 = mismatch; 2 = pipeline or setup failure.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
REF_DIR="$(realpath "${1:-$REPO/GRCh38}")"
IMAGE="${IMAGE:-pgx-suite:latest}"
OUT="${SMOKE_OUT:-$REPO/results/smoke}"
GENES="$(cut -f1 "$HERE/expected.tsv" | grep -v '^#' | tr '\n' ' ')"

[[ -f "$REF_DIR/hg38.fa.fai" ]] || { echo "smoke: no hg38.fa.fai in $REF_DIR" >&2; exit 2; }
for d in "$REPO/StellarPGx/containers" "$REPO/pypgx/pypgx-bundle"; do
    [[ -d "$d" ]] || { echo "smoke: missing $d" >&2; exit 2; }
done

rm -rf "$OUT"; mkdir -p "$OUT"
start=$(date +%s)
docker run --privileged --rm \
    -v "$REPO/StellarPGx:/pgx/stellarpgx" \
    -v "$REPO/StellarPGx/containers:/pgx/containers" \
    -v "$REPO/pypgx/pypgx-bundle:/pgx/bundle" \
    -v "$REF_DIR:/pgx/ref:ro" \
    -v "$HERE/data:/pgx/data:ro" \
    -v "$OUT:/pgx/results" \
    "$IMAGE" \
    snakemake -s /opt/pgx/Snakefile --cores "${CORES:-4}" --quiet rules \
        --config bam=/pgx/data/HG002.smoke.bam outdir=/pgx/results sample=HG002 \
                 genes="$GENES" \
    > "$OUT/smoke.log" 2>&1 || { echo "smoke: pipeline failed, see $OUT/smoke.log" >&2; exit 2; }

# Compare status, diplotype, nAgree, nCalled per gene against expected.tsv
fail=0
while IFS=$'\t' read -r gene want; do
    [[ "$gene" == \#* || -z "$gene" ]] && continue
    got="$(awk -F'\t' -v g="$gene" '$1==g {print $2"\t"$3"\t"$5"\t"$6}' "$OUT/all_genes_summary.tsv")"
    if [[ "$got" == "$want" ]]; then
        printf 'PASS  %-8s %s\n' "$gene" "${want//$'\t'/ }"
    else
        printf 'FAIL  %-8s expected [%s] got [%s]\n' "$gene" "${want//$'\t'/ }" "${got//$'\t'/ }"
        fail=1
    fi
done < "$HERE/expected.tsv"

[[ -s "$OUT/HG002_pgx_report.html" ]] || { echo "FAIL  HTML report not written"; fail=1; }
grep -q '"alt_aware": "yes"' "$OUT/alt_aware.json" && echo "PASS  altcheck alt_aware=yes" \
    || { echo "FAIL  altcheck did not report alt_aware=yes"; fail=1; }
echo "smoke: $(( $(date +%s) - start )) s, output in $OUT"
exit $fail
