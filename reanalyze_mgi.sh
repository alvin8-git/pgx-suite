#!/usr/bin/env bash
# reanalyze_mgi.sh — re-run the full pgx-suite pipeline over a validation cohort
# from the source WGS BAMs, validating the latest pipeline (variant-tier
# reconciliation + Cyrius CYP2D6 + PharmCAT). Sequential (network-mounted BAMs,
# ~120-165 GB each); resumable via a per-sample sentinel.
#
# Usage: BAM_DIR=/path/to/bams ./reanalyze_mgi.sh [JOBS]   (default JOBS=8 genes in parallel per sample)
# Cohort = every <sample>.bwa.sortdup.bqsr.bam in BAM_DIR; LOG_DIR defaults to ${TMPDIR:?set TMPDIR or LOG_DIR}/pgx_reanalysis.
set -uo pipefail

BAM_DIR="${BAM_DIR:?set BAM_DIR to the directory holding the cohort BAMs}"
OUT_PARENT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/results/MGI"
JOBS="${1:-8}"
LOG_DIR="${LOG_DIR:-${TMPDIR:?set TMPDIR or LOG_DIR}/pgx_reanalysis}"
mkdir -p "$LOG_DIR"

# Cohort samples (basename before first dot)
SAMPLES=()
for f in "$BAM_DIR"/*.bwa.sortdup.bqsr.bam; do [[ -e "$f" ]] && SAMPLES+=("$(basename "$f" .bwa.sortdup.bqsr.bam)"); done

echo "=== cohort reanalysis (latest pipeline) — $(date) ==="
echo "BAM dir: $BAM_DIR   Output: $OUT_PARENT   Jobs/sample: $JOBS"

for s in "${SAMPLES[@]}"; do
    bam="$BAM_DIR/${s}.bwa.sortdup.bqsr.bam"
    sentinel="$OUT_PARENT/${s}/.reanalyzed_v2"
    if [[ -f "$sentinel" ]]; then echo ">>> SKIP $s (sentinel present)"; continue; fi
    if [[ ! -f "$bam" ]]; then echo ">>> MISSING BAM $s — skipping"; continue; fi
    echo ">>> $(date '+%H:%M:%S')  START $s"
    if ./run_pgx_suite.sh "$bam" --output "$OUT_PARENT" --jobs "$JOBS" \
         > "$LOG_DIR/${s}.log" 2>&1; then
        touch "$sentinel"
        echo ">>> $(date '+%H:%M:%S')  DONE  $s"
    else
        echo ">>> $(date '+%H:%M:%S')  FAIL  $s (see $LOG_DIR/${s}.log)"
    fi
done
echo "=== cohort reanalysis complete — $(date) ==="
