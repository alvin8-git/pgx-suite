#!/usr/bin/env bash
# reanalyze_mgi.sh — re-run the full pgx-suite pipeline over the TTSH/MGI cohort
# from the source WGS BAMs, validating the latest pipeline (variant-tier
# reconciliation + Cyrius CYP2D6 + PharmCAT). Sequential (network-mounted BAMs,
# ~120-165 GB each); resumable via a per-sample sentinel.
#
# Usage: ./reanalyze_mgi.sh [JOBS]   (default JOBS=8 genes in parallel per sample)
set -uo pipefail

BAM_DIR="/home/alvin/storeData/TTSH/MGI_BAM"
OUT_PARENT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/results/MGI"
JOBS="${1:-8}"
LOG_DIR="/data/alvin/tmp/mgi_reanalysis"
mkdir -p "$LOG_DIR"

# All 20 cohort samples (basename before first dot)
SAMPLES=(C-5839 EQ_2017 EQ_2021 EQ-H05 EQ-H07 EQ-H08 H_342 H_361 \
         P_1143 P_2000 P_2026 P_2483 P_2505 P_2820 P_2912 P_62 P_670 P_93 P_944 P_982)

echo "=== MGI cohort reanalysis (latest pipeline) — $(date) ==="
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
echo "=== MGI cohort reanalysis complete — $(date) ==="
