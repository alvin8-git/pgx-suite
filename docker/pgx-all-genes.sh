#!/usr/bin/env bash
# docker/pgx-all-genes.sh — Run all supported PGx genes for a single sample
#
# Runs pgx-run.sh for every gene in the support matrix, in parallel batches.
# Each gene's stdout/stderr is written to its own log file.
# A summary table is printed at the end aggregating concordance across genes.
#
# Usage:
#   pgx-all-genes.sh <BAM_FILE> [--ref /pgx/ref/hg38.fa] \
#                    [--output /pgx/results] [--jobs 4]
#
# Notes:
#   - PyPGx prepare-depth-of-coverage is run once per SV gene, inside each
#     gene's own output directory.  For a large BAM this means the SV genes
#     each do one depth scan; running at lower --jobs limits concurrency.
#   - StellarPGx Nextflow writes a work/ directory.  Each gene-level call
#     from pgx-run.sh runs Nextflow in its own sub-directory, so they do not
#     collide.  If Nextflow hangs on the first run (JAR download), ensure
#     outbound internet access.
#   - POR and CYPOR are the same locus.  Only POR is run.
#   - GSTT1 is on an alt contig; only StellarPGx is called for it.

set -uo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
REF="/pgx/ref/hg38.fa"
OUTPUT="/pgx/results"
JOBS=4   # genes to run in parallel; raise to 8 on high-memory / high-I/O hosts

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF
Usage: pgx-all-genes.sh <BAM_FILE> [options]

Options:
  --ref PATH       GRCh38 reference FASTA (default: /pgx/ref/hg38.fa)
  --output PATH    Root output directory (default: /pgx/results)
  --jobs N         Genes to run in parallel (default: 4)
  -h, --help       Show this help

Output layout:
  <output>/<GENE>/          per-gene results from each tool
  <output>/logs/<GENE>.log  stdout + stderr for each gene
  <output>/all_genes_summary.tsv   concordance summary across all genes
EOF
}

# ── Argument parsing ───────────────────────────────────────────────────────────
if [[ $# -lt 1 ]]; then usage; exit 1; fi

BAM="$1"; shift

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)    REF="$2";    shift 2 ;;
        --output) OUTPUT="$2"; shift 2 ;;
        --jobs)   JOBS="$2";   shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

# ── Validate BAM ──────────────────────────────────────────────────────────────
if [[ ! -f "$BAM" ]]; then
    echo "ERROR: BAM file not found: $BAM" >&2; exit 1
fi
if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
    echo "ERROR: BAM index not found (expected ${BAM}.bai)" >&2; exit 1
fi
if [[ ! -f "$REF" ]]; then
    echo "ERROR: Reference FASTA not found: $REF" >&2; exit 1
fi

# ── Gene list — all unique supported genes ────────────────────────────────────
# POR/CYPOR are the same locus; only POR is kept here.
GENES=(
    CYP1A1 CYP1A2
    CYP2A6 CYP2B6 CYP2C8 CYP2C9 CYP2C19 CYP2D6 CYP2E1
    CYP3A4 CYP3A5 CYP4F2
    DPYD
    G6PD GSTM1 GSTT1
    IFNL3
    NAT1 NAT2
    NUDT15
    POR
    RYR1
    SLCO1B1
    TPMT
    UGT1A1
    VKORC1
)
TOTAL=${#GENES[@]}

# ── Directories ───────────────────────────────────────────────────────────────
LOG_DIR="${OUTPUT}/logs"
mkdir -p "$LOG_DIR"

# ── Extract sample name from BAM header ──────────────────────────────────────
SAMPLE=$(samtools view -H "$BAM" 2>/dev/null \
         | grep '^@RG' | grep -oP 'SM:\K[^\t]+' | head -1)
SAMPLE="${SAMPLE:-$(basename "$BAM" .bam)}"

# ── Banner ────────────────────────────────────────────────────────────────────
START_TS=$(date '+%Y-%m-%d %H:%M:%S')
cat <<EOF
========================================================================
 PGx Suite — All-genes batch run
 Sample:  ${SAMPLE}
 BAM:     ${BAM}
 Ref:     ${REF}
 Output:  ${OUTPUT}
 Genes:   ${TOTAL}
 Jobs:    ${JOBS} genes in parallel
 Started: ${START_TS}
========================================================================
EOF

# ── Per-gene runner (called in background) ────────────────────────────────────
run_gene() {
    local gene="$1"
    local log="${LOG_DIR}/${gene}.log"

    # Write timestamps to log only — do NOT tee to stdout because this function
    # is called as result=$(run_gene "$GENE") and tee would pollute $result,
    # making the "OK"/"FAILED" string comparison unreliable.
    echo "[$(date '+%H:%M:%S')] START  ${gene}" >> "$log"

    if pgx-run.sh "$gene" "$BAM" \
            --ref "$REF" \
            --output "${OUTPUT}/${gene}" \
            >> "$log" 2>&1; then
        echo "[$(date '+%H:%M:%S')] DONE   ${gene}" >> "$log"
        echo "OK"
    else
        echo "[$(date '+%H:%M:%S')] FAILED ${gene}" >> "$log"
        echo "FAILED"
    fi
}
export -f run_gene
export BAM REF OUTPUT LOG_DIR

# ── Parallel execution with a simple job-pool ─────────────────────────────────
declare -A PIDS STATUS
declare -a RUNNING=()
DONE=0

for GENE in "${GENES[@]}"; do
    # Wait if the pool is full
    while [[ ${#RUNNING[@]} -ge $JOBS ]]; do
        NEW_RUNNING=()
        for pid_gene in "${RUNNING[@]}"; do
            pid="${pid_gene%%:*}"
            g="${pid_gene##*:}"
            if kill -0 "$pid" 2>/dev/null; then
                NEW_RUNNING+=("$pid_gene")
            else
                wait "$pid"
                rc=$?
                if [[ $rc -eq 0 ]]; then STATUS[$g]="OK"; else STATUS[$g]="FAILED"; fi
                DONE=$(( DONE + 1 ))
                echo "    [${DONE}/${TOTAL}] ${g}: ${STATUS[$g]}"
            fi
        done
        RUNNING=("${NEW_RUNNING[@]+"${NEW_RUNNING[@]}"}")
        [[ ${#RUNNING[@]} -ge $JOBS ]] && sleep 5
    done

    # Launch gene
    echo "  --> Queuing ${GENE}"
    (
        result=$(run_gene "$GENE")
        exit $([ "$result" = "OK" ] && echo 0 || echo 1)
    ) &
    PIDS[$GENE]=$!
    RUNNING+=("${PIDS[$GENE]}:${GENE}")
done

# Wait for remaining jobs
for pid_gene in "${RUNNING[@]+"${RUNNING[@]}"}"; do
    pid="${pid_gene%%:*}"
    g="${pid_gene##*:}"
    wait "$pid"
    rc=$?
    if [[ $rc -eq 0 ]]; then STATUS[$g]="OK"; else STATUS[$g]="FAILED"; fi
    DONE=$(( DONE + 1 ))
    echo "    [${DONE}/${TOTAL}] ${g}: ${STATUS[$g]}"
done

# ── Aggregate concordance summary ─────────────────────────────────────────────
SUMMARY_TSV="${OUTPUT}/all_genes_summary.tsv"

echo -e "Gene\tStatus\tTool\tDiplotype\tActivityScore\tPhenotype\tSVMode" \
    > "$SUMMARY_TSV"

OK_COUNT=0
FAIL_COUNT=0

echo ""
echo "========================================================================"
printf " %-10s  %-6s  %-30s  %s\n" "Gene" "Status" "Top diplotype (concordant)" "Log"
echo "------------------------------------------------------------------------"

for GENE in "${GENES[@]}"; do
    TSV_FILE=$(ls "${OUTPUT}/${GENE}/${GENE}_${SAMPLE}_comparison.tsv" 2>/dev/null \
               || ls "${OUTPUT}/${GENE}"/*_comparison.tsv 2>/dev/null | head -1 \
               || true)

    if [[ "${STATUS[$GENE]:-FAILED}" == "OK" && -f "$TSV_FILE" ]]; then
        # Extract concordant diplotype (most frequent in Diplotype column)
        TOP=$(tail -n +2 "$TSV_FILE" \
              | awk -F'\t' '{print $2}' \
              | sort | uniq -c | sort -rn | head -1 \
              | awk '{print $2}')
        printf " %-10s  %-6s  %-30s  %s\n" \
            "$GENE" "OK" "${TOP:--}" "${LOG_DIR}/${GENE}.log"
        # Append rows to summary TSV
        tail -n +2 "$TSV_FILE" | while IFS=$'\t' read -r tool dip score pheno svmode rest; do
            echo -e "${GENE}\tOK\t${tool}\t${dip}\t${score}\t${pheno}\t${svmode}"
        done >> "$SUMMARY_TSV"
        OK_COUNT=$(( OK_COUNT + 1 ))
    else
        printf " %-10s  %-6s  %-30s  %s\n" \
            "$GENE" "FAILED" "-" "${LOG_DIR}/${GENE}.log"
        echo -e "${GENE}\tFAILED\t-\t-\t-\t-\t-" >> "$SUMMARY_TSV"
        FAIL_COUNT=$(( FAIL_COUNT + 1 ))
    fi
done

echo "========================================================================"
echo ""
echo "Genes completed: ${OK_COUNT}/${TOTAL}  |  Failed: ${FAIL_COUNT}/${TOTAL}"
echo "Summary TSV:     ${SUMMARY_TSV}"
echo "Per-gene logs:   ${LOG_DIR}/"
echo "Finished:        $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

[[ $FAIL_COUNT -eq 0 ]] && exit 0 || exit 1
