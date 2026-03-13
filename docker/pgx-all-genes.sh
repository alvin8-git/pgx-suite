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

# ── Validate input alignment file (BAM or CRAM) ───────────────────────────────
if [[ ! -f "$BAM" ]]; then
    echo "ERROR: Input file not found: $BAM" >&2; exit 1
fi
_INPUT_EXT="${BAM##*.}"
case "$_INPUT_EXT" in
    bam)  _IDX_EXT="bai" ;;
    cram) _IDX_EXT="crai" ;;
    *)    echo "ERROR: Input must be .bam or .cram: $BAM" >&2; exit 1 ;;
esac
if [[ ! -f "${BAM}.${_IDX_EXT}" && ! -f "${BAM%.*}.${_IDX_EXT}" ]]; then
    echo "ERROR: Index not found (expected ${BAM}.${_IDX_EXT})" >&2; exit 1
fi
if [[ ! -f "$REF" ]]; then
    echo "ERROR: Reference FASTA not found: $REF" >&2; exit 1
fi

# ── Gene list — all unique supported genes ────────────────────────────────────
# POR/CYPOR are the same locus; only POR is kept here.
GENES=(
    ABCG2
    CACNA1S
    CYP1A1 CYP1A2
    CYP2A6 CYP2B6 CYP2C8 CYP2C9 CYP2C19 CYP2D6 CYP2E1
    CYP3A4 CYP3A5 CYP4F2
    DPYD
    G6PD GSTM1 GSTT1
    HLA-A HLA-B
    IFNL3
    MT-RNR1
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
LOG_DIR="${OUTPUT}/log"
GENES_DIR="${OUTPUT}/Genes"
mkdir -p "$LOG_DIR" "$GENES_DIR"

# ── CRAM → PGx-region BAM conversion ─────────────────────────────────────────
# Aldy, PyPGx, and Stargazer use pysam/internal samtools without an explicit
# --reference flag, so they cannot decode CRAM inside the container (the CRAM
# header's UR: path points to the host filesystem, not /pgx/ref/).
# Solution: extract all PGx gene regions + VDR control + MHC + full chrM into
# a single BAM once per sample run, then hand that BAM to all per-gene calls.
#
# _ORIG_BAM is kept pointing to the original input (CRAM or BAM) so that
# pgx-bamstats.sh can compute whole-genome QC from the full file, not just the
# extracted subset (total reads, genome depth, dup%, mapq% would all be wrong
# if computed from the small region-extracted BAM).
_ORIG_BAM="$BAM"
if [[ "$_INPUT_EXT" == "cram" ]]; then
    _PGX_BAM="${LOG_DIR}/pgx_input.bam"
    # Cap at 8 threads: BGZF block ordering corruption observed with high counts (>16)
    _NCPU=$(( $(nproc 2>/dev/null || echo 4) > 8 ? 8 : $(nproc 2>/dev/null || echo 4) ))
    if [[ ! -f "${_PGX_BAM}.bai" ]]; then
        rm -f "$_PGX_BAM"   # remove any leftover corrupt BAM from a previous failed run
        echo "INFO: CRAM input detected — extracting PGx regions to BAM …"
        echo "      Reference : ${REF}"
        echo "      Threads   : ${_NCPU}"
        echo "      (this is a one-time step; all per-gene calls will reuse ${_PGX_BAM})"

        # Identify the reference the CRAM was encoded with (shown on mismatch errors).
        _CRAM_REF_URI=$(samtools view -H "$BAM" 2>/dev/null \
            | awk '/^@SQ/{for(i=1;i<=NF;i++) if($i~/^UR:/) {sub(/^UR:/,"",$i); print $i; exit}}')
        [[ -n "$_CRAM_REF_URI" ]] && echo "      CRAM encoded with: ${_CRAM_REF_URI}"

        if ! samtools view -b -@ "$_NCPU" -T "$REF" \
                -L /opt/pgx/pgx_cram_regions.bed \
                -o "$_PGX_BAM" "$BAM"; then
            echo "" >&2
            echo "ERROR: CRAM → BAM conversion failed." >&2
            echo "       Most likely cause: reference MD5 mismatch." >&2
            echo "       The CRAM was encoded with a different hg38 build than" >&2
            echo "       the reference currently mounted at: ${REF}" >&2
            if [[ -n "$_CRAM_REF_URI" ]]; then
                echo "" >&2
                echo "       CRAM declares its reference as: ${_CRAM_REF_URI}" >&2
                echo "       Re-run with: ./run_pgx_suite.sh <cram> --ref /path/to/matching/hg38.fa" >&2
            fi
            rm -f "$_PGX_BAM"
            exit 1
        fi

        # Verify BAM integrity before indexing — catches BGZF truncation from write errors.
        if ! samtools quickcheck "$_PGX_BAM" 2>&1; then
            echo "ERROR: CRAM → BAM conversion produced a corrupt BAM (quickcheck failed)." >&2
            echo "       Check disk space (df -h) and re-run." >&2
            rm -f "$_PGX_BAM"
            exit 1
        fi

        if ! samtools index -@ "$_NCPU" "$_PGX_BAM"; then
            echo "ERROR: samtools index failed on ${_PGX_BAM}" >&2
            rm -f "$_PGX_BAM" "${_PGX_BAM}.bai"
            exit 1
        fi
        echo "INFO: CRAM conversion complete → ${_PGX_BAM}"
    else
        echo "INFO: Reusing existing PGx-region BAM: ${_PGX_BAM}"
    fi
    BAM="$_PGX_BAM"
fi

# ── Extract sample name from the original input filename (before the first ".") ─
# Use _ORIG_BAM so that CRAM inputs give the real sample name, not "pgx_input".
# e.g. C-5839.bwa.sortdup.bqsr.cram → C-5839
SAMPLE=$(basename "$_ORIG_BAM" | cut -d'.' -f1)

# ── BAM QC — launched in background, runs in parallel with gene callers ───────
# pgx-bamstats.sh is I/O bound (~6 min on 167 GB WGS at 490 MB/s disk).
# Running it concurrently with gene calling keeps total wall time at
# max(bamstats_time, gene_calling_time) instead of their sum.
# BAMSTATS_PID is waited on after the gene loop, before HTML report generation.
#
# Pass _ORIG_BAM (the full CRAM or BAM), not the extracted pgx_input.bam.
# pgx-bamstats.sh supports CRAM natively — it adds -T/--fasta reference flags
# when the input extension is .cram, so all samtools/mosdepth calls decode
# correctly without a full BAM conversion.
echo "Starting BAM QC in background …"
pgx-bamstats.sh "$_ORIG_BAM" "$SAMPLE" "${LOG_DIR}" "$REF" \
    > "${LOG_DIR}/bamstats.log" 2>&1 &
BAMSTATS_PID=$!

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
            --output "${GENES_DIR}/${gene}" \
            >> "$log" 2>&1; then
        echo "[$(date '+%H:%M:%S')] DONE   ${gene}" >> "$log"
        echo "OK"
    else
        echo "[$(date '+%H:%M:%S')] FAILED ${gene}" >> "$log"
        echo "FAILED"
    fi
}
export -f run_gene
export BAM REF OUTPUT LOG_DIR GENES_DIR

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
SUMMARY_TSV="${LOG_DIR}/all_genes_summary.tsv"

echo -e "Gene\tStatus\tTool\tDiplotype\tActivityScore\tPhenotype\tSVMode" \
    > "$SUMMARY_TSV"

OK_COUNT=0
FAIL_COUNT=0

echo ""
echo "========================================================================"
printf " %-10s  %-6s  %-30s  %s\n" "Gene" "Status" "Top diplotype (concordant)" "Log"
echo "------------------------------------------------------------------------"

for GENE in "${GENES[@]}"; do
    TSV_FILE=$(ls "${GENES_DIR}/${GENE}/${GENE}_${SAMPLE}_comparison.tsv" 2>/dev/null \
               || ls "${GENES_DIR}/${GENE}"/*_comparison.tsv 2>/dev/null | head -1 \
               || true)

    if [[ "${STATUS[$GENE]:-FAILED}" == "OK" && -f "$TSV_FILE" ]]; then
        # Comparison TSV columns: Gene Sample Build Tool Diplotype ActivityScore Phenotype Status SVMode
        # Extract concordant diplotype (most frequent value in Diplotype column = col 5)
        TOP=$(tail -n +2 "$TSV_FILE" \
              | awk -F'\t' '$5 != "-" {print $5}' \
              | sort | uniq -c | sort -rn | head -1 \
              | awk '{print $2}')
        printf " %-10s  %-6s  %-30s  %s\n" \
            "$GENE" "OK" "${TOP:--}" "${LOG_DIR}/${GENE}.log"
        # Append rows to summary TSV
        tail -n +2 "$TSV_FILE" | while IFS=$'\t' read -r _gene _sample _build tool dip score pheno _status svmode; do
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

# ── Wait for background BAM QC to finish ─────────────────────────────────────
echo "Waiting for BAM QC to complete …"
if wait "$BAMSTATS_PID"; then
    echo "BAM QC: OK"
else
    echo "BAM QC: WARN (bam_stats.json may be missing or incomplete)"
fi

# ── HTML report ───────────────────────────────────────────────────────────────
echo "Generating HTML report …"
pgx-report.py \
    --sample    "$SAMPLE" \
    --output    "$OUTPUT" \
    --genes-dir "${GENES_DIR}" \
    --bam       "$BAM" \
    --bam-stats "${LOG_DIR}/bam_stats.json" \
    && echo "HTML report: ${OUTPUT}/${SAMPLE}_pgx_report.html" \
    || echo "WARN: HTML report generation failed (results still available in TSV)"

[[ $FAIL_COUNT -eq 0 ]] && exit 0 || exit 1
