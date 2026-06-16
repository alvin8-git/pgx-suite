#!/usr/bin/env bash
# docker/pgx-run.sh — Phase 2 orchestration: run all supported PGx callers for a single gene
#
# Usage: pgx-run.sh <GENE> <BAM_FILE> [--ref /pgx/ref/hg38.fa] \
#                   [--output /pgx/results] [--sequential]
#
# Execution model (default: parallel):
#
#   Phase 1 — all launched simultaneously (read BAM independently):
#     bcftools mpileup → VCF
#     Aldy genotype
#     StellarPGx nextflow
#     PyPGx prepare-depth-of-coverage  [SV genes only]
#     PyPGx compute-control-statistics [SV genes only]
#     Stargazer bam2gdf                [CYP2A6/2B6/2D6 only]
#
#   Phase 2 — launched once their prerequisites are ready:
#     PyPGx run-ngs-pipeline    (waits for VCF + depth/ctrl zips)
#     Stargazer genotyping      (waits for VCF + GDF)
#
#   Phase 3:
#     pgx-compare.py            (waits for all tools)
#
# Use --sequential to disable parallelism (useful for debugging / low-RAM hosts).
#
# SV handling notes:
#   PyPGx:      SV genes (CYP2A6/2B6/2D6/2E1/4F2/G6PD/GSTM1/GSTT1) need
#               prepare-depth-of-coverage + compute-control-statistics (VDR).
#   Stargazer:  Paralog genes (CYP2A6/2B6/2D6) need a GDF depth profile from
#               the BAM. Falls back to VCF-only mode if GDF creation fails.
#   Aldy:       CN/SV handled automatically by the ILP solver.
#   StellarPGx: SVs detected natively by graphtyper.

set -uo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
REF="/pgx/ref/hg38.fa"
OUTPUT="/pgx/results"
SEQUENTIAL=0
MIN_DEPTH=10   # mean read depth below which the gene is reported NO_CALL

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF
Usage: pgx-run.sh <GENE> <BAM_FILE> [options]

Options:
  --ref PATH       GRCh38 reference FASTA (default: /pgx/ref/hg38.fa)
  --output PATH    Output directory (default: /pgx/results)
  --sequential     Disable parallel execution (for debugging)
  --min-depth N    Min mean depth to attempt a call; below this → NO_CALL (default: 10)
  -h, --help       Show this help

Examples:
  pgx-run.sh CYP2D6  /pgx/data/sample.bam
  pgx-run.sh CYP2C19 /pgx/data/sample.bam --output /pgx/results/cyp2c19
EOF
}

# ── Argument parsing ──────────────────────────────────────────────────────────
if [[ $# -lt 2 ]]; then usage; exit 1; fi

GENE="${1^^}"   # force uppercase
BAM="$2"
shift 2

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)        REF="$2";    shift 2 ;;
        --output)     OUTPUT="$2"; shift 2 ;;
        --sequential) SEQUENTIAL=1; shift ;;
        --min-depth)  MIN_DEPTH="$2"; shift 2 ;;
        -h|--help)    usage; exit 0 ;;
        *) echo "ERROR: Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

# ── Gene config — single source of truth (docker/genes.tsv) ───────────────────
# Parsed into the same associative arrays the rest of this script uses. genes.tsv
# is ALSO read by pgx-compare.py, so the support matrix and coordinates live in
# ONE file instead of being duplicated here (bash) and there (Python).
# Columns: gene region pypgx stargazer aldy stellarpgx optitype mutserve vcf_check pypgx_sv stargazer_sv stargazer_control
GENES_TSV="${PGX_GENES_TSV:-/opt/pgx/genes.tsv}"
if [[ ! -f "$GENES_TSV" ]]; then
    _self="$(readlink -f "${BASH_SOURCE[0]}" 2>/dev/null || echo "${BASH_SOURCE[0]}")"
    GENES_TSV="$(dirname "$_self")/genes.tsv"
fi
if [[ ! -f "$GENES_TSV" ]]; then
    echo "ERROR: gene config not found (looked for /opt/pgx/genes.tsv and beside the script)." >&2
    exit 1
fi

declare -A GENE_COORDS GENE_SUPPORT STARGAZER_CONTROL
PYPGX_SV_GENES=()
STARGAZER_SV_GENES=()
while IFS=$'\t' read -r _g _region _pypgx _stargazer _aldy _stellarpgx \
                        _optitype _mutserve _vcf_check _pypgx_sv _stargazer_sv _control \
                        || [[ -n "$_g" ]]; do
    [[ "$_g" == "gene" || -z "$_g" ]] && continue   # skip header / blank lines
    GENE_COORDS[$_g]="$_region"
    GENE_SUPPORT[$_g]="$_pypgx $_stargazer $_aldy $_stellarpgx"
    [[ "$_pypgx_sv"     == "1" ]] && PYPGX_SV_GENES+=("$_g")
    [[ "$_stargazer_sv" == "1" ]] && STARGAZER_SV_GENES+=("$_g")
    [[ "$_control" != "-" && -n "$_control" ]] && STARGAZER_CONTROL[$_g]="$_control"
done < "$GENES_TSV"

# ── Helpers ───────────────────────────────────────────────────────────────────
in_array() {
    local needle="$1"; shift
    for item in "$@"; do [[ "$item" == "$needle" ]] && return 0; done
    return 1
}

log_status() {
    # Print a timestamped status line to stdout
    echo "[$(date '+%H:%M:%S')] $*"
}

GENE_LOWER="${GENE,,}"

echo "============================================================"
echo " PGx Suite — Orchestration Runner"
echo " Gene:   ${GENE}"
echo " BAM:    ${BAM}"
echo " Ref:    ${REF}"
echo " Output: ${OUTPUT}"
echo " Build:  GRCh38"
echo " Mode:   $([[ $SEQUENTIAL -eq 1 ]] && echo sequential || echo parallel)"
echo "============================================================"
echo ""

# ── Input validation ──────────────────────────────────────────────────────────
log_status "Validating inputs..."

if [[ -z "${GENE_SUPPORT[$GENE]:-}" ]]; then
    echo "ERROR: Gene '${GENE}' is not in the supported gene list." >&2
    echo "Supported: ${!GENE_SUPPORT[*]}" >&2
    exit 1
fi
if [[ ! -f "$BAM" ]]; then
    echo "ERROR: BAM file not found: $BAM" >&2; exit 1
fi
_INPUT_EXT="${BAM##*.}"
case "$_INPUT_EXT" in
    bam)  _IDX_EXT="bai" ;;
    cram) _IDX_EXT="crai" ;;
    *)    echo "ERROR: Input must be .bam or .cram: $BAM" >&2; exit 1 ;;
esac
BAI="${BAM}.${_IDX_EXT}"
[[ ! -f "$BAI" ]] && BAI="${BAM%.*}.${_IDX_EXT}"
if [[ ! -f "$BAI" ]]; then
    echo "ERROR: Index not found (${BAM}.${_IDX_EXT}). Run: samtools index $BAM" >&2; exit 1
fi
if [[ ! -f "$REF" ]]; then
    echo "ERROR: Reference FASTA not found: $REF" >&2; exit 1
fi
if [[ ! -f "${REF}.fai" ]]; then
    echo "ERROR: Reference index not found: ${REF}.fai" >&2; exit 1
fi

echo "  Gene / BAM / Ref — OK"
echo ""

# ── Sample name from BAM read group ──────────────────────────────────────────
SAMPLE=$(samtools view -H "$BAM" | grep '^@RG' | grep -oP 'SM:\K[^\t]+' | head -1 || true)
[[ -z "$SAMPLE" ]] && SAMPLE=$(basename "${BAM%.*}")
log_status "Sample: ${SAMPLE}"
echo ""

# ── Parse tool support flags ──────────────────────────────────────────────────
read -r DO_PYPGX DO_STARGAZER DO_ALDY DO_STELLARPGX <<< "${GENE_SUPPORT[$GENE]}"

IS_PYPGX_SV=0
in_array "$GENE" "${PYPGX_SV_GENES[@]}"    && IS_PYPGX_SV=1

IS_STARGAZER_SV=0
in_array "$GENE" "${STARGAZER_SV_GENES[@]}" && IS_STARGAZER_SV=1

if [[ "$DO_STELLARPGX" -eq 1 ]] && \
   [[ ! -d "/pgx/stellarpgx" || ! -f "/pgx/containers/stellarpgx-dev.sif" ]]; then
    log_status "WARN StellarPGx volumes not mounted — skipping"
    DO_STELLARPGX=0
fi

# ── OptiType: HLA-A and HLA-B genes bypass the standard tool pipeline ─────────
DO_OPTITYPE=0
if [[ "$GENE" =~ ^HLA- ]]; then
    DO_OPTITYPE=1
    DO_PYPGX=0; DO_STARGAZER=0; DO_ALDY=0; DO_STELLARPGX=0
    IS_PYPGX_SV=0; IS_STARGAZER_SV=0
fi
if [[ "$DO_OPTITYPE" -eq 1 ]] && [[ ! -f "/pgx/containers/optitype.sif" ]]; then
    log_status "WARN OptiType SIF not found at /pgx/containers/optitype.sif — skipping HLA typing"
    DO_OPTITYPE=0
fi

# ── mutserve: MT-RNR1 bypasses the standard tool pipeline ─────────────────────
DO_MUTSERVE=0
if [[ "$GENE" == "MT-RNR1" ]]; then
    DO_MUTSERVE=1
    DO_PYPGX=0; DO_STARGAZER=0; DO_ALDY=0; DO_STELLARPGX=0
    IS_PYPGX_SV=0; IS_STARGAZER_SV=0
fi
if [[ "$DO_MUTSERVE" -eq 1 ]] && [[ ! -f "/usr/local/bin/mutserve.jar" ]]; then
    log_status "WARN mutserve.jar not found — skipping MT-RNR1 calling"
    DO_MUTSERVE=0
fi

echo "  PyPGx:      $([[ $DO_PYPGX      -eq 1 ]] && echo YES || echo NO)  (SV preprocessing: $([[ $IS_PYPGX_SV    -eq 1 ]] && echo YES || echo NO))"
echo "  Stargazer:  $([[ $DO_STARGAZER  -eq 1 ]] && echo YES || echo NO)  (GDF/SV mode:      $([[ $IS_STARGAZER_SV -eq 1 ]] && echo YES || echo NO))"
echo "  Aldy:       $([[ $DO_ALDY       -eq 1 ]] && echo YES || echo NO)  (SV: auto via ILP)"
echo "  StellarPGx: $([[ $DO_STELLARPGX -eq 1 ]] && echo YES || echo NO)  (SV: auto via graphtyper)"
echo "  OptiType:   $([[ $DO_OPTITYPE   -eq 1 ]] && echo YES || echo NO)  (HLA Class I typing)"
echo "  mutserve:   $([[ $DO_MUTSERVE   -eq 1 ]] && echo YES || echo NO)  (MT-RNR1 aminoglycoside risk)"
echo ""

# ── Create output directories ─────────────────────────────────────────────────
mkdir -p \
    "${OUTPUT}/pypgx" \
    "${OUTPUT}/stargazer" \
    "${OUTPUT}/aldy" \
    "${OUTPUT}/stellarpgx" \
    "${OUTPUT}/logs"

# Clear stale per-tool exit-code markers from any previous run of this gene.
rm -f "${OUTPUT}/logs/"*.status 2>/dev/null || true

VCF="${OUTPUT}/${GENE}.vcf.gz"
DEPTH_ZIP="${OUTPUT}/depth-of-coverage.zip"
CTRL_ZIP="${OUTPUT}/control-stats-VDR.zip"
COORDS="${GENE_COORDS[$GENE]:-}"

# ── Per-gene coverage (drives pgx-compare's NO_CALL gate) ─────────────────────
# Measured at call time over the gene region; cheap even on CRAM (single seek).
# Without it, a no-coverage region yields a confident wild-type (*1/*1) call.
MEAN_DEPTH=""
if [[ -n "$COORDS" && "$COORDS" != "ALT_CONTIG" ]]; then
    MEAN_DEPTH=$(samtools coverage -r "$COORDS" --reference "$REF" "$BAM" 2>/dev/null \
                 | awk 'NR==2 {print $7}')   # column 7 = meandepth
    if [[ -z "$MEAN_DEPTH" ]]; then
        log_status "WARN  could not measure depth for ${COORDS}; coverage gate disabled"
    else
        log_status "INFO  mean depth over ${GENE} (${COORDS}): ${MEAN_DEPTH}x (gate: ${MIN_DEPTH}x)"
    fi
fi

# ── Tool runner functions ─────────────────────────────────────────────────────
# Each writes its own log to ${OUTPUT}/logs/<tool>.log
# Returns 0 on success, non-zero on failure.

run_bcftools() {
    local log="${OUTPUT}/logs/bcftools.log"
    log_status "START  bcftools mpileup  (region: ${COORDS})"
    if bcftools mpileup \
            -r "$COORDS" \
            -f "$REF" \
            -a AD,DP \
            --max-depth 500 \
            -o - \
            "$BAM" \
        | bcftools call \
            -m -v \
            --output-type z \
            -o "$VCF" \
        >> "$log" 2>&1 \
       && tabix -p vcf "$VCF" >> "$log" 2>&1; then
        log_status "DONE   bcftools  →  ${VCF}"
        return 0
    else
        log_status "FAILED bcftools  (see ${log})"
        return 1
    fi
}

run_pypgx_sv_preprocessing() {
    local log="${OUTPUT}/logs/pypgx_sv_prep.log"
    log_status "START  PyPGx SV preprocessing  (prepare-depth-of-coverage + control-statistics VDR)"
    {
        pypgx prepare-depth-of-coverage \
            "$DEPTH_ZIP" "$BAM" --assembly GRCh38 \
        && pypgx compute-control-statistics \
            VDR "$CTRL_ZIP" "$BAM" --assembly GRCh38
    } >> "$log" 2>&1
    local rc=$?
    if [[ $rc -eq 0 ]]; then
        log_status "DONE   PyPGx SV preprocessing"
    else
        log_status "FAILED PyPGx SV preprocessing  (see ${log})"
    fi
    return $rc
}

run_stargazer_gdf() {
    local log="${OUTPUT}/logs/stargazer_gdf.log"
    # NOTE: gdf_dir must NOT be inside the stargazer output dir because Stargazer's
    # genotyping step does shutil.rmtree(output_dir) on startup, which would delete
    # any GDF placed under ${OUTPUT}/stargazer/.
    local gdf_dir="${OUTPUT}/stargazer_gdf"
    mkdir -p "$gdf_dir"
    log_status "START  Stargazer bam2gdf  (control=VDR)"
    if stargazer \
            -G "${GENE_LOWER}.gdf" \
            -t "$GENE_LOWER" \
            -c vdr \
            -B "$BAM" \
            -o "$gdf_dir" \
            -a grc38 \
            -d wgs \
        >> "$log" 2>&1; then
        log_status "DONE   Stargazer bam2gdf  →  ${gdf_dir}/${GENE_LOWER}.gdf (outside stargazer/ to survive rmtree)"
        return 0
    else
        log_status "FAILED Stargazer bam2gdf  (see ${log})"
        return 1
    fi
}

run_aldy() {
    local log="${OUTPUT}/logs/aldy.log"
    local rc=0
    log_status "START  Aldy"
    if aldy genotype \
            -g "$GENE" \
            -p illumina \
            "$BAM" \
            -o "${OUTPUT}/aldy/${GENE}.aldy" \
        >> "$log" 2>&1; then
        log_status "DONE   Aldy"
    else
        rc=$?
        log_status "FAILED Aldy  (rc=${rc}, see ${log})"
    fi
    # Failure is non-fatal to the run, but the real exit code is recorded so the
    # aggregator can tell "ran and crashed" from "never ran" (no set -e, so a
    # non-zero return does not abort the background job pool).
    echo "$rc" > "${OUTPUT}/logs/aldy.status"
    return "$rc"
}

run_stellarpgx() {
    local log="${OUTPUT}/logs/stellarpgx.log"
    local bam_dir bam_base
    bam_dir="$(dirname "$BAM")"
    bam_base="$(basename "${BAM%.*}")"   # strips .bam or .cram
    # StellarPGx main.nf uses 'cypor' not 'por' as the gene identifier.
    local stellar_gene="${GENE_LOWER}"
    local rc=0
    [[ "$stellar_gene" == "por" ]] && stellar_gene="cypor"
    log_status "START  StellarPGx"
    if nextflow run /pgx/stellarpgx/main.nf \
            --gene "${stellar_gene}" \
            --in_bam "${bam_dir}/${bam_base}*{${_INPUT_EXT},${_IDX_EXT}}" \
            --ref_file "$REF" \
            --out_dir "${OUTPUT}/stellarpgx" \
            --res_init /pgx/stellarpgx/resources \
            --db_init /pgx/stellarpgx/database \
            --caller_init /pgx/stellarpgx/scripts \
            -work-dir "${OUTPUT}/stellarpgx/.work" \
        >> "$log" 2>&1; then
        log_status "DONE   StellarPGx"
    else
        rc=$?
        log_status "FAILED StellarPGx  (rc=${rc}, see ${log})"
    fi
    echo "$rc" > "${OUTPUT}/logs/stellarpgx.status"
    return "$rc"
}

run_hla() {
    local log="${OUTPUT}/logs/hla.log"
    local rc=0
    log_status "START  OptiType HLA typing  (gene: ${GENE})"
    if pgx-hla.sh "$GENE" "$BAM" "$SAMPLE" "$OUTPUT" \
        >> "$log" 2>&1; then
        log_status "DONE   OptiType HLA typing"
    else
        rc=$?
        log_status "FAILED OptiType HLA typing  (rc=${rc}, see ${log})"
    fi
    echo "$rc" > "${OUTPUT}/logs/hla.status"
    return "$rc"
}

run_mt() {
    local log="${OUTPUT}/logs/mutserve.log"
    local rc=0
    log_status "START  mutserve MT-RNR1 calling"
    if pgx-mt.sh "$BAM" "$SAMPLE" "$OUTPUT" --ref "$REF" \
        >> "$log" 2>&1; then
        log_status "DONE   mutserve MT-RNR1 calling"
    else
        rc=$?
        log_status "FAILED mutserve MT-RNR1 calling  (rc=${rc}, see ${log})"
    fi
    echo "$rc" > "${OUTPUT}/logs/mutserve.status"
    return "$rc"
}

run_pypgx_pipeline() {
    local log="${OUTPUT}/logs/pypgx.log"
    local pypgx_args=("$GENE" "${OUTPUT}/pypgx" --assembly GRCh38 --force)
    [[ "$NEED_VCF" -eq 1 ]] && pypgx_args+=(--variants "$VCF")
    if [[ "$IS_PYPGX_SV" -eq 1 ]]; then
        pypgx_args+=(--depth-of-coverage "$DEPTH_ZIP" --control-statistics "$CTRL_ZIP")
        log_status "START  PyPGx  (SV mode)"
    else
        log_status "START  PyPGx  (standard)"
    fi
    local rc=0
    if pypgx run-ngs-pipeline "${pypgx_args[@]}" >> "$log" 2>&1; then
        log_status "DONE   PyPGx"
    else
        rc=$?
        log_status "FAILED PyPGx  (rc=${rc}, see ${log})"
    fi
    echo "$rc" > "${OUTPUT}/logs/pypgx.status"
    return "$rc"
}

run_stargazer_genotype() {
    local log="${OUTPUT}/logs/stargazer.log"
    local gdf_file="${OUTPUT}/stargazer_gdf/${GENE_LOWER}.gdf"
    local stargazer_args=(-t "$GENE_LOWER" -d wgs -a grc38 -i "$VCF" -o "${OUTPUT}/stargazer")
    local control="${STARGAZER_CONTROL[$GENE]:-}"

    if [[ "$IS_STARGAZER_SV" -eq 1 && -f "$gdf_file" ]]; then
        stargazer_args+=(-c vdr -g "$gdf_file")
        log_status "START  Stargazer  (SV/GDF mode)"
    else
        [[ -n "$control" ]] && stargazer_args+=(-c "$control")
        log_status "START  Stargazer  (VCF-only mode)"
    fi
    local rc=0
    if stargazer "${stargazer_args[@]}" >> "$log" 2>&1; then
        log_status "DONE   Stargazer"
    else
        rc=$?
        log_status "FAILED Stargazer  (rc=${rc}, see ${log})"
    fi
    echo "$rc" > "${OUTPUT}/logs/stargazer.status"
    return "$rc"
}

# ── VCF availability check ────────────────────────────────────────────────────
NEED_VCF=$(( (DO_PYPGX + DO_STARGAZER) > 0 ? 1 : 0 ))
if [[ "$COORDS" == "ALT_CONTIG" ]]; then
    log_status "INFO   ${GENE} is on an alt contig — bcftools skipped; VCF-dependent tools may be limited"
    NEED_VCF=0
elif [[ -z "$COORDS" ]]; then
    echo "ERROR: No GRCh38 coordinates defined for ${GENE}" >&2; exit 1
fi

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 1 — Launch all BAM-reading tasks simultaneously
# ─────────────────────────────────────────────────────────────────────────────
echo "------------------------------------------------------------"
log_status "Phase 1: launching all BAM-reading tasks"
echo "------------------------------------------------------------"

BCFTOOLS_PID=""
DEPTH_PID=""
GDF_PID=""
ALDY_PID=""
STELLARPGX_PID=""
HLA_PID=""
MT_PID=""

if [[ "$NEED_VCF" -eq 1 ]]; then
    if [[ "$SEQUENTIAL" -eq 0 ]]; then
        run_bcftools & BCFTOOLS_PID=$!
    else
        run_bcftools
    fi
fi

if [[ "$DO_PYPGX" -eq 1 && "$IS_PYPGX_SV" -eq 1 ]]; then
    if [[ "$SEQUENTIAL" -eq 0 ]]; then
        run_pypgx_sv_preprocessing & DEPTH_PID=$!
    else
        run_pypgx_sv_preprocessing
    fi
fi

if [[ "$DO_STARGAZER" -eq 1 && "$IS_STARGAZER_SV" -eq 1 ]]; then
    if [[ "$SEQUENTIAL" -eq 0 ]]; then
        run_stargazer_gdf & GDF_PID=$!
    else
        run_stargazer_gdf || IS_STARGAZER_SV=0
    fi
fi

if [[ "$DO_ALDY" -eq 1 ]]; then
    if [[ "$SEQUENTIAL" -eq 0 ]]; then
        run_aldy & ALDY_PID=$!
    else
        run_aldy
    fi
fi

if [[ "$DO_STELLARPGX" -eq 1 ]]; then
    if [[ "$SEQUENTIAL" -eq 0 ]]; then
        run_stellarpgx & STELLARPGX_PID=$!
    else
        run_stellarpgx
    fi
fi

if [[ "$DO_OPTITYPE" -eq 1 ]]; then
    if [[ "$SEQUENTIAL" -eq 0 ]]; then
        run_hla & HLA_PID=$!
    else
        run_hla
    fi
fi

if [[ "$DO_MUTSERVE" -eq 1 ]]; then
    if [[ "$SEQUENTIAL" -eq 0 ]]; then
        run_mt & MT_PID=$!
    else
        run_mt
    fi
fi

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 2 — Wait for prerequisites, then launch VCF-dependent tools
# ─────────────────────────────────────────────────────────────────────────────
echo "------------------------------------------------------------"
log_status "Phase 2: waiting for VCF and SV prerequisites"
echo "------------------------------------------------------------"

PYPGX_PID=""
STARGAZER_PID=""

if [[ "$SEQUENTIAL" -eq 0 ]]; then
    # Wait for bcftools before launching VCF-dependent tools
    if [[ -n "$BCFTOOLS_PID" ]]; then
        wait "$BCFTOOLS_PID" || log_status "WARN  bcftools exited non-zero — VCF may be missing"
    fi

    # PyPGx: wait for SV preprocessing to also complete, then launch pipeline
    if [[ "$DO_PYPGX" -eq 1 ]]; then
        [[ -n "$DEPTH_PID" ]] && wait "$DEPTH_PID" || true
        run_pypgx_pipeline & PYPGX_PID=$!
    fi

    # Stargazer: wait for GDF to also complete, then launch genotyping
    if [[ "$DO_STARGAZER" -eq 1 && "$NEED_VCF" -eq 1 ]]; then
        if [[ -n "$GDF_PID" ]]; then
            wait "$GDF_PID" || { log_status "WARN  bam2gdf failed — falling back to VCF-only mode"; IS_STARGAZER_SV=0; }
        fi
        run_stargazer_genotype & STARGAZER_PID=$!
    fi
else
    # Sequential fallback
    [[ "$DO_PYPGX"      -eq 1 ]]                          && run_pypgx_pipeline
    [[ "$DO_STARGAZER"  -eq 1 && "$NEED_VCF" -eq 1 ]]    && run_stargazer_genotype
fi

# ─────────────────────────────────────────────────────────────────────────────
# Wait for all remaining background jobs
# ─────────────────────────────────────────────────────────────────────────────
echo "------------------------------------------------------------"
log_status "Waiting for all tools to finish..."
echo "------------------------------------------------------------"

[[ -n "$PYPGX_PID"      ]] && wait "$PYPGX_PID"      || true
[[ -n "$STARGAZER_PID"  ]] && wait "$STARGAZER_PID"  || true
[[ -n "$ALDY_PID"       ]] && wait "$ALDY_PID"       || true
[[ -n "$STELLARPGX_PID" ]] && wait "$STELLARPGX_PID" || true
[[ -n "$HLA_PID"        ]] && wait "$HLA_PID"        || true
[[ -n "$MT_PID"         ]] && wait "$MT_PID"         || true

echo ""
log_status "All tools finished."

# ── Tool execution status (real exit codes, no longer masked as success) ──────
# Each runner writes <tool>.status with its exit code; a missing file means the
# tool was not selected for this gene. Failures are reported to pgx-compare so a
# crashed caller is marked 'failed', not the misleading 'not_run'.
FAILED_TOOLS=()
echo "  Tool execution status:"
for _pair in "PyPGx:pypgx" "Stargazer:stargazer" "Aldy:aldy" \
             "StellarPGx:stellarpgx" "OptiType:hla" "mutserve:mutserve"; do
    _name="${_pair%%:*}"; _key="${_pair##*:}"
    _sf="${OUTPUT}/logs/${_key}.status"
    [[ -f "$_sf" ]] || continue
    _rc="$(cat "$_sf" 2>/dev/null || echo '?')"
    if [[ "$_rc" == "0" ]]; then
        echo "    ${_name}: OK"
    else
        echo "    ${_name}: FAILED (rc=${_rc})"
        FAILED_TOOLS+=("$_name")
    fi
done
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 3 — Comparison report
# ─────────────────────────────────────────────────────────────────────────────
echo "------------------------------------------------------------"
log_status "Phase 3: generating comparison report"
echo "------------------------------------------------------------"
COMPARE_ARGS=(--gene "$GENE" --sample "$SAMPLE" --output-dir "$OUTPUT" --min-depth "$MIN_DEPTH")
[[ -n "$MEAN_DEPTH" ]] && COMPARE_ARGS+=(--region-depth "$MEAN_DEPTH")
if [[ ${#FAILED_TOOLS[@]} -gt 0 ]]; then
    COMPARE_ARGS+=(--failed-tools "$(IFS=,; echo "${FAILED_TOOLS[*]}")")
fi
python3 /opt/pgx/pgx-compare.py "${COMPARE_ARGS[@]}"
