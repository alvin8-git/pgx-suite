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

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF
Usage: pgx-run.sh <GENE> <BAM_FILE> [options]

Options:
  --ref PATH       GRCh38 reference FASTA (default: /pgx/ref/hg38.fa)
  --output PATH    Output directory (default: /pgx/results)
  --sequential     Disable parallel execution (for debugging)
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
        -h|--help)    usage; exit 0 ;;
        *) echo "ERROR: Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

# ── GRCh38 gene coordinate table ──────────────────────────────────────────────
# Sourced from pypgx gene-table.csv (authoritative regions with upstream/downstream buffer).
# GSTT1 is on chr22_KI270879v1_alt — bcftools mpileup is skipped for it.
declare -A GENE_COORDS
GENE_COORDS=(
    [CYP1A1]="chr15:74716541-74728528"
    [CYP1A2]="chr15:74745844-74759607"
    [CYP2A6]="chr19:40833540-40890447"
    [CYP2B6]="chr19:40921281-41028398"
    [CYP2C8]="chr10:95033771-95072497"
    [CYP2C9]="chr10:94935657-94993091"
    [CYP2C19]="chr10:94759680-94858547"
    [CYP2D6]="chr22:42116498-42155810"
    [CYP2E1]="chr10:133517362-133549123"
    [CYP3A4]="chr7:99753966-99787184"
    [CYP3A5]="chr7:99645193-99682996"
    [CYP4F2]="chr19:15863022-15913074"
    [DPYD]="chr1:97074742-97924034"
    [G6PD]="chrX:154528389-154550018"
    [GSTM1]="chr1:109684816-109696745"
    [GSTT1]="ALT_CONTIG"   # chr22_KI270879v1_alt — skip bcftools
    [IFNL3]="chr19:39240552-39253525"
    [NAT1]="chr8:18207108-18226689"
    [NAT2]="chr8:18388281-18404218"
    [NUDT15]="chr13:48034725-48050221"
    [POR]="chr7:75912154-75989855"
    [CYPOR]="chr7:75912154-75989855"
    [RYR1]="chr19:38430690-38590564"
    [SLCO1B1]="chr12:21128193-21242796"
    [TPMT]="chr6:18125310-18158169"
    [UGT1A1]="chr2:233754269-233779300"
    [VKORC1]="chr16:31087853-31097797"
    [ABCG2]="chr4:88085265-88236626"
    [HLA-A]="chr6:28510020-33480577"   # MHC region — used for read extraction by pgx-hla.sh
    [HLA-B]="chr6:28510020-33480577"   # MHC region — same extraction as HLA-A
    [CACNA1S]="chr1:201006956-201083927"
    [MT-RNR1]="chrM:648-1601"          # mitochondrial 12S rRNA; mutserve uses full chrM
)

# ── Gene support matrix (pypgx stargazer aldy stellarpgx) ─────────────────────
# 1=supported, 0=not supported
declare -A GENE_SUPPORT
GENE_SUPPORT=(
    [CYP2D6]="1 1 1 1"
    [CYP2C19]="1 1 1 1"
    [CYP2C9]="1 1 1 1"
    [CYP2B6]="1 1 1 1"
    [CYP2C8]="1 1 1 1"
    [CYP3A4]="1 1 1 1"
    [CYP3A5]="1 1 1 1"
    [CYP4F2]="1 1 1 1"
    [NUDT15]="1 1 1 1"
    [TPMT]="1 1 1 1"
    [UGT1A1]="1 1 1 1"
    [SLCO1B1]="1 1 1 1"
    [DPYD]="1 1 1 0"
    [NAT1]="1 1 1 1"
    [NAT2]="1 1 1 1"
    [G6PD]="1 1 1 0"
    [GSTM1]="1 1 1 1"
    [GSTT1]="1 0 0 1"
    [POR]="1 1 0 1"
    [CYPOR]="1 1 0 1"
    [VKORC1]="1 1 1 0"
    [CYP1A1]="1 1 1 1"
    [CYP1A2]="1 1 1 1"
    [CYP2A6]="1 1 1 1"
    [CYP2E1]="1 1 1 1"
    [IFNL3]="1 1 1 0"
    [RYR1]="1 1 1 0"
    [ABCG2]="0 0 1 1"
    [HLA-A]="0 0 0 0"   # OptiType only — handled by run_hla() special case
    [HLA-B]="0 0 0 0"   # OptiType only — handled by run_hla() special case
    [CACNA1S]="1 0 0 1" # PyPGx + StellarPGx; Stargazer/Aldy do not support CACNA1S
    [MT-RNR1]="0 0 0 0" # mutserve only — handled by run_mt() special case
)

# ── PyPGx SV genes — need depth-of-coverage + control-statistics (VDR) ───────
PYPGX_SV_GENES=(CYP2A6 CYP2B6 CYP2D6 CYP2E1 CYP4F2 G6PD GSTM1 GSTT1)

# ── Stargazer SV genes — need GDF depth profile created from BAM ──────────────
STARGAZER_SV_GENES=(CYP2A6 CYP2B6 CYP2D6)

# ── Stargazer control gene mapping ────────────────────────────────────────────
declare -A STARGAZER_CONTROL
STARGAZER_CONTROL=(
    [CYP2D6]="vdr"   [CYP2C19]="vdr"  [CYP2C9]="vdr"   [CYP2B6]="vdr"
    [CYP2C8]="vdr"   [CYP3A4]="vdr"   [CYP3A5]="vdr"   [CYP4F2]="vdr"
    [CYP1A1]="vdr"   [CYP1A2]="vdr"   [CYP2A6]="vdr"   [CYP2E1]="vdr"
    [SLCO1B1]="vdr"  [G6PD]="vdr"     [VKORC1]="vdr"
    [NUDT15]="vdr"   [TPMT]="vdr"     [UGT1A1]="vdr"   [DPYD]="vdr"
    [NAT1]="vdr"     [NAT2]="vdr"     [GSTM1]="vdr"    [IFNL3]="vdr"
    [RYR1]="vdr"     [POR]="vdr"
)

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

VCF="${OUTPUT}/${GENE}.vcf.gz"
DEPTH_ZIP="${OUTPUT}/depth-of-coverage.zip"
CTRL_ZIP="${OUTPUT}/control-stats-VDR.zip"
COORDS="${GENE_COORDS[$GENE]:-}"

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
    log_status "START  Aldy"
    if aldy genotype \
            -g "$GENE" \
            -p illumina \
            "$BAM" \
            -o "${OUTPUT}/aldy/${GENE}.aldy" \
        >> "$log" 2>&1; then
        log_status "DONE   Aldy"
    else
        log_status "FAILED Aldy  (see ${log})"
    fi
    # Aldy failure is non-fatal; always return 0 so it does not kill a background job pool
    return 0
}

run_stellarpgx() {
    local log="${OUTPUT}/logs/stellarpgx.log"
    local bam_dir bam_base
    bam_dir="$(dirname "$BAM")"
    bam_base="$(basename "${BAM%.*}")"   # strips .bam or .cram
    # StellarPGx main.nf uses 'cypor' not 'por' as the gene identifier.
    local stellar_gene="${GENE_LOWER}"
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
        log_status "FAILED StellarPGx  (see ${log})"
    fi
    return 0
}

run_hla() {
    local log="${OUTPUT}/logs/hla.log"
    log_status "START  OptiType HLA typing  (gene: ${GENE})"
    if pgx-hla.sh "$GENE" "$BAM" "$SAMPLE" "$OUTPUT" \
        >> "$log" 2>&1; then
        log_status "DONE   OptiType HLA typing"
    else
        log_status "FAILED OptiType HLA typing  (see ${log})"
    fi
    return 0
}

run_mt() {
    local log="${OUTPUT}/logs/mutserve.log"
    log_status "START  mutserve MT-RNR1 calling"
    if pgx-mt.sh "$BAM" "$SAMPLE" "$OUTPUT" --ref "$REF" \
        >> "$log" 2>&1; then
        log_status "DONE   mutserve MT-RNR1 calling"
    else
        log_status "FAILED mutserve MT-RNR1 calling  (see ${log})"
    fi
    return 0
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
    if pypgx run-ngs-pipeline "${pypgx_args[@]}" >> "$log" 2>&1; then
        log_status "DONE   PyPGx"
    else
        log_status "FAILED PyPGx  (see ${log})"
    fi
    return 0
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
    if stargazer "${stargazer_args[@]}" >> "$log" 2>&1; then
        log_status "DONE   Stargazer"
    else
        log_status "FAILED Stargazer  (see ${log})"
    fi
    return 0
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

[[ -n "$PYPGX_PID"      ]] && wait "$PYPGX_PID"
[[ -n "$STARGAZER_PID"  ]] && wait "$STARGAZER_PID"
[[ -n "$ALDY_PID"       ]] && wait "$ALDY_PID"
[[ -n "$STELLARPGX_PID" ]] && wait "$STELLARPGX_PID"
[[ -n "$HLA_PID"        ]] && wait "$HLA_PID"
[[ -n "$MT_PID"         ]] && wait "$MT_PID"

echo ""
log_status "All tools finished."
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# PHASE 3 — Comparison report
# ─────────────────────────────────────────────────────────────────────────────
echo "------------------------------------------------------------"
log_status "Phase 3: generating comparison report"
echo "------------------------------------------------------------"
python3 /opt/pgx/pgx-compare.py \
    --gene "$GENE" \
    --sample "$SAMPLE" \
    --output-dir "$OUTPUT"
