#!/usr/bin/env bash
# docker/pgx-run.sh — Phase 2 orchestration: run all supported PGx callers for a single gene
#
# Usage: pgx-run.sh <GENE> <BAM_FILE> [--ref /pgx/ref/hg38.fa] [--output /pgx/results]
# Example: pgx-run.sh CYP2D6 /pgx/data/NA12878.bam
#
# SV handling notes:
#   PyPGx:      SV genes (CYP2A6/2B6/2D6/2E1/4F2/G6PD/GSTM1/GSTT1) need an extra
#               prepare-depth-of-coverage + compute-control-statistics (VDR) step.
#               These intermediate files feed --depth-of-coverage and
#               --control-statistics flags in run-ngs-pipeline.
#   Stargazer:  Paralog genes (CYP2A6/2B6/2D6) need a GDF depth-profile
#               generated from the BAM (-G flag). Without the GDF, Stargazer
#               runs in VCF-only mode and cannot detect gene duplications/deletions.
#   Aldy:       CN/SV detection is automatic — handled internally by the ILP solver.
#               No special parameters required.
#   StellarPGx: SV detection is automatic — graphtyper calls SVs natively.
#               No special parameters required.
set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
REF="/pgx/ref/hg38.fa"
OUTPUT="/pgx/results"

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF
Usage: pgx-run.sh <GENE> <BAM_FILE> [options]

Options:
  --ref PATH       GRCh38 reference FASTA (default: /pgx/ref/hg38.fa)
  --output PATH    Output directory (default: /pgx/results)
  -h, --help       Show this help

Examples:
  pgx-run.sh CYP2D6  /pgx/data/sample.bam
  pgx-run.sh CYP2C19 /pgx/data/sample.bam --output /pgx/results
EOF
}

# ── Argument parsing ──────────────────────────────────────────────────────────
if [[ $# -lt 2 ]]; then usage; exit 1; fi

GENE="${1^^}"   # force uppercase
BAM="$2"
shift 2

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)     REF="$2";    shift 2 ;;
        --output)  OUTPUT="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

# ── GRCh38 gene coordinate table ──────────────────────────────────────────────
# Sourced from pypgx gene-table.csv (authoritative regions including upstream/downstream).
# GSTT1 is on an alt contig (chr22_KI270879v1_alt) — bcftools mpileup is skipped for it.
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
)

# ── Gene support matrix (pypgx stargazer aldy stellarpgx) ────────────────────
# 1=supported, 0=not supported by that tool
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
    [NUDT15]="1 0 1 1"
    [TPMT]="1 0 1 1"
    [UGT1A1]="1 0 1 1"
    [SLCO1B1]="1 1 1 1"
    [DPYD]="1 0 1 0"
    [NAT1]="1 0 0 1"
    [NAT2]="1 0 0 1"
    [G6PD]="1 1 0 0"
    [GSTM1]="1 0 0 1"
    [GSTT1]="0 0 0 1"
    [POR]="0 0 0 1"
    [CYPOR]="0 0 0 1"
    [VKORC1]="1 1 0 0"
    [CYP1A1]="0 1 1 1"
    [CYP1A2]="0 1 1 1"
    [CYP2A6]="0 1 1 1"
    [CYP2E1]="0 1 0 1"
    [IFNL3]="1 0 0 0"
    [RYR1]="1 0 0 0"
)

# ── PyPGx SV genes — need depth-of-coverage + control-statistics (VDR) ───────
# Source: pypgx gene-table.csv SV=TRUE column, filtered to our supported set.
# Control gene is always VDR (the only one demonstrated in pypgx tutorials).
PYPGX_SV_GENES=(CYP2A6 CYP2B6 CYP2D6 CYP2E1 CYP4F2 G6PD GSTM1 GSTT1)

# ── Stargazer SV genes — need GDF depth profile created from BAM ──────────────
# Only the three paralog genes benefit from GDF-based CN normalization in Stargazer.
# Source: Stargazer gene_table.tsv, genes where paralog != '.'
# Control gene: vdr (Stargazer accepts: egfr, ryr1, vdr)
STARGAZER_SV_GENES=(CYP2A6 CYP2B6 CYP2D6)

# ── Stargazer control gene mapping (for VCF-only mode, used for non-SV genes) ─
declare -A STARGAZER_CONTROL
STARGAZER_CONTROL=(
    [CYP2D6]="vdr"
    [CYP2C19]="vdr"
    [CYP2C9]="vdr"
    [CYP2B6]="vdr"
    [CYP2C8]="vdr"
    [CYP3A4]="vdr"
    [CYP3A5]="vdr"
    [CYP4F2]="vdr"
    [CYP1A1]="vdr"
    [CYP1A2]="vdr"
    [CYP2A6]="vdr"
    [CYP2E1]="vdr"
    [SLCO1B1]="vdr"
    [G6PD]="vdr"
    [VKORC1]="vdr"
)

# ── Helper: check if value is in array ───────────────────────────────────────
in_array() {
    local needle="$1"; shift
    for item in "$@"; do [[ "$item" == "$needle" ]] && return 0; done
    return 1
}

GENE_LOWER="${GENE,,}"

echo "============================================================"
echo " PGx Suite — Orchestration Runner"
echo " Gene:   ${GENE}"
echo " BAM:    ${BAM}"
echo " Ref:    ${REF}"
echo " Output: ${OUTPUT}"
echo " Build:  GRCh38"
echo "============================================================"
echo ""

# ── Input validation ──────────────────────────────────────────────────────────
echo "[Validate] Checking inputs..."

if [[ -z "${GENE_SUPPORT[$GENE]:-}" ]]; then
    echo "ERROR: Gene '${GENE}' is not in the supported gene list." >&2
    echo "Supported: ${!GENE_SUPPORT[*]}" >&2
    exit 1
fi

if [[ ! -f "$BAM" ]]; then
    echo "ERROR: BAM file not found: $BAM" >&2; exit 1
fi

BAI="${BAM}.bai"
[[ ! -f "$BAI" ]] && BAI="${BAM%.bam}.bai"
if [[ ! -f "$BAI" ]]; then
    echo "ERROR: BAM index not found. Run: samtools index $BAM" >&2; exit 1
fi

if [[ ! -f "$REF" ]]; then
    echo "ERROR: Reference FASTA not found: $REF" >&2
    echo "       Mount GRCh38 FASTA with: -v /path/to/ref:/pgx/ref" >&2; exit 1
fi
if [[ ! -f "${REF}.fai" ]]; then
    echo "ERROR: Reference index not found: ${REF}.fai — run: samtools faidx $REF" >&2; exit 1
fi

echo "  Gene:      ${GENE} — OK"
echo "  BAM:       ${BAM} — OK"
echo "  BAM index: ${BAI} — OK"
echo "  Reference: ${REF} — OK"
echo ""

# ── Sample name from BAM read group ──────────────────────────────────────────
SAMPLE=$(samtools view -H "$BAM" | grep '^@RG' | grep -oP 'SM:\K[^\t]+' | head -1 || true)
[[ -z "$SAMPLE" ]] && SAMPLE=$(basename "$BAM" .bam)
echo "[Info] Sample name: ${SAMPLE}"
echo ""

# ── Parse tool support flags ──────────────────────────────────────────────────
read -r DO_PYPGX DO_STARGAZER DO_ALDY DO_STELLARPGX <<< "${GENE_SUPPORT[$GENE]}"

IS_PYPGX_SV=0
in_array "$GENE" "${PYPGX_SV_GENES[@]}" && IS_PYPGX_SV=1

IS_STARGAZER_SV=0
in_array "$GENE" "${STARGAZER_SV_GENES[@]}" && IS_STARGAZER_SV=1

# StellarPGx requires mounted volumes
if [[ "$DO_STELLARPGX" -eq 1 ]] && \
   [[ ! -d "/pgx/stellarpgx" || ! -f "/pgx/containers/stellarpgx-dev.sif" ]]; then
    echo "[Warn] StellarPGx volumes not mounted — skipping"
    echo "       Requires: -v \$(pwd)/StellarPGx:/pgx/stellarpgx"
    echo "                 -v \$(pwd)/StellarPGx/containers:/pgx/containers"
    DO_STELLARPGX=0
fi

echo "[Info] Tool support for ${GENE}:"
echo "  PyPGx:      $([[ $DO_PYPGX -eq 1 ]] && echo YES || echo NO)  (SV mode: $([[ $IS_PYPGX_SV -eq 1 ]] && echo YES || echo NO))"
echo "  Stargazer:  $([[ $DO_STARGAZER -eq 1 ]] && echo YES || echo NO)  (SV/GDF mode: $([[ $IS_STARGAZER_SV -eq 1 ]] && echo YES || echo NO))"
echo "  Aldy:       $([[ $DO_ALDY -eq 1 ]] && echo YES || echo NO)  (SV: auto from BAM)"
echo "  StellarPGx: $([[ $DO_STELLARPGX -eq 1 ]] && echo YES || echo NO)  (SV: auto via graphtyper)"
echo ""

# ── Create output directories ─────────────────────────────────────────────────
mkdir -p \
    "${OUTPUT}/pypgx" \
    "${OUTPUT}/stargazer" \
    "${OUTPUT}/aldy" \
    "${OUTPUT}/stellarpgx"

VCF="${OUTPUT}/${GENE}.vcf.gz"
COORDS="${GENE_COORDS[$GENE]:-}"

# ── Shared preprocessing: bcftools mpileup → gene-region VCF ─────────────────
# Used by PyPGx (variant input) and Stargazer (variant input).
# Skipped for GSTT1 (alt contig) and when neither tool is run.
NEED_VCF=$(( (DO_PYPGX + DO_STARGAZER) > 0 ? 1 : 0 ))

if [[ "$NEED_VCF" -eq 1 ]]; then
    if [[ "$COORDS" == "ALT_CONTIG" ]]; then
        echo "[Step 1] Skipping bcftools — ${GENE} is on an alt contig (${GENE_COORDS[$GENE]})"
        echo "         PyPGx will use depth-of-coverage only; Stargazer requires manual VCF"
        echo ""
        NEED_VCF=0
    elif [[ -z "$COORDS" ]]; then
        echo "ERROR: No GRCh38 coordinates defined for ${GENE}" >&2; exit 1
    else
        echo "[Step 1] Variant calling — region: ${COORDS}"
        bcftools mpileup \
            -r "$COORDS" \
            -f "$REF" \
            -a AD,DP \
            --max-depth 500 \
            -o - \
            "$BAM" \
        | bcftools call \
            -m -v \
            --output-type z \
            -o "$VCF"
        tabix -p vcf "$VCF"
        echo "  VCF written: ${VCF}"
        echo ""
    fi
fi

# ── PyPGx SV pre-processing ───────────────────────────────────────────────────
# For SV genes: compute depth-of-coverage and control statistics (VDR) from BAM.
# These are required inputs for run-ngs-pipeline to call copy number variations.
DEPTH_ZIP="${OUTPUT}/depth-of-coverage.zip"
CTRL_ZIP="${OUTPUT}/control-stats-VDR.zip"

if [[ "$DO_PYPGX" -eq 1 && "$IS_PYPGX_SV" -eq 1 ]]; then
    echo "[Step 2a-pre] PyPGx SV preprocessing for ${GENE}"
    echo "  Computing depth-of-coverage..."
    pypgx prepare-depth-of-coverage \
        "$DEPTH_ZIP" \
        "$BAM" \
        --assembly GRCh38
    echo "  Computing VDR control statistics..."
    pypgx compute-control-statistics \
        VDR \
        "$CTRL_ZIP" \
        "$BAM" \
        --assembly GRCh38
    echo "  PyPGx SV inputs ready."
    echo ""
fi

# ── PyPGx ─────────────────────────────────────────────────────────────────────
if [[ "$DO_PYPGX" -eq 1 ]]; then
    echo "[Step 2a] PyPGx"
    PYPGX_ARGS=(
        "$GENE"
        "${OUTPUT}/pypgx"
        --assembly GRCh38
        --force
    )
    [[ "$NEED_VCF" -eq 1 ]] && PYPGX_ARGS+=(--variants "$VCF")
    if [[ "$IS_PYPGX_SV" -eq 1 ]]; then
        PYPGX_ARGS+=(
            --depth-of-coverage "$DEPTH_ZIP"
            --control-statistics "$CTRL_ZIP"
        )
        echo "  Mode: SV (depth-of-coverage + VDR control statistics)"
    else
        echo "  Mode: standard (VCF-only)"
    fi
    pypgx run-ngs-pipeline "${PYPGX_ARGS[@]}" \
        && echo "  PyPGx: OK" \
        || echo "  PyPGx: FAILED (continuing)"
    echo ""
fi

# ── Stargazer ─────────────────────────────────────────────────────────────────
# For SV genes (CYP2A6/2B6/2D6): create GDF depth profile from BAM first,
# then run with -c vdr -g <gdf_file> to enable CN-based SV detection.
# For all other genes: VCF-only mode (no CN calling).
if [[ "$DO_STARGAZER" -eq 1 ]]; then
    echo "[Step 2b] Stargazer"
    if [[ "$IS_STARGAZER_SV" -eq 1 ]]; then
        echo "  Mode: SV (creating GDF depth profile from BAM, control=VDR)"
        GDF_DIR="${OUTPUT}/stargazer/gdf"
        GDF_FILE="${GDF_DIR}/${GENE_LOWER}.gdf"
        mkdir -p "$GDF_DIR"
        echo "  Step 2b-pre: bam2gdf → ${GDF_FILE}"
        stargazer \
            -G "${GENE_LOWER}.gdf" \
            -t "$GENE_LOWER" \
            -c vdr \
            -B "$BAM" \
            -o "$GDF_DIR" \
            -a grc38 \
            -d wgs \
            && echo "  GDF created: ${GDF_FILE}" \
            || { echo "  GDF creation FAILED — falling back to VCF-only mode"; IS_STARGAZER_SV=0; }

        if [[ "$IS_STARGAZER_SV" -eq 1 && -f "$GDF_FILE" && "$NEED_VCF" -eq 1 ]]; then
            stargazer \
                -t "$GENE_LOWER" \
                -d wgs \
                -a grc38 \
                -i "$VCF" \
                -c vdr \
                -g "$GDF_FILE" \
                -o "${OUTPUT}/stargazer" \
                && echo "  Stargazer: OK (SV mode)" \
                || echo "  Stargazer: FAILED (continuing)"
        elif [[ "$NEED_VCF" -eq 1 ]]; then
            echo "  Falling back to VCF-only mode"
            CONTROL="${STARGAZER_CONTROL[$GENE]:-}"
            STARGAZER_EXTRA=()
            [[ -n "$CONTROL" ]] && STARGAZER_EXTRA+=(-c "$CONTROL")
            stargazer \
                -t "$GENE_LOWER" \
                -d wgs \
                -a grc38 \
                -i "$VCF" \
                -o "${OUTPUT}/stargazer" \
                "${STARGAZER_EXTRA[@]+"${STARGAZER_EXTRA[@]}"}" \
                && echo "  Stargazer: OK (VCF-only)" \
                || echo "  Stargazer: FAILED (continuing)"
        else
            echo "  Stargazer: SKIPPED — no VCF available"
        fi
    else
        echo "  Mode: VCF-only (gene has no paralog; no GDF needed)"
        if [[ "$NEED_VCF" -eq 1 ]]; then
            CONTROL="${STARGAZER_CONTROL[$GENE]:-}"
            STARGAZER_EXTRA=()
            [[ -n "$CONTROL" ]] && STARGAZER_EXTRA+=(-c "$CONTROL")
            stargazer \
                -t "$GENE_LOWER" \
                -d wgs \
                -a grc38 \
                -i "$VCF" \
                -o "${OUTPUT}/stargazer" \
                "${STARGAZER_EXTRA[@]+"${STARGAZER_EXTRA[@]}"}" \
                && echo "  Stargazer: OK" \
                || echo "  Stargazer: FAILED (continuing)"
        else
            echo "  Stargazer: SKIPPED — no VCF available (alt contig gene)"
        fi
    fi
    echo ""
fi

# ── Aldy ──────────────────────────────────────────────────────────────────────
# Aldy reads the BAM directly and handles CN/SV detection automatically via ILP.
# No special parameters needed for SV genes.
if [[ "$DO_ALDY" -eq 1 ]]; then
    echo "[Step 2c] Aldy (SV/CN detection: automatic from BAM)"
    aldy genotype \
        -g "$GENE" \
        -p hg38 \
        "$BAM" \
        -o "${OUTPUT}/aldy/${GENE}.aldy" \
        && echo "  Aldy: OK" \
        || echo "  Aldy: FAILED (continuing)"
    echo ""
fi

# ── StellarPGx ────────────────────────────────────────────────────────────────
# StellarPGx uses graphtyper for variant calling, which handles SVs natively.
# No special parameters needed for SV genes.
if [[ "$DO_STELLARPGX" -eq 1 ]]; then
    echo "[Step 2d] StellarPGx (SV detection: automatic via graphtyper)"
    BAM_DIR="$(dirname "$BAM")"
    BAM_BASE="$(basename "$BAM" .bam)"
    nextflow run /pgx/stellarpgx/main.nf \
        --gene "${GENE_LOWER}" \
        --in_bam "${BAM_DIR}/${BAM_BASE}*{bam,bai}" \
        --ref_file "$REF" \
        --out_dir "${OUTPUT}/stellarpgx" \
        -work-dir "${OUTPUT}/stellarpgx/.work" \
        && echo "  StellarPGx: OK" \
        || echo "  StellarPGx: FAILED (continuing)"
    echo ""
fi

# ── Comparison table ──────────────────────────────────────────────────────────
echo "[Step 3] Generating comparison report..."
python3 /opt/pgx/pgx-compare.py \
    --gene "$GENE" \
    --sample "$SAMPLE" \
    --output-dir "$OUTPUT"
