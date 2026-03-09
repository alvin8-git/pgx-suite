#!/usr/bin/env bash
# docker/pgx-hla.sh — HLA-A/HLA-B typing using OptiType (via Apptainer)
#
# Usage: pgx-hla.sh <GENE> <BAM> <SAMPLE> <OUTPUT_DIR>
#   GENE:       HLA-A or HLA-B
#   BAM:        Input BAM/CRAM file (GRCh38)
#   SAMPLE:     Sample name
#   OUTPUT_DIR: Gene output directory (e.g. /pgx/results/HLA-A)
#
# Output:
#   <OUTPUT_DIR>/optitype/<SAMPLE>_result.tsv
#   <OUTPUT_DIR>/optitype/<SAMPLE>_coverage_plot.pdf
#
# Method:
#   1. samtools view extracts reads from the full MHC region (chr6:28,510,020-33,480,577)
#      covering HLA-A (chr6:29,910,247-29,913,661) and HLA-B (chr6:31,321,649-31,324,666).
#   2. samtools fastq converts extracted reads to FASTQ (R1+R2 for paired, single for SE).
#   3. OptiType (Apptainer container) types HLA Class I (A, B, C) via ILP optimisation.
#      Internally OptiType runs razers3 against its bundled IMGT/HLA reference.
#   4. pgx-compare.py reads the result TSV and extracts A1/A2 (for HLA-A) or B1/B2
#      (for HLA-B) to produce the per-gene comparison TSV and detail JSON.
#
# Container requirement:
#   The OptiType Apptainer SIF must be present at /pgx/containers/optitype.sif
#   Pull once with:
#     apptainer pull --name optitype.sif \
#       docker://quay.io/biocontainers/optitype:1.3.5--hdfd78af_1
#
# Timing profile (30X WGS, warm disk cache):
#   samtools view + sort + fastq (MHC region, ~50 MB)   ~30-60 sec
#   OptiType ILP (razers3 + optimisation)               ~2-5 min
#   Total                                               ~3-6 min
#
# Notes:
#   - Both HLA-A and HLA-B pgx-run.sh calls run full OptiType independently;
#     each writes results to its own OUTPUT_DIR/optitype/. This is slightly
#     redundant but keeps each gene call self-contained.
#   - MHC extraction uses the primary assembly region. HLA-A and HLA-B genes
#     are both within chr6:28.5 Mb-33.5 Mb (GRCh38).

set -uo pipefail

if [[ $# -lt 4 ]]; then
    echo "Usage: pgx-hla.sh <GENE> <BAM> <SAMPLE> <OUTPUT_DIR>" >&2
    exit 1
fi

GENE="$1"
BAM="$2"
SAMPLE="$3"
OUTPUT_DIR="$4"

# Full MHC primary-assembly region — covers HLA-A, HLA-B, HLA-C
MHC_REGION="chr6:28510020-33480577"
NCPU=$(nproc 2>/dev/null || echo 4)
OPTITYPE_SIF="/pgx/containers/optitype.sif"
OPTITYPE_OUT="${OUTPUT_DIR}/optitype"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

mkdir -p "$OPTITYPE_OUT"

# ── Timing helper ────────────────────────────────────────────────────────────
_t0()      { date +%s%3N; }
_elapsed() { echo "$(( $(date +%s%3N) - $1 )) ms"; }

# ── Dependency checks ────────────────────────────────────────────────────────
for _cmd in samtools apptainer; do
    if ! command -v "$_cmd" &>/dev/null; then
        echo "ERROR: $_cmd not found in PATH" >&2; exit 1
    fi
done

if [[ ! -f "$OPTITYPE_SIF" ]]; then
    echo "ERROR: OptiType SIF not found at $OPTITYPE_SIF" >&2
    echo "Pull with:" >&2
    echo "  apptainer pull --name optitype.sif \\" >&2
    echo "    docker://quay.io/biocontainers/optitype:1.3.5--hdfd78af_1" >&2
    exit 1
fi

echo "[hla] Gene:   ${GENE}"
echo "[hla] BAM:    ${BAM}"
echo "[hla] Region: ${MHC_REGION}"

# ── Step 1: Extract MHC reads ─────────────────────────────────────────────────
echo "[hla] Extracting reads from MHC region …"
_T=$(_t0)

FASTQ_1="${TMP}/mhc_1.fastq"
FASTQ_2="${TMP}/mhc_2.fastq"
FASTQ_S="${TMP}/mhc_s.fastq"

# Detect paired-end by checking flag 0x1 (paired read) in MHC region
paired_count=$(samtools view -f 0x1 -c -@ "$NCPU" "$BAM" "$MHC_REGION" 2>/dev/null || echo 0)

if [[ "$paired_count" -gt 0 ]]; then
    # Paired-end: extract, name-sort, then split to R1/R2
    samtools view -b -@ "$NCPU" "$BAM" "$MHC_REGION" 2>/dev/null \
        | samtools sort -n -@ "$NCPU" \
        | samtools fastq -@ "$NCPU" \
            -1 "$FASTQ_1" \
            -2 "$FASTQ_2" \
            -s "$FASTQ_S" \
            -0 /dev/null \
            -n 2>/dev/null
    echo "[hla] MHC extraction (paired-end): $(_elapsed $_T)"
    OPTITYPE_INPUT="-i $FASTQ_1 $FASTQ_2"
else
    # Single-end: stream directly to FASTQ
    samtools view -b -@ "$NCPU" "$BAM" "$MHC_REGION" 2>/dev/null \
        | samtools fastq -@ "$NCPU" - > "$FASTQ_S" 2>/dev/null
    echo "[hla] MHC extraction (single-end): $(_elapsed $_T)"
    OPTITYPE_INPUT="-i $FASTQ_S"
fi

n_reads=$(wc -l < "${FASTQ_1:-$FASTQ_S}" 2>/dev/null || echo 0)
n_reads=$(( n_reads / 4 ))
echo "[hla] MHC reads extracted: ${n_reads}"

if [[ "$n_reads" -eq 0 ]]; then
    echo "[hla] WARN: No reads extracted from ${MHC_REGION} — BAM may not have chr-prefix contigs" >&2
    # Try without chr prefix
    MHC_REGION_NOCHR="6:28510020-33480577"
    samtools view -b -@ "$NCPU" "$BAM" "$MHC_REGION_NOCHR" 2>/dev/null \
        | samtools sort -n -@ "$NCPU" \
        | samtools fastq -@ "$NCPU" \
            -1 "$FASTQ_1" -2 "$FASTQ_2" -s "$FASTQ_S" -0 /dev/null -n 2>/dev/null \
        || true
    n_reads=$(wc -l < "${FASTQ_1:-$FASTQ_S}" 2>/dev/null || echo 0)
    n_reads=$(( n_reads / 4 ))
    echo "[hla] Retry (no-chr prefix): ${n_reads} reads"
fi

# ── Step 2: Run OptiType via Apptainer ───────────────────────────────────────
echo "[hla] Running OptiType …"
_T=$(_t0)

apptainer exec \
    --bind "${TMP}:${TMP}" \
    --bind "${OPTITYPE_OUT}:${OPTITYPE_OUT}" \
    "$OPTITYPE_SIF" \
    OptiTypePipeline.py \
        $OPTITYPE_INPUT \
        --dna \
        -o "$OPTITYPE_OUT" \
        -p "$SAMPLE" \
        -v

echo "[hla] OptiType: $(_elapsed $_T)"

# ── Step 3: Verify output ─────────────────────────────────────────────────────
RESULT_TSV="${OPTITYPE_OUT}/${SAMPLE}_result.tsv"
if [[ -f "$RESULT_TSV" ]]; then
    echo "[hla] HLA typing complete. Result: ${RESULT_TSV}"
    cat "$RESULT_TSV"
else
    echo "[hla] WARN: Expected result TSV not found: ${RESULT_TSV}" >&2
    # Check for timestamped filename (older OptiType versions)
    RESULT_TSV_GLOB=$(ls "${OPTITYPE_OUT}"/*_result.tsv 2>/dev/null | tail -1 || true)
    if [[ -n "$RESULT_TSV_GLOB" ]]; then
        echo "[hla] Found timestamped result: ${RESULT_TSV_GLOB}"
        # Symlink to the expected filename for pgx-compare.py
        ln -sf "$(basename "$RESULT_TSV_GLOB")" "$RESULT_TSV"
    else
        echo "[hla] ERROR: No OptiType result TSV found in ${OPTITYPE_OUT}" >&2
        exit 1
    fi
fi
