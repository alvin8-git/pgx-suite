#!/usr/bin/env bash
# pgx-hla2.sh — HLA class-II typing (DQA1/DQB1/DRB1 + class I) via arcasHLA
#
# Usage: pgx-hla2.sh <BAM|CRAM> <SAMPLE> <HLA2_OUTDIR> [--ref REF] [--sif SIF] [--threads N]
#
# Runs ONCE per sample (like PharmCAT). OptiType (pgx-hla.sh) types class I only;
# arcasHLA adds class II, enabling celiac (DQ2.5/DQ8), narcolepsy (DQB1*06:02),
# and T1D susceptibility reporting. Per-gene compare (parse_arcashla in
# pgx-compare.py) reads the JSON below and emits the SAME *_comparison.tsv schema
# as OptiType, so OmniGen needs no new reader.
#
# Output: <HLA2_OUTDIR>/<SAMPLE>.genotype.json   (arcasHLA genotype)
#
# Requires: arcasHLA Apptainer SIF at /pgx/containers/arcashla.sif, built with its
# IMGT/HLA reference DB. Pull once:
#   apptainer pull --name arcashla.sif docker://quay.io/biocontainers/arcas-hla:0.6.0--hdfd78af_1
#   apptainer exec arcashla.sif arcasHLA reference --version 3.55.0   # build DB (one-off)

set -euo pipefail

SCRIPT=$(basename "$0")
log_status() { echo "[$(date '+%H:%M:%S')] ${SCRIPT}: $*" >&2; }

REF="/pgx/ref/hg38.fa"
SIF="/pgx/containers/arcashla.sif"
THREADS=4
# HLA genes to type: class I (cross-check vs OptiType) + class II (the new value).
GENES="A,B,C,DQA1,DQB1,DRB1"

BAM="${1:-}"; SAMPLE="${2:-}"; OUTDIR="${3:-}"; shift 3 || true
while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)     REF="$2";     shift 2 ;;
        --sif)     SIF="$2";     shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --genes)   GENES="$2";   shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done
[[ -n "$BAM" && -n "$SAMPLE" && -n "$OUTDIR" ]] || {
    echo "Usage: pgx-hla2.sh <BAM> <SAMPLE> <HLA2_OUTDIR> [--ref REF] [--sif SIF]" >&2; exit 1; }

mkdir -p "$OUTDIR"
TMP=$(mktemp -d "${OUTDIR}/tmp_XXXXXX")
trap 'rm -rf "$TMP"' EXIT

if [[ ! -f "$SIF" ]]; then
    log_status "ERROR arcasHLA SIF not found at ${SIF} — see header for pull command"
    exit 1
fi

# ── 1. Extract MHC reads (chr6 MHC; arcasHLA re-partitions internally) ────────
# arcasHLA extract works from a name-sorted / paired BAM of MHC reads.
MHC="chr6:28510020-33480577"
log_status "Extracting MHC reads (${MHC})"
samtools view -b -@ "$THREADS" "$BAM" "$MHC" 2>/dev/null \
    | samtools sort -n -@ "$THREADS" -o "${TMP}/mhc.bam" - 2>/dev/null
samtools index -@ "$THREADS" "${TMP}/mhc.bam" 2>/dev/null || true

APPTAINER_EXEC=(apptainer exec --cleanenv --bind "${TMP}:${TMP}" --bind "${OUTDIR}:${OUTDIR}" "$SIF")

# ── 2. arcasHLA extract (BAM -> FASTQ of HLA reads) ──────────────────────────
log_status "arcasHLA extract"
"${APPTAINER_EXEC[@]}" arcasHLA extract "${TMP}/mhc.bam" \
    -o "${TMP}/extract" -t "$THREADS" --temp "${TMP}/at" \
    > "${OUTDIR}/arcashla_extract.log" 2>&1 || log_status "WARN arcasHLA extract non-zero"

FASTQS=$(ls "${TMP}/extract"/*.fq.gz 2>/dev/null | tr '\n' ' ' || true)
if [[ -z "$FASTQS" ]]; then
    log_status "ERROR arcasHLA extract produced no FASTQ"; exit 1
fi

# ── 3. arcasHLA genotype (class I + II) ──────────────────────────────────────
log_status "arcasHLA genotype (${GENES})"
"${APPTAINER_EXEC[@]}" arcasHLA genotype ${FASTQS} \
    -g "$GENES" -o "${TMP}/geno" -t "$THREADS" --temp "${TMP}/gt" \
    > "${OUTDIR}/arcashla_genotype.log" 2>&1 || log_status "WARN arcasHLA genotype non-zero"

# arcasHLA writes <prefix>.genotype.json; normalise its filename to <SAMPLE>.genotype.json.
SRC=$(ls "${TMP}/geno"/*.genotype.json 2>/dev/null | head -1 || true)
if [[ -z "$SRC" ]]; then
    log_status "ERROR no arcasHLA genotype.json produced"; exit 1
fi
cp "$SRC" "${OUTDIR}/${SAMPLE}.genotype.json"
log_status "Wrote ${OUTDIR}/${SAMPLE}.genotype.json"
