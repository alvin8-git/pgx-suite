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
#   apptainer pull/build the biocontainer, then copy a built IMGT/HLA reference into it
#   (the biocontainer ships NO usable DB; `arcasHLA reference --version X` needs a git
#   commit hash + GitHub LFS, which fails). Exact steps: docs/ToolsDocumentation.md
#   "Building arcashla.sif". The SIF is a required asset: StellarPGx/containers/arcashla.sif.

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
    log_status "ERROR arcasHLA SIF not found at ${SIF} -- HLA class II (DQA1/DQB1/DRB1) NOT typed; build it (see docs/ToolsDocumentation.md / README 'Building the arcasHLA SIF')"
    exit 1
fi

# ── 1. MHC + chr6-ALT + HLA-contig reads -> paired FASTQ ─────────────────────
# Same read set as the validated host-side route (docs/ToolsDocumentation.md §5b):
# the MHC window plus every chr6 *_alt and HLA-* contig in the BAM header. The raw
# reads go straight to `arcasHLA genotype` (NOT `arcasHLA extract`: extract keeps too
# few reads at WGS depth and silently drops DQA1/DQB1, leaving only DRB1).
MHC="chr6:28510020-33480577"
log_status "Extracting MHC + alt/HLA-contig reads (${MHC})"
{ printf '%s\n' "$MHC"
  samtools view -H "$BAM" | awk -F'\t' '$1=="@SQ"{sub(/^SN:/,"",$2); if ($2 ~ /^chr6_.*_alt$/ || $2 ~ /^HLA-/) print $2}'
} > "${TMP}/regions.txt"
# shellcheck disable=SC2046
samtools view -b -@ "$THREADS" "$BAM" $(cat "${TMP}/regions.txt") 2>/dev/null \
    | samtools sort -n -@ "$THREADS" -o "${TMP}/mhc.nsort.bam" - 2>/dev/null
samtools fastq -@ "$THREADS" -n -1 "${TMP}/hla.1.fq.gz" -2 "${TMP}/hla.2.fq.gz" \
    -s "${TMP}/hla.single.fq.gz" -0 /dev/null "${TMP}/mhc.nsort.bam" 2>/dev/null

APPTAINER_EXEC=(apptainer exec --cleanenv --bind "${TMP}:${TMP}" --bind "${OUTDIR}:${OUTDIR}" "$SIF")

FASTQS="${TMP}/hla.1.fq.gz ${TMP}/hla.2.fq.gz"
if [[ ! -s "${TMP}/hla.1.fq.gz" || ! -s "${TMP}/hla.2.fq.gz" ]]; then
    log_status "ERROR no paired MHC reads extracted from ${BAM}"; exit 1
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
