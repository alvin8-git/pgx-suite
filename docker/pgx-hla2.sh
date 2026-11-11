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
# Output: <HLA2_OUTDIR>/<SAMPLE>.genotype.json          (arcasHLA genotype; unchanged contract)
#         <HLA2_OUTDIR>/arcashla_genotype.log           (arcasHLA --verbose: read counts,
#                                                        abundance, alternative pairs)
#         <HLA2_OUTDIR>/<SAMPLE>.arcashla_genes.json    (reads aligned per gene)
#         <HLA2_OUTDIR>/<SAMPLE>.t1k_genotype.tsv       (T1K second typer, DRB1/DQA1/DQB1 + all)
#         <HLA2_OUTDIR>/<SAMPLE>.hla2_consensus.json    (two-typer consensus + agreement flag)
# The second typer (T1K, MIT, same IPD-IMGT/HLA release) needs `run-t1k` on PATH (baked
# into the pgx-suite image) and the index dir (--t1k-idx, default
# /pgx/containers/t1k_hlaidx, built once from hla.dat with t1k-build.pl -- see
# docs/ToolsDocumentation.md "HLA class II second typer"). If either is missing the run
# still completes with arcasHLA only and the consensus says single_typer.
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
T1K_IDX="/pgx/containers/t1k_hlaidx"
THREADS=4
# HLA genes to type: class I (cross-check vs OptiType) + class II (the new value).
GENES="A,B,C,DQA1,DQB1,DRB1"

BAM="${1:-}"; SAMPLE="${2:-}"; OUTDIR="${3:-}"; shift 3 || true
while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)     REF="$2";     shift 2 ;;
        --sif)     SIF="$2";     shift 2 ;;
        --t1k-idx) T1K_IDX="$2"; shift 2 ;;
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
# samtools writes a 28-byte BGZF EOF block even for zero reads, so -s is never false;
# test for decompressed content instead (|| true: head closing the pipe sends zcat SIGPIPE).
if [[ -z "$(zcat "${TMP}/hla.1.fq.gz" | head -c1 || true)" \
   || -z "$(zcat "${TMP}/hla.2.fq.gz" | head -c1 || true)" ]]; then
    log_status "ERROR no paired MHC reads extracted from ${BAM}"; exit 1
fi

# ── 3. arcasHLA genotype (class I + II) ──────────────────────────────────────
log_status "arcasHLA genotype (${GENES})"
"${APPTAINER_EXEC[@]}" arcasHLA genotype ${FASTQS} \
    -g "$GENES" -o "${TMP}/geno" -t "$THREADS" --temp "${TMP}/gt" --verbose \
    > "${OUTDIR}/arcashla_genotype.log" 2>&1 || log_status "WARN arcasHLA genotype non-zero"

# arcasHLA writes <prefix>.genotype.json; normalise its filename to <SAMPLE>.genotype.json.
SRC=$(ls "${TMP}/geno"/*.genotype.json 2>/dev/null | head -1 || true)
if [[ -z "$SRC" ]]; then
    log_status "ERROR no arcasHLA genotype.json produced"; exit 1
fi
cp "$SRC" "${OUTDIR}/${SAMPLE}.genotype.json"
log_status "Wrote ${OUTDIR}/${SAMPLE}.genotype.json"
GJ=$(ls "${TMP}/geno"/*.genes.json 2>/dev/null | head -1 || true)
[[ -n "$GJ" ]] && cp "$GJ" "${OUTDIR}/${SAMPLE}.arcashla_genes.json"

# ── 4. Second typer: T1K (independent algorithm, same IMGT/HLA release) ──────
T1K_ARGS=()
if command -v run-t1k >/dev/null 2>&1 && [[ -f "${T1K_IDX}/hlaidx_dna_seq.fa" ]]; then
    log_status "T1K genotype (second typer, --preset hla-wgs)"
    mkdir -p "${TMP}/t1k"
    if run-t1k -1 "${TMP}/hla.1.fq.gz" -2 "${TMP}/hla.2.fq.gz" --preset hla-wgs \
           -f "${T1K_IDX}/hlaidx_dna_seq.fa" -t "$THREADS" -o t1k --od "${TMP}/t1k" \
           > "${OUTDIR}/t1k_genotype.log" 2>&1 && [[ -s "${TMP}/t1k/t1k_genotype.tsv" ]]; then
        cp "${TMP}/t1k/t1k_genotype.tsv" "${OUTDIR}/${SAMPLE}.t1k_genotype.tsv"
        T1K_ARGS=(--t1k "${OUTDIR}/${SAMPLE}.t1k_genotype.tsv")
    else
        log_status "WARN T1K failed (see t1k_genotype.log); consensus will be single-typer"
    fi
else
    log_status "WARN T1K or its index (${T1K_IDX}) not available; consensus will be single-typer"
fi

# ── 5. Consensus + agreement flag (never fatal; genotype.json is already written) ──
HC="$(dirname "$(readlink -f "$0")")/hla2_consensus.py"
python3 "$HC" --sample "$SAMPLE" --arcas-json "${OUTDIR}/${SAMPLE}.genotype.json" \
    --arcas-verbose "${OUTDIR}/arcashla_genotype.log" ${T1K_ARGS[@]+"${T1K_ARGS[@]}"} \
    --out "${OUTDIR}/${SAMPLE}.hla2_consensus.json" \
    || log_status "WARN hla2_consensus.py failed; hla2_consensus.json not written"
