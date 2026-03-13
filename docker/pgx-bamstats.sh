#!/usr/bin/env bash
# docker/pgx-bamstats.sh — Compute BAM QC statistics for PGx reporting
#
# Tool usage:
#   samtools idxstats  — total + mapped reads from .bai index (no BAM scan)
#   samtools stats     — read length, insert size, error rate (chr1 4 Mb window)
#   mosdepth           — per-gene mean depth + ≥20×/≥30× coverage fractions
#   samtools view -c   — MAPQ≥20 + duplicate fraction (chr1 4 Mb window)
#
# Note on flagstat alternatives tested:
#   samtools flagstat (~580 sec) and sambamba flagstat (~563 sec) on a 167 GB WGS
#   BAM are both I/O bound at the disk's 490 MB/s ceiling — no tool parallelism
#   helps. idxstats reads from the .bai index in 28 ms and is the right choice.
#   Duplicate% is estimated from the chr1 4 Mb window (accurate to ±1–2%).
#   sambamba remains installed for optional use outside this script.
#
# Dependencies: mosdepth ≥0.3.12, samtools ≥1.10
#
# Usage (called from pgx-all-genes.sh, typically run in background):
#   pgx-bamstats.sh <BAM|CRAM> <SAMPLE> <OUTPUT_DIR> [<REF_FASTA>]
#
# Output:
#   <OUTPUT_DIR>/bam_stats.json
#
# Timing profile (167 GB WGS ~30× depth, observed, cold disk cache):
#   samtools idxstats            ~30–50 ms   (reads .bai index only)
#   samtools stats (chr1 4 Mb)  ~2 sec
#   mosdepth (14 regions)       ~550–570 sec (I/O bound at 490 MB/s HDD)
#   samtools view -c ×3         ~400–550 ms
#   Total                       ~9–10 min
#
# mosdepth dominates (>99% of runtime) due to random BAM access across 14 gene
# regions in a 167 GB file. No tool or threading choice changes this — it is a
# fundamental disk bandwidth limit. The script is designed to run in the
# background (pgx-all-genes.sh launches it with &) so it overlaps with gene
# calling (~10 min) and does not add to total pipeline wall time.

set -uo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: pgx-bamstats.sh <BAM> <SAMPLE> <OUTPUT_DIR> [<REF_FASTA>]" >&2
    exit 1
fi

BAM="$1"
SAMPLE="$2"
OUT_DIR="$3"
REF="${4:-}"

JSON_OUT="${OUT_DIR}/bam_stats.json"
mkdir -p "$OUT_DIR"

NCPU=$(nproc 2>/dev/null || echo 4)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# ── Timing helper ──────────────────────────────────────────────────────────────
_t0()      { date +%s%3N; }
_elapsed() { local t0=$1; echo "$(( $(date +%s%3N) - t0 )) ms"; }

# ── Validate dependencies ──────────────────────────────────────────────────────
for _cmd in mosdepth samtools; do
    if ! command -v "$_cmd" &>/dev/null; then
        echo "ERROR: $_cmd not found in PATH." >&2; exit 1
    fi
done

# ── PGx gene regions BED (0-based, half-open, 4 columns) ─────────────────────
# Primary GRCh38 loci only. GSTT1 (alt contig) and genes without a defined
# primary locus are omitted — their depth is not clinically actionable here.
PGX_BED="${TMP}/pgx_regions.bed"
# 27 primary-assembly PGx genes (0-based half-open BED, GRCh38).
# Coordinates from pgx-run.sh GENE_COORDS (1-based) with start−1.
# GSTT1 is on chr22_KI270879v1_alt (alt contig); omitted from BED —
# its depth entry in bam_stats.json is set to null with note "alt_contig".
cat > "$PGX_BED" << 'BEDEOF'
chr4	88085264	88236626	ABCG2
chr1	201006955	201095000	CACNA1S
chrM	647	1601	MT-RNR1
chr15	74716540	74728528	CYP1A1
chr15	74745843	74759607	CYP1A2
chr19	40833539	40890447	CYP2A6
chr19	40921280	41028398	CYP2B6
chr10	95033770	95072497	CYP2C8
chr10	94935656	94993091	CYP2C9
chr10	94759679	94858547	CYP2C19
chr22	42116497	42155810	CYP2D6
chr10	133517361	133549123	CYP2E1
chr7	99753965	99787184	CYP3A4
chr7	99645192	99682996	CYP3A5
chr19	15863021	15913074	CYP4F2
chr1	97074741	97924034	DPYD
chrX	154528388	154550018	G6PD
chr1	109684815	109696745	GSTM1
chr6	29905000	29920000	HLA-A
chr6	31316000	31330000	HLA-B
chr19	39240551	39253525	IFNL3
chr8	18207107	18226689	NAT1
chr8	18388280	18404218	NAT2
chr13	48034724	48050221	NUDT15
chr7	75912153	75989855	POR
chr19	38430689	38590564	RYR1
chr12	21128192	21242796	SLCO1B1
chr6	18125309	18158169	TPMT
chr2	233754268	233779300	UGT1A1
chr16	31087852	31097797	VKORC1
BEDEOF

STATS_REGION="chr1:10000000-14000000"

# ── CRAM reference flags ────────────────────────────────────────────────────────
# samtools and mosdepth need the reference to decode CRAM data.
# idxstats reads from the .crai index and does NOT need a reference.
#
# NOTE: samtools uses DIFFERENT reference flags per sub-command:
#   samtools stats  → -r FILE   (not -T; -T is silently ignored by stats)
#   samtools view   → -T FILE
#   mosdepth        → --fasta FILE
_EXT="${BAM##*.}"
_ST_STATS_REF_FLAG=()   # for samtools stats
_ST_VIEW_REF_FLAG=()    # for samtools view
_MSD_REF_FLAG=()        # for mosdepth
if [[ "$_EXT" == "cram" ]]; then
    if [[ -n "$REF" && -f "$REF" ]]; then
        _ST_STATS_REF_FLAG=(-r "$REF")
        _ST_VIEW_REF_FLAG=(-T "$REF")
        _MSD_REF_FLAG=(--fasta "$REF")
        echo "[bamstats] CRAM input — will use reference: ${REF}"
    else
        echo "[bamstats] WARN: CRAM input but no reference supplied; samtools/mosdepth may fail." >&2
    fi
fi

# ── 1. samtools idxstats — total + mapped reads from .bai index (28 ms) ───────
# No BAM scan. Columns: chrom  chrom_len  mapped  unmapped.
# Also used for chrX/Y depth estimation (sex inference) and genome size.
echo "[bamstats] Running idxstats …"
_T=$(_t0)
IDXSTATS=$(samtools idxstats "$BAM" 2>/dev/null)
echo "[bamstats] idxstats: $(_elapsed $_T)"

mapped_reads=$(  awk '{m+=$3} END{print m+0}' <<< "$IDXSTATS")
unmapped_reads=$(awk '{u+=$4} END{print u+0}' <<< "$IDXSTATS")
total_reads=$(( mapped_reads + unmapped_reads ))

if [[ "$total_reads" -gt 0 ]]; then
    mapped_pct=$(awk "BEGIN{printf \"%.2f\", ${mapped_reads}/${total_reads}*100}")
else
    mapped_pct="0.00"
fi

# paired_reads: check one read in the chr1 window; WGS is always paired
paired_flag=$(samtools view "${_ST_VIEW_REF_FLAG[@]}" -f 0x1 -c -@ "$NCPU" "$BAM" "$STATS_REGION" 2>/dev/null || echo 0)
paired_reads=$(( paired_flag > 0 ? mapped_reads : 0 ))

# ── 2. samtools stats — 4 Mb representative chr1 window (~15 sec) ─────────────
# Note: in samtools stats the region is a positional argument after the BAM,
# NOT the -r flag (which specifies a reference FASTA, not a genomic region).
# A 4 Mb window at chr1:10M-14M avoids telomere/centromere and gives stable
# estimates of read length, insert size, and error rate much faster than chr1
# in full (~248 Mb scan reduced to a ~4 Mb scan).
echo "[bamstats] Running samtools stats on chr1:10M-14M …"
_T=$(_t0)
STATS=$(samtools stats "${_ST_STATS_REF_FLAG[@]}" -@ "$NCPU" "$BAM" "$STATS_REGION" 2>/dev/null)
if ! grep -q "^SN" <<< "$STATS" 2>/dev/null; then
    # Fallback: non-chr-prefixed reference (e.g. "1" instead of "chr1")
    STATS=$(samtools stats "${_ST_STATS_REF_FLAG[@]}" -@ "$NCPU" "$BAM" "1:10000000-14000000" 2>/dev/null)
fi
echo "[bamstats] samtools stats: $(_elapsed $_T)"

read_length=$(     grep "^SN	average length:"                  <<< "$STATS" | awk '{print int($NF)}')
insert_size_mean=$(grep "^SN	insert size average:"             <<< "$STATS" | awk '{print int($NF)}')
insert_size_sd=$(  grep "^SN	insert size standard deviation:"  <<< "$STATS" | awk '{print int($NF)}')
error_rate=$(      grep "^SN	error rate:"                      <<< "$STATS" | awk '{printf "%.5f", $NF}')

read_length="${read_length:-0}"
insert_size_mean="${insert_size_mean:-0}"
insert_size_sd="${insert_size_sd:-0}"
error_rate="${error_rate:-0}"

# ── 3. Global depth + sex inference — from idxstats + read_length ────────────
# Avoids a genome-wide mosdepth scan entirely.
# genome_depth = (mapped_reads × read_length) / sum(chrom_lengths)
# chrX/Y depth  = (chrX/Y_mapped × read_length) / chrX/Y_length
genome_size=$(awk '{s+=$2} END{print s+0}' <<< "$IDXSTATS")
chrX_mapped=$( awk '($1=="chrX"||$1=="X")  {print $3+0; exit}' <<< "$IDXSTATS")
chrY_mapped=$( awk '($1=="chrY"||$1=="Y")  {print $3+0; exit}' <<< "$IDXSTATS")
chrX_len=$(    awk '($1=="chrX"||$1=="X")  {print $2+0; exit}' <<< "$IDXSTATS")
chrY_len=$(    awk '($1=="chrY"||$1=="Y")  {print $2+0; exit}' <<< "$IDXSTATS")

mean_depth_genome=$(awk "BEGIN{
    if (${genome_size}+0 > 0) printf \"%.1f\", ${mapped_reads}*${read_length}/${genome_size};
    else print \"0.0\"}")
chrX_depth=$(awk "BEGIN{if (${chrX_len}+0>0) printf \"%.3f\", ${chrX_mapped}*${read_length}/${chrX_len}; else print \"0\"}")
chrY_depth=$(awk "BEGIN{if (${chrY_len}+0>0) printf \"%.3f\", ${chrY_mapped}*${read_length}/${chrY_len}; else print \"0\"}")
chrX_depth="${chrX_depth:-0}"; chrY_depth="${chrY_depth:-0}"

if awk "BEGIN{exit (${chrY_depth}+0 > 0) ? 0 : 1}"; then
    xy_ratio=$(   awk "BEGIN{printf \"%.2f\", ${chrX_depth}/${chrY_depth}}")
    inferred_sex=$(awk "BEGIN{print (${chrX_depth}/${chrY_depth} > 5) ? \"FEMALE\" : \"MALE\"}")
else
    xy_ratio="inf"
    inferred_sex="FEMALE"
fi

# ── 4. mosdepth — per-gene depth + thresholds ────────────────────────────────
# --by regions.bed     : index-guided access to each PGx region (no genome scan)
# --thresholds 0,20,30 : bases covered ≥N× per region → thresholds.bed.gz
# --fast-mode          : skip cigar-op parsing + mate-overlap correction (~20% faster)
# --no-per-base        : skip per-base BED output (not needed; large file)
# Note: pre-extracting a subset BAM was tested but slower overall — samtools view
# re-compression (~558 sec) outweighs the mosdepth speedup on the subset (4 sec).
MSD="${TMP}/pgx"

echo "[bamstats] Running mosdepth …"
_T=$(_t0)
mosdepth \
    "${_MSD_REF_FLAG[@]}"  \
    --by        "$PGX_BED" \
    --thresholds 0,20,30   \
    --no-per-base          \
    --fast-mode            \
    --threads   "$NCPU"    \
    "$MSD"                 \
    "$BAM"
echo "[bamstats] mosdepth: $(_elapsed $_T)"

# Per-gene mean depth from regions BED
# Columns: chrom  start  end  name  mean_depth
declare -A GENE_MEAN=()
while IFS=$'\t' read -r chrom start end name mean; do
    GENE_MEAN["$name"]="$mean"
done < <(zcat "${MSD}.regions.bed.gz" 2>/dev/null)

# Per-gene ≥20× / ≥30× coverage fractions from thresholds BED
# Columns: chrom  start  end  name  bases_ge_0x  bases_ge_20x  bases_ge_30x
# First line is a header (#chrom…) — skip it; fractions = bases_ge_Nx / region_length
declare -A GENE_PCT20=() GENE_PCT30=()
while IFS=$'\t' read -r chrom start end name b0 b20 b30; do
    [[ "$chrom" == \#* ]] && continue
    region_len=$(( end - start ))
    if [[ $region_len -gt 0 ]]; then
        GENE_PCT20["$name"]=$(awk "BEGIN{printf \"%.2f\", ${b20:-0}/${region_len}*100}")
        GENE_PCT30["$name"]=$(awk "BEGIN{printf \"%.2f\", ${b30:-0}/${region_len}*100}")
    fi
done < <(zcat "${MSD}.thresholds.bed.gz" 2>/dev/null)

# ── 5. MAPQ≥20 + duplicate estimate — chr1 4 Mb window (~500 ms) ─────────────
# Three fast samtools view -c calls on the same small region.
# Duplicate rate estimated from chr1 window (FLAG 0x400 set by MarkDuplicates).
echo "[bamstats] Computing MAPQ≥20 and duplicate fraction from chr1:10M-14M …"
_T=$(_t0)
chr1_total=$(  samtools view "${_ST_VIEW_REF_FLAG[@]}" -c          -@ "$NCPU" "$BAM" "$STATS_REGION" 2>/dev/null || echo 0)
chr1_mapq20=$( samtools view "${_ST_VIEW_REF_FLAG[@]}" -c -q 20    -@ "$NCPU" "$BAM" "$STATS_REGION" 2>/dev/null || echo 0)
chr1_dups=$(   samtools view "${_ST_VIEW_REF_FLAG[@]}" -c -f 0x400 -@ "$NCPU" "$BAM" "$STATS_REGION" 2>/dev/null || echo 0)
echo "[bamstats] samtools view -c (×3): $(_elapsed $_T)"
chr1_total="${chr1_total:-0}"; chr1_mapq20="${chr1_mapq20:-0}"; chr1_dups="${chr1_dups:-0}"
if [[ "$chr1_total" -gt 0 ]]; then
    mapq20_pct=$(awk "BEGIN{printf \"%.2f\", ${chr1_mapq20}/${chr1_total}*100}")
    dup_pct=$(   awk "BEGIN{printf \"%.2f\", ${chr1_dups}/${chr1_total}*100}")
else
    mapq20_pct="0.00"; dup_pct="0.00"
fi
duplicates=$(awk "BEGIN{printf \"%d\", ${dup_pct}/100*${mapped_reads}}")

# ── 6. Build gene_depth JSON object ───────────────────────────────────────────
# GSTT1 is on an alt contig — inject a null entry with a note so the report
# can display "Alt contig" rather than a missing value.
gene_depth_json="{"
first=1
for gene in "${!GENE_MEAN[@]}"; do
    mean="${GENE_MEAN[$gene]:-0}"
    pct20="${GENE_PCT20[$gene]:-0.00}"
    pct30="${GENE_PCT30[$gene]:-0.00}"
    [[ $first -eq 0 ]] && gene_depth_json+=","
    gene_depth_json+="\"${gene}\":{\"mean\":${mean},\"pct_ge_20x\":${pct20},\"pct_ge_30x\":${pct30}}"
    first=0
done
[[ $first -eq 0 ]] && gene_depth_json+=","
gene_depth_json+="\"GSTT1\":{\"mean\":null,\"pct_ge_20x\":null,\"pct_ge_30x\":null,\"note\":\"alt_contig\"}"
gene_depth_json+="}"

# ── 7. Write JSON ──────────────────────────────────────────────────────────────
cat > "$JSON_OUT" << EOF
{
  "sample":            "${SAMPLE}",
  "bam":               "${BAM}",
  "total_reads":       ${total_reads},
  "mapped_reads":      ${mapped_reads},
  "mapped_pct":        ${mapped_pct},
  "duplicate_reads":   ${duplicates},
  "duplicate_pct":     ${dup_pct},
  "paired_reads":      ${paired_reads},
  "read_length":       ${read_length},
  "insert_size_mean":  ${insert_size_mean},
  "insert_size_sd":    ${insert_size_sd},
  "mean_depth_genome": ${mean_depth_genome},
  "mapq20_pct":        ${mapq20_pct},
  "error_rate":        ${error_rate},
  "inferred_sex":      "${inferred_sex}",
  "xy_depth_ratio":    "${xy_ratio}",
  "gene_depth":        ${gene_depth_json}
}
EOF

echo "[bamstats] BAM stats written to: ${JSON_OUT}"
