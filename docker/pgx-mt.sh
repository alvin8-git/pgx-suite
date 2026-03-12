#!/usr/bin/env bash
# pgx-mt.sh — MT-RNR1 mitochondrial variant calling using mutserve
#
# Usage: pgx-mt.sh <BAM|CRAM> <SAMPLE> <OUTPUT_DIR> [--ref REF] [--threads N]
#
# Calls mitochondrial variants with mutserve, then reports CPIC Level A
# aminoglycoside-ototoxicity variants: m.1555A>G and m.1494C>T.
#
# Output: <OUTPUT_DIR>/mt-rnr1/<SAMPLE>_mtrna1_result.json
#
# Requires:
#   - mutserve.jar at /usr/local/bin/mutserve.jar
#   - Java 21 JRE (already in image)
#   - samtools (already in image)
#   - /pgx/ref/hg38.fa (volume-mounted)

set -euo pipefail

SCRIPT=$(basename "$0")
log_status() { echo "[$(date '+%H:%M:%S')] ${SCRIPT}: $*" >&2; }

# ── Defaults ─────────────────────────────────────────────────────────────────
REF="/pgx/ref/hg38.fa"
THREADS=4

# ── Arguments ────────────────────────────────────────────────────────────────
BAM="${1:-}"
SAMPLE="${2:-}"
OUTPUT="${3:-/pgx/results}"
shift 3 || true

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)     REF="$2";     shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$BAM" || -z "$SAMPLE" ]]; then
    echo "Usage: pgx-mt.sh <BAM|CRAM> <SAMPLE> [OUTPUT] [--ref REF] [--threads N]" >&2
    exit 1
fi

MUTSERVE_JAR="/usr/local/bin/mutserve.jar"
if [[ ! -f "$MUTSERVE_JAR" ]]; then
    log_status "ERROR mutserve.jar not found at ${MUTSERVE_JAR}"
    exit 1
fi

OUTDIR="${OUTPUT}/mt-rnr1"
mkdir -p "$OUTDIR"
TMP=$(mktemp -d "${OUTDIR}/tmp_XXXXXX")
trap 'rm -rf "$TMP"' EXIT

# ── 1. Extract chrM reads ─────────────────────────────────────────────────────
log_status "Extracting chrM reads from ${BAM}"
samtools view -b -@ "$THREADS" "$BAM" chrM \
    | samtools sort -@ "$THREADS" -o "${TMP}/chrM_sorted.bam" -
samtools index "${TMP}/chrM_sorted.bam"

MAPPED=$(samtools flagstat "${TMP}/chrM_sorted.bam" \
    | awk '/^[0-9]+ \+ [0-9]+ mapped/{print $1; exit}')
log_status "chrM mapped reads: ${MAPPED:-unknown}"

# ── 2. Extract chrM reference ─────────────────────────────────────────────────
log_status "Extracting chrM reference from ${REF}"
samtools faidx "$REF" chrM > "${TMP}/chrM_ref.fa"
samtools faidx "${TMP}/chrM_ref.fa"

# ── 3. Run mutserve ───────────────────────────────────────────────────────────
log_status "Running mutserve v2 (--level 0.01, --threads ${THREADS})"
MUTSERVE_RAW="${OUTDIR}/${SAMPLE}_mutserve_raw.txt"

java -Xmx2g -jar "$MUTSERVE_JAR" call \
    --input     "${TMP}/chrM_sorted.bam" \
    --output    "$MUTSERVE_RAW" \
    --reference "${TMP}/chrM_ref.fa" \
    --level     0.01 \
    --threads   "$THREADS" \
    > "${OUTDIR}/mutserve.log" 2>&1 \
    || log_status "WARN mutserve exited non-zero — check ${OUTDIR}/mutserve.log"

# ── 4. Parse output and write JSON result ────────────────────────────────────
log_status "Parsing mutserve output for m.1555A>G (pos 1555) and m.1494C>T (pos 1494)"

RESULT_JSON="${OUTDIR}/${SAMPLE}_mtrna1_result.json"

SAMPLE="$SAMPLE" MUTSERVE_RAW="$MUTSERVE_RAW" \
python3 - <<'PYEOF' > "$RESULT_JSON"
import json, os

sample      = os.environ["SAMPLE"]
raw_file    = os.environ["MUTSERVE_RAW"]

# CPIC Level A target variants (chrM rCRS positions)
TARGETS = {
    (1555, "A", "G"): "m.1555A>G",
    (1494, "C", "T"): "m.1494C>T",
}

found = []
if os.path.exists(raw_file):
    with open(raw_file) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            # First non-comment line: detect header by checking if Pos column is text
            if header is None and not parts[1].isdigit():
                header = {v.strip(): i for i, v in enumerate(parts)}
                continue
            try:
                pos   = int(parts[1])
                ref   = parts[2].upper()
                alt   = parts[3].upper()
                pct_v = float(parts[5]) if len(parts) > 5 else 0.0
                af    = round(pct_v / 100.0, 4)
                vtype = parts[7].lower() if len(parts) > 7 else "unknown"
            except (ValueError, IndexError):
                continue
            label = TARGETS.get((pos, ref, alt))
            if label:
                found.append({"pos": pos, "ref": ref, "alt": alt,
                               "af": af, "type": vtype, "label": label})

if found:
    diplotype      = ";".join(v["label"] for v in found)
    phenotype      = "Aminoglycoside-ototoxicity risk"
    classification = "carrier"
    notes = (
        f"CPIC Level A variant(s) detected: {diplotype}. "
        "Avoid all aminoglycosides (gentamicin, tobramycin, amikacin, "
        "streptomycin, neomycin) unless no safe alternatives exist."
    )
else:
    diplotype      = "Reference"
    phenotype      = "Standard aminoglycoside risk"
    classification = "non-carrier"
    notes = ("No CPIC Level A MT-RNR1 variants detected (m.1555A>G, m.1494C>T). "
             "Standard aminoglycoside prescribing applies.")

print(json.dumps({
    "sample":         sample,
    "gene":           "MT-RNR1",
    "tool":           "mutserve",
    "variants":       found,
    "diplotype":      diplotype,
    "phenotype":      phenotype,
    "classification": classification,
    "notes":          notes,
}, indent=2))
PYEOF

DIPLOTYPE=$(python3 -c "
import json, sys
try:
    print(json.load(open('${RESULT_JSON}'))['diplotype'])
except: print('?')
" 2>/dev/null)
log_status "MT-RNR1 result: ${DIPLOTYPE}"
log_status "Result written to ${RESULT_JSON}"
