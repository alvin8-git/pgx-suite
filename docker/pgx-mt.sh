#!/usr/bin/env bash
# pgx-mt.sh — MT-RNR1 mitochondrial variant calling using mutserve
#
# Usage: pgx-mt.sh <BAM|CRAM> <SAMPLE> <OUTPUT_DIR> [--ref REF] [--threads N]
#
# Calls mitochondrial variants with mutserve, then reports every CPIC MT-RNR1
# allele (docker/vcf_check_variants.json, from the CPIC aminoglycoside/MT-RNR1
# allele table): increased risk = m.1095T>C, m.1494C>T, m.1555A>G; uncertain-risk
# SNVs are listed but never escalate. Heteroplasmy level is reported.
# ERROR is never NEGATIVE: a mutserve failure aborts (exit 1, no result JSON),
# and insufficient chrM coverage at an increased-risk site reads "not typed".
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
# mutserve 2.x CLI: BAM/CRAM inputs are POSITIONAL (there is no --input option;
# passing one makes mutserve exit with "Unknown option" and write no output).
log_status "Running mutserve v2 (--level 0.01, --threads ${THREADS})"
MUTSERVE_RAW="${OUTDIR}/${SAMPLE}_mutserve_raw.txt"
RESULT_JSON="${OUTDIR}/${SAMPLE}_mtrna1_result.json"
rm -f "$MUTSERVE_RAW" "$RESULT_JSON"

if ! java -Xmx2g -jar "$MUTSERVE_JAR" call "${TMP}/chrM_sorted.bam" \
        --output    "$MUTSERVE_RAW" \
        --reference "${TMP}/chrM_ref.fa" \
        --level     0.01 \
        --threads   "$THREADS" \
        > "${OUTDIR}/mutserve.log" 2>&1; then
    log_status "ERROR mutserve failed — see ${OUTDIR}/mutserve.log (no result written; MT-RNR1 stays untyped)"
    exit 1
fi
if [[ ! -s "$MUTSERVE_RAW" ]]; then
    log_status "ERROR mutserve produced no output file (${MUTSERVE_RAW})"
    exit 1
fi

# ── 4. Per-site coverage at every table position (same filters as mutserve) ──
VARTABLE="$(dirname "$(readlink -f "$0")")/vcf_check_variants.json"
[[ -f "$VARTABLE" ]] || { log_status "ERROR variant table not found: ${VARTABLE}"; exit 1; }
python3 -c "
import json
for v in json.load(open('${VARTABLE}'))['MT-RNR1']['variants']:
    print('chrM\t%d\t%d' % (v['pos']-1, v['pos']))" > "${TMP}/sites.bed"
samtools depth -a -q 20 -Q 20 -b "${TMP}/sites.bed" "${TMP}/chrM_sorted.bam" > "${TMP}/site_depth.tsv"

# ── 5. Parse output and write JSON result ────────────────────────────────────
log_status "Reporting CPIC MT-RNR1 alleles (increased: m.1095T>C, m.1494C>T, m.1555A>G)"

SAMPLE="$SAMPLE" MUTSERVE_RAW="$MUTSERVE_RAW" VARTABLE="$VARTABLE" DEPTH="${TMP}/site_depth.tsv" \
python3 - <<'PYEOF' > "${RESULT_JSON}.tmp"
import json, os

sample = os.environ["SAMPLE"]
tbl    = json.load(open(os.environ["VARTABLE"]))["MT-RNR1"]
min_cov = int(tbl.get("min_cov", 100))
by_key = {(v["pos"], v["ref"], v["alt"]): v for v in tbl["variants"]}  # 807 has two alts
INCREASED = sorted(v["pos"] for v in tbl["variants"] if v["cpic_class"] == "increased risk")

depth = {}
for line in open(os.environ["DEPTH"]):
    c, p, d = line.split("\t")
    depth[int(p)] = int(d)

# mutserve 2.x raw table: header row names the columns (ID, Filter, Pos, Ref,
# Variant, VariantLevel, ..., Coverage, ...).  Parse by name, never by index.
calls = []
with open(os.environ["MUTSERVE_RAW"]) as fh:
    rows = [l.rstrip("\n").split("\t") for l in fh if l.strip() and not l.startswith("#")]
if not rows or "Pos" not in rows[0]:
    raise SystemExit("mutserve output has no recognisable header: %r" % (rows[:1],))
hdr = {name: i for i, name in enumerate(rows[0])}
for need in ("Pos", "Ref", "Variant", "VariantLevel"):
    if need not in hdr:
        raise SystemExit("mutserve output lacks column %s" % need)
for r in rows[1:]:
    try:
        calls.append({"pos": int(r[hdr["Pos"]]), "ref": r[hdr["Ref"]].upper(),
                      "alt": r[hdr["Variant"]].upper(),
                      "level": float(r[hdr["VariantLevel"]]),
                      "filter": r[hdr["Filter"]] if "Filter" in hdr else "PASS",
                      "cov": int(r[hdr["Coverage"]]) if "Coverage" in hdr else None})
    except (ValueError, IndexError) as e:
        raise SystemExit("mutserve data row could not be parsed (%s): %r" % (e, r))

found, filtered = [], []
for c in calls:
    v = by_key.get((c["pos"], c["ref"], c["alt"]))
    if v and c["filter"] != "PASS":   # seen but filtered: never silently read as reference
        filtered.append((v, c))
    if v and c["filter"] == "PASS":
        found.append({"pos": c["pos"], "ref": c["ref"], "alt": c["alt"], "label": v["label"],
                      "af": round(c["level"], 4), "coverage": c["cov"] or depth.get(c["pos"]),
                      "type": "homoplasmic" if c["level"] >= 0.90 else "heteroplasmic",
                      "cpic_class": v["cpic_class"]})
found.sort(key=lambda x: x["pos"])
inc = [f for f in found if f["cpic_class"] == "increased risk"]
unc = [f for f in found if f["cpic_class"] == "uncertain risk"]
untyped = sorted(set([p for p in INCREASED if depth.get(p, 0) < min_cov] +
                     [v["pos"] for v, c in filtered if v["cpic_class"] == "increased risk"]))
filt_txt = "; ".join("%s filtered by mutserve (%s, level %.2f)" % (v["label"], c["filter"], c["level"])
                     for v, c in filtered)
cov_txt = ", ".join("m.%d %dx" % (p, depth.get(p, 0)) for p in INCREASED)

def fmt(f):
    return "%s (AF=%.1f%%, %s)" % (f["label"], f["af"] * 100, f["type"])

if inc:
    diplotype, phenotype, classification = ";".join(f["label"] for f in inc), \
        "Aminoglycoside-ototoxicity risk", "carrier"
    notes = ("CPIC increased-risk MT-RNR1 variant(s) detected: " + "; ".join(fmt(f) for f in inc) +
             ". CPIC: avoid aminoglycosides (gentamicin, tobramycin, amikacin, streptomycin, "
             "neomycin) unless the high risk of ototoxicity is outweighed by the severity of "
             "infection and safe alternatives are unavailable. Maternally inherited.")
    if untyped:
        notes += " Other increased-risk sites not typed (coverage < %dx): %s." % (
            min_cov, ", ".join("m.%d" % p for p in untyped))
elif untyped:
    diplotype, phenotype, classification = "not typed", \
        "Not typed - chrM coverage below %dx or a filtered call at increased-risk site(s) (%s)" % (min_cov, cov_txt), "not_typed"
    notes = ("MT-RNR1 could not be typed: depth at the CPIC increased-risk sites is below %dx, or a call "
             "there was filtered (%s). No result is not a negative result." % (min_cov, cov_txt))
else:
    diplotype, classification = "Reference", "non-carrier"
    phenotype = ("No CPIC increased-risk MT-RNR1 variant detected (m.1095T>C, m.1494C>T, m.1555A>G "
                 "typed at %s; heteroplasmy >= 1%%)" % cov_txt)
    notes = phenotype + ". Variants outside the CPIC table are not assessed."
if filt_txt:
    notes += " " + filt_txt + "."
if unc:
    notes += " CPIC uncertain-risk variant(s) observed (not actionable): " + "; ".join(fmt(f) for f in unc) + "."

print(json.dumps({
    "sample": sample, "gene": "MT-RNR1", "tool": "mutserve",
    "variants": [f for f in found if f["cpic_class"] != "normal risk"],
    "site_depth": {"m.%d" % p: depth.get(p, 0) for p in INCREASED},
    "min_cov": min_cov, "sites_checked": len({k[0] for k in by_key}),
    "diplotype": diplotype, "phenotype": phenotype,
    "classification": classification, "notes": notes,
}, indent=2))
PYEOF
mv "${RESULT_JSON}.tmp" "$RESULT_JSON"   # parser abort (set -e) leaves no result JSON

DIPLOTYPE=$(python3 -c "
import json, sys
try:
    print(json.load(open('${RESULT_JSON}'))['diplotype'])
except: print('?')
" 2>/dev/null)
log_status "MT-RNR1 result: ${DIPLOTYPE}"
log_status "Result written to ${RESULT_JSON}"
