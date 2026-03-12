#!/usr/bin/env bash
# run_pgx_suite.sh — Host-side launcher for the pgx-suite Docker pipeline
#
# Runs the full 31-gene PGx star-allele + HLA typing + mtDNA pipeline for a single
# WGS sample by launching the pgx-suite Docker container with appropriate
# volume mounts, forwarding to pgx-all-genes.sh inside.
#
# Usage:
#   run_pgx_suite.sh <BAM|CRAM> [options]
#
# Examples:
#   ./run_pgx_suite.sh /data/samples/NA12878.bam
#   ./run_pgx_suite.sh /data/samples/NA12878.cram --output /data/pgx_out
#   ./run_pgx_suite.sh /data/samples/NA12878.bam  --jobs 8 --output /data/out
#
# Requirements:
#   - Docker with --privileged support
#   - pgx-suite:latest image (docker build -t pgx-suite:latest .)
#   - Resource directories under this repo (see PATHS section below)
#
# Fixed resource directories (relative to this script):
#   StellarPGx/            StellarPGx Nextflow pipeline + databases
#   StellarPGx/containers/ SIF images (stellarpgx-dev.sif, optitype.sif)
#   pypgx/pypgx-bundle/    PyPGx Beagle haplotype panel (~500 MB)
#   GRCh38/                GRCh38 reference FASTA (hg38.fa + hg38.fa.fai)
# mutserve.jar (MT-RNR1) is baked into the Docker image — no extra mount needed.

set -euo pipefail

# ── ANSI colours (disabled automatically when stdout is not a terminal) ────────
if [[ -t 1 ]]; then
    C_RESET='\033[0m'
    C_BOLD='\033[1m'
    C_DIM='\033[2m'
    C_RED='\033[0;31m'
    C_GREEN='\033[0;32m'
    C_YELLOW='\033[0;33m'
    C_BLUE='\033[0;34m'
    C_CYAN='\033[0;36m'
    C_WHITE='\033[0;37m'
    C_BOLD_RED='\033[1;31m'
    C_BOLD_GREEN='\033[1;32m'
    C_BOLD_YELLOW='\033[1;33m'
    C_BOLD_CYAN='\033[1;36m'
else
    C_RESET='' C_BOLD='' C_DIM=''
    C_RED='' C_GREEN='' C_YELLOW='' C_BLUE='' C_CYAN='' C_WHITE=''
    C_BOLD_RED='' C_BOLD_GREEN='' C_BOLD_YELLOW='' C_BOLD_CYAN=''
fi

# ── Helpers ───────────────────────────────────────────────────────────────────
_err()   { echo -e "${C_BOLD_RED}ERROR:${C_RESET} $*" >&2; }
_warn()  { echo -e "${C_BOLD_YELLOW}WARN:${C_RESET}  $*"; }
_ok()    { echo -e "${C_BOLD_GREEN}OK${C_RESET}"; }
_stage() { echo -e "\n${C_BOLD_CYAN}▶ $*${C_RESET}"; }

_elapsed() {
    local secs=$(( $(date +%s) - START_EPOCH ))
    printf '%02d:%02d' $(( secs / 60 )) $(( secs % 60 ))
}

# ── Locate repo root (directory containing this script) ───────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Fixed resource paths (host-side) ─────────────────────────────────────────
STELLARPGX_DIR="${SCRIPT_DIR}/StellarPGx"
CONTAINERS_DIR="${SCRIPT_DIR}/StellarPGx/containers"
BUNDLE_DIR="${SCRIPT_DIR}/pypgx/pypgx-bundle"
REF_DIR="${SCRIPT_DIR}/GRCh38"

# ── Defaults ──────────────────────────────────────────────────────────────────
JOBS=4
IMAGE="pgx-suite:latest"
OUTPUT_HOST=""

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF
Usage: $(basename "$0") <BAM|CRAM> [options]

  <BAM|CRAM>         Input alignment file (must have .bai/.crai index alongside)

Options:
  --output DIR       Parent output directory on host; results go into DIR/<SAMPLE>/
                     (default parent: <repo>/results/)
  --jobs N           Genes to run in parallel  (default: 4)
  --image NAME       Docker image to use       (default: pgx-suite:latest)
  -h, --help         Show this help

Outputs written to <output>/<SAMPLE>/:
  <SAMPLE>_pgx_report.html      Standalone single-file HTML report (all genes embedded)
  Genes/<GENE>/                 Per-gene tool outputs, comparison TSV, detail JSON
  log/all_genes_summary.tsv     Concordance summary across all 31 genes
  log/bam_stats.json            Coverage statistics (mosdepth per-gene depth)
  log/<GENE>.log                Per-gene stdout/stderr logs

Fixed resource directories (expected under ${SCRIPT_DIR}/):
  StellarPGx/              StellarPGx Nextflow pipeline
  StellarPGx/containers/   SIF images (stellarpgx-dev.sif, optitype.sif)
  pypgx/pypgx-bundle/      PyPGx Beagle haplotype panel
  GRCh38/                  GRCh38 reference FASTA (hg38.fa + hg38.fa.fai)

Note: --privileged is required for Apptainer used by StellarPGx and OptiType.
      mutserve.jar (MT-RNR1 chrM caller) is baked into the image — no extra mount needed.
EOF
}

# ── Argument parsing ───────────────────────────────────────────────────────────
if [[ $# -lt 1 ]]; then usage; exit 1; fi

INPUT_HOST="$1"; shift

while [[ $# -gt 0 ]]; do
    case "$1" in
        --output) OUTPUT_HOST="$2"; shift 2 ;;
        --jobs)   JOBS="$2";        shift 2 ;;
        --image)  IMAGE="$2";       shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) _err "Unknown argument: $1"; usage; exit 1 ;;
    esac
done

# ── Resolve to absolute path ──────────────────────────────────────────────────
INPUT_HOST="$(realpath "$INPUT_HOST")"

# ── Detect BAM vs CRAM ────────────────────────────────────────────────────────
EXT="${INPUT_HOST##*.}"
case "$EXT" in
    bam)  INDEX_EXT="bai" ;;
    cram) INDEX_EXT="crai" ;;
    *)    _err "Input must be a .bam or .cram file: ${INPUT_HOST}"; exit 1 ;;
esac

# ── Validate inputs ───────────────────────────────────────────────────────────
_stage "Validating inputs"

if [[ ! -f "$INPUT_HOST" ]]; then
    _err "Input file not found: ${INPUT_HOST}"; exit 1
fi
echo -e "  Input file   ${C_GREEN}found${C_RESET}   ${C_DIM}${INPUT_HOST}${C_RESET}"

# Accept <file>.bam.bai or <file>.bai
if [[ ! -f "${INPUT_HOST}.${INDEX_EXT}" && ! -f "${INPUT_HOST%.${EXT}}.${INDEX_EXT}" ]]; then
    _err "Index not found (expected ${INPUT_HOST}.${INDEX_EXT})"
    echo "       Run: samtools index ${INPUT_HOST}" >&2
    exit 1
fi
echo -e "  Index file   ${C_GREEN}found${C_RESET}"

for dir in "$STELLARPGX_DIR" "$CONTAINERS_DIR" "$BUNDLE_DIR" "$REF_DIR"; do
    if [[ ! -d "$dir" ]]; then
        _err "Required resource directory not found: ${dir}"; exit 1
    fi
    echo -e "  Resource dir ${C_GREEN}found${C_RESET}   ${C_DIM}${dir}${C_RESET}"
done

if ! command -v docker &>/dev/null; then
    _err "'docker' not found in PATH"; exit 1
fi
echo -e "  Docker       ${C_GREEN}found${C_RESET}   ${C_DIM}$(docker --version)${C_RESET}"

if ! docker image inspect "$IMAGE" &>/dev/null 2>&1; then
    _err "Docker image '${IMAGE}' not found."
    echo "       Build it with: docker build -t ${IMAGE} ${SCRIPT_DIR}/" >&2
    exit 1
fi
echo -e "  Image        ${C_GREEN}found${C_RESET}   ${C_DIM}${IMAGE}${C_RESET}"

# ── Derive sample name and output path ────────────────────────────────────────
SAMPLE="$(basename "$INPUT_HOST" | cut -d'.' -f1)"

# Always place results in a <SAMPLE> subdirectory so that:
#   --output ./results/  →  ./results/HG005/
#   (no --output)        →  <repo>/results/HG005/
if [[ -z "$OUTPUT_HOST" ]]; then
    OUTPUT_HOST="${SCRIPT_DIR}/results/${SAMPLE}"
else
    OUTPUT_HOST="${OUTPUT_HOST%/}/${SAMPLE}"
fi

mkdir -p "$OUTPUT_HOST"
OUTPUT_HOST="$(realpath "$OUTPUT_HOST")"

# ── Container-side paths ──────────────────────────────────────────────────────
INPUT_DIR_HOST="$(dirname "$INPUT_HOST")"
INPUT_FILENAME="$(basename "$INPUT_HOST")"
INPUT_CONTAINER="/pgx/data/${INPUT_FILENAME}"

# ── Start timer ───────────────────────────────────────────────────────────────
START_EPOCH=$(date +%s)
START_TS=$(date '+%Y-%m-%d %H:%M:%S')

# ── Banner ────────────────────────────────────────────────────────────────────
echo -e "
${C_BOLD_CYAN}╔══════════════════════════════════════════════════════════════════════╗
║              PGx Suite — 31-gene pharmacogenomics pipeline           ║
╚══════════════════════════════════════════════════════════════════════╝${C_RESET}
  ${C_BOLD}Sample :${C_RESET}  ${C_YELLOW}${SAMPLE}${C_RESET}
  ${C_BOLD}Input  :${C_RESET}  ${C_DIM}${INPUT_HOST}${C_RESET}
  ${C_BOLD}Output :${C_RESET}  ${C_DIM}${OUTPUT_HOST}${C_RESET}
  ${C_BOLD}Image  :${C_RESET}  ${C_DIM}${IMAGE}${C_RESET}
  ${C_BOLD}Jobs   :${C_RESET}  ${JOBS} genes in parallel
  ${C_BOLD}Started:${C_RESET}  ${START_TS}"

# ── Run the Docker container ──────────────────────────────────────────────────
_stage "Launching pipeline (pgx-all-genes.sh inside container)"
echo -e "  ${C_DIM}Container output follows — gene results stream in as they complete${C_RESET}\n"

docker run \
    --privileged \
    --rm \
    -v "${STELLARPGX_DIR}:/pgx/stellarpgx" \
    -v "${CONTAINERS_DIR}:/pgx/containers" \
    -v "${BUNDLE_DIR}:/pgx/bundle" \
    -v "${REF_DIR}:/pgx/ref:ro" \
    -v "${INPUT_DIR_HOST}:/pgx/data:ro" \
    -v "${OUTPUT_HOST}:/pgx/results" \
    "$IMAGE" \
    pgx-all-genes.sh \
        "$INPUT_CONTAINER" \
        --ref    /pgx/ref/hg38.fa \
        --output /pgx/results \
        --jobs   "$JOBS"

RC=$?

# ── Final summary ─────────────────────────────────────────────────────────────
ELAPSED=$(_elapsed)
END_TS=$(date '+%Y-%m-%d %H:%M:%S')

echo ""
if [[ $RC -eq 0 ]]; then
    STATUS_LINE="${C_BOLD_GREEN}✔  Pipeline completed successfully${C_RESET}"
else
    STATUS_LINE="${C_BOLD_RED}✘  Pipeline finished with errors  (exit code: ${RC})${C_RESET}"
fi

echo -e "${C_BOLD_CYAN}╔══════════════════════════════════════════════════════════════════════╗
║                          Run Summary                                 ║
╚══════════════════════════════════════════════════════════════════════╝${C_RESET}
  ${STATUS_LINE}

  ${C_BOLD}Sample  :${C_RESET}  ${C_YELLOW}${SAMPLE}${C_RESET}
  ${C_BOLD}Duration:${C_RESET}  ${C_BOLD_YELLOW}${ELAPSED}${C_RESET}  (${START_TS} → ${END_TS})

  ${C_BOLD}Outputs :${C_RESET}
    ${C_BOLD_GREEN}▸${C_RESET} HTML report  ${C_DIM}${OUTPUT_HOST}/${SAMPLE}_pgx_report.html${C_RESET}
    ${C_BOLD_GREEN}▸${C_RESET} Summary TSV  ${C_DIM}${OUTPUT_HOST}/log/all_genes_summary.tsv${C_RESET}
    ${C_BOLD_GREEN}▸${C_RESET} Gene results ${C_DIM}${OUTPUT_HOST}/Genes/<GENE>/${C_RESET}
    ${C_BOLD_GREEN}▸${C_RESET} Logs         ${C_DIM}${OUTPUT_HOST}/log/<GENE>.log${C_RESET}"

if [[ $RC -ne 0 ]]; then
    echo -e "\n  ${C_BOLD_YELLOW}Check per-gene logs for details:${C_RESET}"
    echo -e "    ${C_DIM}ls ${OUTPUT_HOST}/log/*.log${C_RESET}"
fi

echo ""
exit $RC
