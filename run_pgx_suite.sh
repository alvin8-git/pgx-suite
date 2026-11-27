#!/usr/bin/env bash
# run_pgx_suite.sh — Host-side launcher for the pgx-suite Docker pipeline
#
# Runs every gene in docker/genes.tsv for one WGS sample. Default engine: Nextflow on the
# host (main.nf), one container per step. --engine snakemake runs the Snakefile inside one
# container instead; both give the same results (docs/nextflow-port-plan.md).
#
# Usage:
#   run_pgx_suite.sh <BAM|CRAM> [options]
#
# Examples:
#   ./run_pgx_suite.sh /data/samples/NA12878.bam
#   ./run_pgx_suite.sh /data/samples/NA12878.cram --output /data/pgx_out
#   ./run_pgx_suite.sh /data/samples/NA12878.bam  --jobs 8 --output /data/out
#   ./run_pgx_suite.sh /data/samples/NA12878.bam  --engine snakemake
#
# Requirements:
#   - Docker with --privileged support
#   - Java 17+ for the default Nextflow engine (./nextflow downloads itself on first use)
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
MEM_GB="${PGX_MEM_GB:-}"   # empty: the Nextflow default (--max_memory, 32 GB)
IMAGE="pgx-suite:latest"
OUTPUT_HOST=""
ENGINE="nextflow"
RESUME=""

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF
Usage: $(basename "$0") <BAM|CRAM> [options]

  <BAM|CRAM>         Input alignment file (must have .bai/.crai index alongside)

Options:
  --output DIR       Parent output directory on host; results go into DIR/<SAMPLE>/
                     (default parent: <repo>/results/)
  --jobs N           CPUs to use               (default: 4)
  --mem GB           nextflow: memory budget for all concurrent steps (default: 32; env PGX_MEM_GB)
  --image NAME       Docker image to use       (default: pgx-suite:latest)
  --engine NAME      nextflow (default) or snakemake
  --resume           nextflow only: reuse finished steps of an earlier run of the same sample
  -h, --help         Show this help

Outputs written to <output>/<SAMPLE>/:
  <SAMPLE>_pgx_report.html      Standalone single-file HTML report (all genes embedded)
  all_genes_summary.tsv         One line per gene: verdict, consensus diplotype and phenotype
  genes/<GENE>/                 Per-gene tool outputs, comparison TSV, detail JSON, logs/
  bam_stats.json                Coverage statistics
  provenance.json               Image, tool versions, reference
  pipeline_info/                nextflow only: trace, timeline and resource report

Fixed resource directories (expected under ${SCRIPT_DIR}/):
  StellarPGx/              StellarPGx Nextflow pipeline
  StellarPGx/containers/   SIF images (stellarpgx-dev.sif, optitype.sif)
  pypgx/pypgx-bundle/      PyPGx Beagle haplotype panel
  GRCh38/                  GRCh38 reference FASTA (hg38.fa + hg38.fa.fai)

Note: --privileged is required for Apptainer used by StellarPGx and OptiType.
      Every container runs as the invoking user (--user UID:GID), never as root, with a generated passwd/group (this user plus root, from getent) mounted read-only (Apptainer refuses
      an unknown uid) and HOME=/tmp inside the container. OG_DOCKER_USER overrides: "--user ..." as-is,
      "uid[:gid]" gets --user prepended; an empty string runs as root again (not recommended; outputs
      become root-owned). Files left root-owned by an earlier run stop the launch with the chown to run.
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
        --mem)    MEM_GB="$2";      shift 2 ;;
        --image)  IMAGE="$2";       shift 2 ;;
        --engine) ENGINE="$2";      shift 2 ;;
        --resume) RESUME="-resume"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) _err "Unknown argument: $1"; usage; exit 1 ;;
    esac
done

case "$ENGINE" in
    nextflow|snakemake) ;;
    *) _err "--engine must be nextflow or snakemake, not '${ENGINE}'"; exit 1 ;;
esac

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

if [[ "$ENGINE" == "nextflow" ]]; then
    if ! command -v java &>/dev/null; then
        _err "'java' not found in PATH (Nextflow needs Java 17+). Install it, or use --engine snakemake."; exit 1
    fi
    echo -e "  Java         ${C_GREEN}found${C_RESET}   ${C_DIM}$(java -version 2>&1 | head -1)${C_RESET}"
fi

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
║                PGx Suite — pharmacogenomics pipeline                 ║
╚══════════════════════════════════════════════════════════════════════╝${C_RESET}
  ${C_BOLD}Sample :${C_RESET}  ${C_YELLOW}${SAMPLE}${C_RESET}
  ${C_BOLD}Input  :${C_RESET}  ${C_DIM}${INPUT_HOST}${C_RESET}
  ${C_BOLD}Output :${C_RESET}  ${C_DIM}${OUTPUT_HOST}${C_RESET}
  ${C_BOLD}Image  :${C_RESET}  ${C_DIM}${IMAGE}${C_RESET}
  ${C_BOLD}Engine :${C_RESET}  ${ENGINE}
  ${C_BOLD}CPUs   :${C_RESET}  ${JOBS}
  ${C_BOLD}Started:${C_RESET}  ${START_TS}"

# ── Run containers as the invoking user, never root ───────────────────────────
# docker/run_as_user.sh is the source of truth (OG_DOCKER_USER handling, generated passwd/group, HOME).
# The generated id files live in a mktemp folder this run owns (under TMPDIR), removed on exit.
# shellcheck source=docker/run_as_user.sh
source "${SCRIPT_DIR}/docker/run_as_user.sh"
pgx_user_opts || { echo "ERROR: cannot set up the run-as-user options" >&2; exit 2; }
[[ -n "${PGX_IDS_DIR:-}" ]] && trap 'rm -rf "$PGX_IDS_DIR"' EXIT
PGX_USER_OPTS_VAL="$PGX_USER_OPTS"
export PGX_DOCKER_USER_OPTS="$PGX_USER_OPTS_VAL"
echo -e "  ${C_BOLD}User   :${C_RESET}  ${C_DIM}${PGX_USER_OPTS_VAL:-root (OG_DOCKER_USER empty)}${C_RESET}"

# ── Fail fast if earlier root runs left files this user cannot overwrite ──────
# Only when running as a user: with OG_DOCKER_USER empty (root) ownership is not a problem.
if [[ -n "$PGX_USER_OPTS_VAL" ]]; then
    _chk=("$OUTPUT_HOST")
    [[ "$ENGINE" == "nextflow" ]] && _chk+=("$(dirname "$OUTPUT_HOST")/.nextflow-run/${SAMPLE}")
    for _d in "${_chk[@]}"; do
        [[ -d "$_d" ]] || continue
        if [[ ! -w "$_d" || ! -x "$_d" ]]; then
            echo "ERROR: $_d is not writable by $(id -un) (uid $(id -u))." >&2
            echo "  One-off fix: sudo chown -R $(id -u):$(id -g) $_d" >&2
            exit 2
        fi
        _bad="$(find "$_d" -maxdepth 3 ! -user "$(id -u)" -print -quit 2>/dev/null)"
        if [[ -n "$_bad" ]]; then
            echo "ERROR: $_bad is not owned by uid $(id -u) (left by an earlier run as root?)." >&2
            echo "  The container now runs as you and cannot overwrite it." >&2
            echo "  One-off fix: sudo chown -R $(id -u):$(id -g) $_d" >&2
            exit 2
        fi
    done
fi

# ── Run the pipeline ──────────────────────────────────────────────────────────
RC=0
if [[ "$ENGINE" == "nextflow" ]]; then
    _stage "Launching pipeline (Nextflow, one container per step)"
    # Launch folder and work folder per sample, beside the results: -resume finds its cache there,
    # and concurrent samples do not share one.
    NF_HOME="$(dirname "$OUTPUT_HOST")/.nextflow-run/${SAMPLE}"
    mkdir -p "$NF_HOME"
    ( cd "$NF_HOME" && "${SCRIPT_DIR}/nextflow" run "${SCRIPT_DIR}/main.nf" -profile docker $RESUME \
        --bam "$INPUT_HOST" --sample "$SAMPLE" --ref "${REF_DIR}/hg38.fa" \
        --pgx_image "$IMAGE" --max_cpus "$JOBS" ${MEM_GB:+--max_mem_total "$MEM_GB"} --outdir "$(dirname "$OUTPUT_HOST")" \
        -work-dir "${NF_HOME}/work" ) || RC=$?
else
    _stage "Launching pipeline (Snakemake inside one container)"
    echo -e "  ${C_DIM}Container output follows — gene results stream in as they complete${C_RESET}\n"
    # shellcheck disable=SC2086  # PGX_USER_OPTS_VAL is deliberately word-split
    docker run \
        --privileged \
        $PGX_USER_OPTS_VAL \
        --tmpfs /pgx:mode=0777 \
        --rm \
        -v "${STELLARPGX_DIR}:/pgx/stellarpgx" \
        -v "${CONTAINERS_DIR}:/pgx/containers" \
        -v "${BUNDLE_DIR}:/pgx/bundle" \
        -v "${REF_DIR}:/pgx/ref:ro" \
        -v "${INPUT_DIR_HOST}:/pgx/data:ro" \
        -v "${OUTPUT_HOST}:/pgx/results" \
        "$IMAGE" \
        snakemake \
            -s /opt/pgx/Snakefile \
            --cores "$JOBS" \
            --config bam="$INPUT_CONTAINER" outdir=/pgx/results ref=/pgx/ref/hg38.fa sample="$SAMPLE" || RC=$?
fi

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
    ${C_BOLD_GREEN}▸${C_RESET} Summary TSV  ${C_DIM}${OUTPUT_HOST}/all_genes_summary.tsv${C_RESET}
    ${C_BOLD_GREEN}▸${C_RESET} Gene results ${C_DIM}${OUTPUT_HOST}/genes/<GENE>/${C_RESET}"

if [[ $RC -ne 0 ]]; then
    echo -e "\n  ${C_BOLD_YELLOW}Check the logs for details:${C_RESET}"
    if [[ "$ENGINE" == "nextflow" ]]; then
        echo -e "    ${C_DIM}${NF_HOME}/.nextflow.log${C_RESET}"
    else
        echo -e "    ${C_DIM}ls ${OUTPUT_HOST}/genes/*/logs/${C_RESET}"
    fi
fi

echo ""
exit $RC
