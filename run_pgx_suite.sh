#!/usr/bin/env bash
# run_pgx_suite.sh — Host-side launcher for the pgx-suite Docker pipeline
#
# Runs the full 29-gene PGx star-allele + HLA typing pipeline for a single
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
#   - pgx-suite:latest image (docker build -t pgx-suite:latest docker/)
#   - Resource directories under this repo (see PATHS section below)
#
# Fixed resource directories (relative to this script):
#   StellarPGx/            StellarPGx Nextflow pipeline + databases
#   StellarPGx/containers/ SIF images (stellarpgx-dev.sif, optitype.sif)
#   pypgx/pypgx-bundle/    PyPGx Beagle haplotype panel (~500 MB)
#   GRCh38/                GRCh38 reference FASTA (hg38.fa) + index

set -euo pipefail

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
  --output DIR       Output directory on host  (default: <repo>/results/<SAMPLE>)
  --jobs N           Genes to run in parallel  (default: 4)
  --image NAME       Docker image to use       (default: pgx-suite:latest)
  -h, --help         Show this help

Outputs written to <output>/:
  <GENE>/                  Per-gene results (comparison TSV + tool subdirs)
  logs/                    Per-gene stdout/stderr logs
  bam_stats.json           Coverage statistics across all 27 primary-assembly genes
  all_genes_summary.tsv    Concordance summary across all 29 genes
  html_reports/            HTML report (<SAMPLE>.pgx.html + per-gene pages)

Fixed resource directories (expected under ${SCRIPT_DIR}/):
  StellarPGx/              StellarPGx Nextflow pipeline
  StellarPGx/containers/   SIF images (stellarpgx-dev.sif, optitype.sif)
  pypgx/pypgx-bundle/      PyPGx Beagle haplotype panel
  GRCh38/                  GRCh38 reference FASTA (hg38.fa)

Note: --privileged is required for Apptainer used by StellarPGx and OptiType.
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
        *) echo "ERROR: Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

# ── Resolve to absolute path ──────────────────────────────────────────────────
INPUT_HOST="$(realpath "$INPUT_HOST")"

# ── Detect BAM vs CRAM ────────────────────────────────────────────────────────
EXT="${INPUT_HOST##*.}"
case "$EXT" in
    bam)  INDEX_EXT="bai" ;;
    cram) INDEX_EXT="crai" ;;
    *)    echo "ERROR: Input must be a .bam or .cram file: ${INPUT_HOST}" >&2; exit 1 ;;
esac

# ── Validate inputs ───────────────────────────────────────────────────────────
if [[ ! -f "$INPUT_HOST" ]]; then
    echo "ERROR: Input file not found: ${INPUT_HOST}" >&2; exit 1
fi

# Accept <file>.bam.bai or <file>.bai
if [[ ! -f "${INPUT_HOST}.${INDEX_EXT}" && ! -f "${INPUT_HOST%.${EXT}}.${INDEX_EXT}" ]]; then
    echo "ERROR: Index not found (expected ${INPUT_HOST}.${INDEX_EXT})" >&2
    echo "       Run: samtools index ${INPUT_HOST}" >&2
    exit 1
fi

for dir in "$STELLARPGX_DIR" "$CONTAINERS_DIR" "$BUNDLE_DIR" "$REF_DIR"; do
    if [[ ! -d "$dir" ]]; then
        echo "ERROR: Required resource directory not found: ${dir}" >&2
        exit 1
    fi
done

if ! command -v docker &>/dev/null; then
    echo "ERROR: 'docker' not found in PATH" >&2; exit 1
fi

if ! docker image inspect "$IMAGE" &>/dev/null 2>&1; then
    echo "ERROR: Docker image '${IMAGE}' not found." >&2
    echo "       Build it with: docker build -t ${IMAGE} ${SCRIPT_DIR}/docker/" >&2
    exit 1
fi

# ── Derive sample name and output path ────────────────────────────────────────
SAMPLE="$(basename "$INPUT_HOST" | cut -d'.' -f1)"

if [[ -z "$OUTPUT_HOST" ]]; then
    OUTPUT_HOST="${SCRIPT_DIR}/results/${SAMPLE}"
fi

mkdir -p "$OUTPUT_HOST"
OUTPUT_HOST="$(realpath "$OUTPUT_HOST")"

# ── Container-side paths ──────────────────────────────────────────────────────
INPUT_DIR_HOST="$(dirname "$INPUT_HOST")"
INPUT_FILENAME="$(basename "$INPUT_HOST")"
INPUT_CONTAINER="/pgx/data/${INPUT_FILENAME}"

# ── Banner ────────────────────────────────────────────────────────────────────
cat <<EOF
========================================================================
 PGx Suite — Host launcher
 Sample:  ${SAMPLE}
 Input:   ${INPUT_HOST}
 Output:  ${OUTPUT_HOST}
 Image:   ${IMAGE}
 Jobs:    ${JOBS} genes in parallel
 Started: $(date '+%Y-%m-%d %H:%M:%S')
========================================================================
EOF

# ── Run the Docker container ──────────────────────────────────────────────────
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

cat <<EOF

========================================================================
 Run complete — exit code: ${RC}
 HTML report: ${OUTPUT_HOST}/html_reports/${SAMPLE}.pgx.html
 Summary TSV: ${OUTPUT_HOST}/all_genes_summary.tsv
 Finished:    $(date '+%Y-%m-%d %H:%M:%S')
========================================================================
EOF

exit $RC
