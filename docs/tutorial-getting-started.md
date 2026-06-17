# Tutorial: from a BAM to your first PGx report

By the end of this you will have built the pgx-suite image and produced a
self-contained HTML pharmacogenomics report for one WGS sample — all 31 genes,
every caller, reconciled into one verdict per gene.

This is the learning path. For the reasoning behind each piece see the
[README](../README.md); for a specific task (e.g. adding a gene) see the
[how-to guides](howto-add-a-gene.md).

## What you'll need

- **Docker** with `--privileged` support (a local workstation, not a locked-down CI runner).
  StellarPGx and OptiType run Apptainer inside the container, which needs it.
- The repository, with submodules cloned:
  ```bash
  git clone --recurse-submodules https://github.com/alvin8-git/pgx-suite.git
  cd pgx-suite
  ```
- A **GRCh38 reference FASTA** with a `.fai` index, `chr`-prefixed contigs, and the
  `chr22_KI270879v1_alt` contig (needed for GSTT1). e.g. `hg38.fa` + `hg38.fa.fai`.
- A **WGS BAM or CRAM** aligned to that reference, with its `.bai`/`.crai` index.
- The bundled large data, already in the repo after cloning:
  `StellarPGx/containers/stellarpgx-dev.sif`, `pypgx/pypgx-bundle/`. Pull the OptiType
  SIF once if you want HLA typing (see the README "Installation" section).
- ~3 GB free disk for the image; the first build takes ~15-25 min.

## Step 1: Build the image

```bash
docker build -t pgx-suite:latest .
```

You'll see a multi-stage build: a builder stage compiling the Python venv (pysam,
pandas, Snakemake, …), then a lean runtime stage. When it finishes you have a
`pgx-suite:latest` image with all six callers + the Snakemake orchestrator baked in.

Confirm it:

```bash
docker run --rm pgx-suite:latest snakemake --version
```

You should see a Snakemake version number (e.g. `9.x`).

## Step 2: Run the pipeline on your sample

The host launcher wires up every Docker mount for you. Point it at your BAM:

```bash
./run_pgx_suite.sh /path/to/sample.bam \
    --ref /path/to/GRCh38/hg38.fa \
    --output results/my_sample \
    --jobs 4
```

What happens: the launcher starts the container, Snakemake expands one rule graph
over all 31 genes (`--jobs` controls how many run at once), each caller runs as a
non-fatal rule, `pgx-compare.py` reconciles them into a verdict per gene, BAM QC runs,
and `pgx-report.py` writes the HTML. On a 30-50× WGS this takes a few minutes; the
gene results stream to the console as they complete.

> A single caller failing (StellarPGx is the usual one) does **not** stop the run —
> that gene is still reported, with the failed tool marked `failed`.

## Step 3: Open the report

```bash
# the report is named after the sample (from the BAM @RG SM tag)
ls results/my_sample/*_pgx_report.html
```

Open that file in any browser. The landing page shows:

- a one-line **headline strip** — "N genes assessed · M need review · …",
- the **⚠ Needs review** group (discordant / NO_CALL / partial) first, then **Normal — concordant**,
- click any gene card for the full **drill-down**: every tool's diplotype side by side,
  the supporting variants, and the CPIC dosing recommendations.

It is a single self-contained file — no server, no external assets. Email it, archive
it, or print it to PDF (severity colour is preserved in print).

## Verify it worked

A quick sanity check on the reference sample (GIAB HG002) is CYP2D6 → `*2/*4`,
Activity Score 1.0, Intermediate Metabolizer, 4/4 tools concordant.

For one gene only (faster, good for a smoke run), pass a `genes=` subset to the
container directly:

```bash
docker run --privileged --rm \
  -v "$(pwd)/StellarPGx:/pgx/stellarpgx" -v "$(pwd)/StellarPGx/containers:/pgx/containers" \
  -v "$(pwd)/pypgx/pypgx-bundle:/pgx/bundle" -v "/path/to/GRCh38:/pgx/ref" \
  -v "/path/to/data:/pgx/data" -v "$(pwd)/results/one_gene:/pgx/results" \
  pgx-suite:latest \
  snakemake -s /opt/pgx/Snakefile --cores 4 \
    --config bam=/pgx/data/sample.bam outdir=/pgx/results genes="CYP2D6"
```

## Troubleshooting

| Symptom | Fix |
|---------|-----|
| `FATAL: could not open image` | Mount the SIF: `-v $(pwd)/StellarPGx/containers:/pgx/containers` |
| `FATAL: kernel too old` | Run with `--privileged` |
| A gene reports `NO_CALL` | Mean depth over that gene is below the floor (`--config min_depth=10`). Real result, not a bug — that region lacks coverage to call. |
| A gene reports `DISCORDANT` | The callers disagree; open the drill-down to see each tool's call. The pipeline refuses to pick a silent winner. |
| Nextflow hangs on first run | StellarPGx's Nextflow needs outbound internet on first launch. |

## What you built

A reproducible, multi-caller PGx report for a whole-genome sample: 31 genes, six
callers, one verdict each, with a coverage gate that refuses to guess and a report
a clinician can scan in seconds or drill into for review. Next:

- [How to add or change a gene](howto-add-a-gene.md)
- [Per-tool reference](ToolsDocumentation.md) · [per-gene clinical reference](PGxDocumentation.md)
