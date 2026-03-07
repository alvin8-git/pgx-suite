# PGx Suite — Changelog

All notable changes to this project are documented here.
Format: reverse-chronological, grouped by phase/milestone.

---

## 2026-03-07 — Phase 2: Orchestration & SV-Aware Calling

### New files

#### `docker/pgx-run.sh`
Main orchestration entry point. Accepts `<GENE> <BAM> [--ref ...] [--output ...]`
and runs all supported tools for the requested gene, then calls `pgx-compare.py`.

Key behaviours:
- Validates BAM + `.bai` index, GRCh38 reference FASTA + `.fai`, and gene support
  before doing any work.
- Extracts sample name from the BAM `@RG SM:` tag; falls back to filename stem.
- Generates a shared gene-region VCF via `bcftools mpileup | bcftools call | tabix`
  for tools that need it (PyPGx, Stargazer). Region bounds taken from the pypgx
  `gene-table.csv` (includes upstream/downstream buffer).
- Runs each applicable tool sequentially; individual failures are logged and the
  script continues to the next tool rather than aborting.
- Calls `pgx-compare.py` at the end to produce the comparison table and TSV.

**SV handling — PyPGx** (`CYP2A6 CYP2B6 CYP2D6 CYP2E1 CYP4F2 G6PD GSTM1 GSTT1`):
These eight genes have copy number variation (SVs) that cannot be detected from SNV
calls alone. For them, two pre-processing steps run before `run-ngs-pipeline`:

1. `pypgx prepare-depth-of-coverage <output>/depth-of-coverage.zip <bam> --assembly GRCh38`
   — computes per-base read depth for all SV target gene loci.
2. `pypgx compute-control-statistics VDR <output>/control-stats-VDR.zip <bam> --assembly GRCh38`
   — computes normalisation statistics for the VDR control locus.

`run-ngs-pipeline` is then called with `--depth-of-coverage` and
`--control-statistics` so PyPGx can model copy number and call SV alleles
(e.g. CYP2D6\*5 deletion, \*2x2 duplication).

Source: pypgx `gene-table.csv` (`SV=TRUE` column); VDR confirmed as the standard
control gene in the pypgx GeT-RM WGS tutorial.

**SV handling — Stargazer** (`CYP2A6 CYP2B6 CYP2D6` — genes with pseudogene paralogs):
These three genes have pseudogene paralogs (CYP2A7, CYP2B7, CYP2D7) whose read
depth must be normalised against a control locus to call copy number. Without this,
Stargazer runs in VCF-only mode and cannot detect deletions or duplications.

For these genes a GDF depth-profile file is created first:
```
stargazer -G <gene>.gdf -t <gene> -c vdr -B <bam> -o <gdf_dir> -a grc38 -d wgs
```
The main Stargazer call then adds `-c vdr -g <gdf_file>` to enable SV detection.
If GDF creation fails, the script falls back to VCF-only mode and logs a warning.

Source: Stargazer `gene_table.tsv` (`paralog` column); CLI `--help` confirms the
three supported control genes (`egfr`, `ryr1`, `vdr`).

**SV handling — Aldy**: No special parameters. Aldy's integer linear programming
solver detects copy number configurations automatically from the BAM coverage signal.
Supports CYP2D6, CYP2B6, CYP4F2, and others natively. Minimum 40× coverage
recommended for reliable CN calls.

Source: Aldy documentation — "Copy number and structural variation supported";
CN detection confirmed automatic with no extra CLI flags.

**SV handling — StellarPGx**: No special parameters. StellarPGx uses graphtyper
internally, which performs SV-aware variant calling as part of its standard pipeline.

**GSTT1 special case**: The GRCh38 locus for GSTT1 resides on an alternate contig
(`chr22_KI270879v1_alt`). bcftools mpileup is skipped for this gene; Stargazer is
also inapplicable (no Stargazer support). Only StellarPGx runs for GSTT1.

**Gene coordinate table**: Updated from manually-curated approximate coordinates to
the authoritative regions from `pypgx/pypgx/api/data/gene-table.csv` (`GRCh38Region`
column). These wider windows (including upstream/downstream buffer) ensure all
relevant variants are captured:

| Gene | Old coords (approx) | New coords (pypgx) |
|------|--------------------|--------------------|
| CYP2D6 | chr22:42126499-42130806 | chr22:42116498-42155810 |
| CYP2B6 | chr19:40991352-41018693 | chr19:40921281-41028398 |
| CYP2C9 | chr10:94938683-94994564 | chr10:94935657-94993091 |
| CYP3A4 | chr7:99756967-99800931  | chr7:99753966-99787184  |
| CYP4F2 | chr19:15966804-16009734 | chr19:15863022-15913074 |
| SLCO1B1 | chr12:21282608-21396854 | chr12:21128193-21242796 |
| *(all others similarly updated)* | | |

---

#### `docker/pgx-compare.py`
Parses each tool's output directory and produces a side-by-side comparison table
plus a `<GENE>_<SAMPLE>_comparison.tsv` file.

Parsers implemented:
- **PyPGx**: opens `pypgx/results.zip`, extracts `data.tsv`, reads `Diplotype`,
  `ActivityScore`, `Phenotype` columns.
- **Stargazer**: reads `stargazer/genotype.txt` (handles `#`-prefixed header line),
  extracts `Genotype`, `Activity_score`, `Phenotype`.
- **Aldy**: reads `aldy/<GENE>.aldy`; tries named columns from the tab-separated
  header first, then falls back to a heuristic star-allele pattern scan.
- **StellarPGx**: scans `stellarpgx/star_allele_calls/` for files containing
  `*N/*N` diplotype patterns.

Each parser returns a `CallerResult` with `status` ∈ `{not_run, ok, failed}`.
Tools that were not run for the gene (per the support matrix) are silently omitted.

Output table example:
```
========================================================================
 PGx Star Allele Results
 Gene:   CYP2D6
 Sample: NA12878
 Build:  GRCh38
 SV:     SV mode — PyPGx: depth-of-coverage + VDR control stats; Stargazer: GDF from BAM (paralog CN normalisation)
========================================================================
Tool          Diplotype         Activity Score    Phenotype
────────────────────────────────────────────────────────────────────────
PyPGx         *1/*4             1.0               Normal Metabolizer
Stargazer     *1/*4             1.0               Normal Metabolizer
Aldy          *2/*4             1.0               Normal Metabolizer
StellarPGx    *1/*4             -                 -
────────────────────────────────────────────────────────────────────────
Concordance: 3/4 tools agree on *1/*4
========================================================================

Full results saved to: /pgx/results/CYP2D6_NA12878_comparison.tsv
```

Added in this session (SV annotation update):
- `PYPGX_SV_GENES` and `STARGAZER_SV_GENES` frozensets for annotation.
- `SV:` line in the printed header shows which SV mode was used per tool.
- TSV output includes a `SVMode` column.

---

### Modified files

#### `Dockerfile`
Added Phase 2 scripts to the image (step 12):
```dockerfile
COPY docker/pgx-run.sh     /opt/pgx/pgx-run.sh
COPY docker/pgx-compare.py /opt/pgx/pgx-compare.py
RUN chmod +x /opt/pgx/test.sh /opt/pgx/pgx-run.sh \
    && ln -s /opt/pgx/pgx-run.sh /usr/local/bin/pgx-run.sh
```
`pgx-run.sh` is symlinked onto `$PATH` so it can be called as `pgx-run.sh CYP2D6 /pgx/data/sample.bam`
without specifying the full path. Image requires a rebuild to pick up these changes.

#### `TODO.md`
Phase 2 tasks marked complete (2a–2d). Phase 2e (integration test with real BAM)
and Phase 3 (multi-gene batch mode) remain open.

---

## 2026-03-06 — Phase 1: Docker Container & Smoke Tests ✅

Initial commit (`077dd3c`). Built and validated a single `pgx-suite:latest`
Docker image containing all four PGx callers for GRCh38.

### Files created

| File | Description |
|------|-------------|
| `Dockerfile` | Ubuntu 22.04 + Python 3.11 + Java 21 + Apptainer + Nextflow |
| `.dockerignore` | Excludes large data volumes from build context |
| `docker/test.sh` | Smoke tests for all 4 tools (no BAM required) |
| `docker/docker-run.sh` | Convenience `docker run` wrapper with volume mounts |
| `docker/README.md` | Container-specific usage documentation |
| `docs/plans/2026-03-06-pgx-docker-design.md` | Architecture design document |
| `docs/plans/2026-03-06-pgx-docker-implementation.md` | Step-by-step implementation plan |
| `README.md` | Top-level project README |
| `TODO.md` | Phase tracking document |

### Key design decisions recorded in Phase 1

- **Apptainer before Python 3.11**: the Apptainer PPA step must run while the system
  Python is still 3.10, because `add-apt-repository` uses `python3-apt` compiled
  against 3.10.
- **Stargazer wrapper script**: Stargazer uses Python 2-style bare imports that only
  resolve when `__main__.py` is the script entry point (not when installed via pip).
  A `/usr/local/bin/stargazer` wrapper calls `python3 /opt/stargazer/stargazer/__main__.py`.
- **Java 21**: Nextflow 25.x requires Java 21; the original design specified Java 11.
- **`--privileged` for Apptainer**: running Singularity/Apptainer inside Docker
  requires kernel namespace access; `--privileged` is the supported approach for
  local workstations.
- **Volume mounts at runtime**: `pypgx-bundle` (500 MB VCFs) and the StellarPGx
  pipeline + SIF container are not baked into the image to keep the image size
  manageable and to allow updates without rebuilding.

### Smoke test results (all passing)

| Tool | Test | Result |
|------|------|--------|
| PyPGx 0.26.0 | `pypgx --version` + import + `run-ngs-pipeline CYP4F2` | PASS |
| Stargazer 2.0.3 | `stargazer --version` + WGS example (CYP2D6/VDR GDF) | PASS |
| Aldy 4.8.3 | `aldy --version` + `aldy test` (built-in test suite) | PASS |
| StellarPGx 1.2.7 | `nextflow -version` + `apptainer --version` | PASS |

---

## Known limitations

- **GSTT1**: Located on `chr22_KI270879v1_alt` (alternate contig in GRCh38).
  bcftools mpileup and Stargazer are not applicable. Only StellarPGx runs for GSTT1.
- **StellarPGx requires `--privileged`**: Apptainer inside Docker needs `SYS_ADMIN`.
- **pypgx-bundle not in image**: 1KGP VCFs (~500 MB) are volume-mounted at runtime
  (`-v ./pypgx/pypgx-bundle:/pgx/bundle`). Required for haplotype phasing.
- **Stargazer / Aldy licensing**: Non-commercial academic use only. This image must
  not be pushed to any public registry (Docker Hub, GHCR, etc.).
- **Phase 2e pending**: End-to-end test with a real BAM (e.g. GeT-RM NA12878 chr22
  region) has not yet been run. Awaiting a test BAM file.
