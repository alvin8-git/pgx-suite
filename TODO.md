# PGx Suite — TODO

## Phase 1: Docker Container & Smoke Tests ✅ COMPLETE

- [x] Design document (`docs/plans/2026-03-06-pgx-docker-design.md`)
- [x] Implementation plan (`docs/plans/2026-03-06-pgx-docker-implementation.md`)
- [x] `Dockerfile` — Ubuntu 22.04, Python 3.11, Java 21, Apptainer, Nextflow
- [x] `.dockerignore` — exclude large data volumes from build context
- [x] Install PyPGx 0.26.0 (from source)
- [x] Install Aldy 4.8.3 (via pip)
- [x] Install Stargazer 2.0.3 GRCh38 build (wrapper script, bare-import compatible)
- [x] Install StellarPGx 1.2.7 via Nextflow + Apptainer inside container
- [x] `docker/test.sh` — bundled smoke tests for all 4 tools (no BAM required)
- [x] `docker/docker-run.sh` — convenience run wrapper with correct volume mounts
- [x] `docker/README.md` — container-specific documentation

---

## Phase 2: `pgx-run.sh` Orchestrator (awaiting test BAM file)

### 2a. Input validation script

- [x] `docker/pgx-run.sh` — main entry point
  - [x] Accept `<GENE> <BAM_FILE> [--ref /pgx/ref/hg38.fa] [--output /pgx/results]`
  - [x] Validate: BAM file exists and is indexed (`.bai` present)
  - [x] Validate: gene is supported by at least one tool
  - [x] Check reference FASTA is mounted at `/pgx/ref/`

### 2b. Shared preprocessing

- [x] bcftools mpileup → gene-region VCF (shared input for PyPGx and Stargazer)
  - [x] Extract gene chromosomal coordinates from a lookup table (28 genes in pgx-run.sh)
  - [x] `bcftools mpileup -r chr{N}:{start}-{end} -f /pgx/ref/hg38.fa <bam>`
  - [x] `bcftools call -m -v -o results/<gene>.vcf.gz`
  - [x] `tabix -p vcf results/<gene>.vcf.gz`

### 2c. Per-tool callers

- [x] PyPGx caller block
- [x] Stargazer caller block (with per-gene control gene mapping)
- [x] Aldy caller block
- [x] StellarPGx caller block

### 2d. Result aggregation

- [x] `docker/pgx-compare.py` — parse all 4 tool outputs
  - [x] Gene support matrix (27 genes)
  - [x] Extract: diplotype, activity score, phenotype per tool
  - [x] Calculate concordance (how many tools agree on top diplotype)
  - [x] Save `results/<GENE>_<SAMPLE>_comparison.tsv`
- [x] Dockerfile updated to install pgx-run.sh + pgx-compare.py + pgx-all-genes.sh + pgx-bamstats.sh + pgx-report.py → rebuild image needed

### 2e. Integration test with real BAM ✅ COMPLETE

- [x] Test BAM: `T7_NA24385.bwa.sortdup.bqsr.bam` (HG002/NA24385, 167 GB WGS, GRCh38)
- [x] Single-gene: `pgx-run.sh CYP2D6` — 4/4 concordant (`*2/*4`, AS 1.0, Intermediate Metabolizer)
- [x] All-gene batch: `pgx-all-genes.sh` — 26/26 genes completed; HTML report generated
- [x] Wall time for CYP2D6: ~53 sec (bottleneck: PyPGx sklearn CNV caller ~49 sec)

---

## Phase 3: Multi-gene batch mode & HTML reports ✅ COMPLETE

- [x] `docker/pgx-all-genes.sh` — run all 27 supported genes for a single BAM
  - [x] Parallel execution with a configurable job pool (`--jobs N`, default 4)
  - [x] Per-gene stdout/stderr captured to `<output>/logs/<GENE>.log`
  - [x] Aggregate `all_genes_summary.tsv` written at the end
  - [x] BAM QC launched in background (`pgx-bamstats.sh &`) to overlap with gene calling
- [x] `docker/pgx-bamstats.sh` — BAM QC pipeline → `bam_stats.json`
  - [x] `samtools idxstats` (28–46 ms) — total/mapped reads from `.bai` index; no BAM scan
  - [x] `samtools stats` on chr1 4 Mb window (~2 sec) — read length, insert size, error rate
  - [x] Genome depth + sex inference computed from idxstats × read_length (no extra scan)
  - [x] `mosdepth v0.3.12` with `--by pgx.bed --thresholds 0,20,30 --fast-mode` (~550 sec) — per-gene mean depth + ≥20×/≥30× fractions for 14 PGx genes
  - [x] `samtools view -c` ×3 on chr1 4 Mb window (~500 ms) — MAPQ≥20% and duplicate% estimate
  - [x] Timing instrumentation via `date +%s%3N` on every command
- [x] `docker/pgx-compare.py` extended to emit `<GENE>_<SAMPLE>_detail.json` (17-field × 4-tool rich output)
- [x] `docker/pgx-report.py` — HTML report generation
  - [x] Landing page: 26-gene card grid colour-coded by tool concordance
  - [x] Per-gene depth displayed on gene cards (mean depth × and ≥20× fraction) from `bam_stats.json`
  - [x] Per-gene detail page: 17×4 table (diplotype, variants, phenotype, SV, phasing, …)
  - [x] BAM QC panel when `bam_stats.json` present
  - [x] HTML written to `<output>/html_reports/` via `--html-dir` flag
- [x] `Dockerfile` — added mosdepth v0.3.12 and sambamba v1.0.1
- [x] Validated end-to-end on `T7_NA24385` (HG002/NA24385, 167 GB WGS): 26/26 genes OK; CYP2D6 4/4 concordant

---

---

## Phase 4: ABCG2 + HLA typing ✅ COMPLETE

- [x] ABCG2 (rosuvastatin, CPIC Level A) added to pipeline:
  - [x] Aldy + StellarPGx callers (PyPGx and Stargazer have no ABCG2 support)
  - [x] pgx-run.sh: GENE_COORDS + GENE_SUPPORT (0 0 1 1)
  - [x] pgx-all-genes.sh: GENES list
  - [x] pgx-bamstats.sh: BED (chr4:88085264-88236626)
  - [x] pgx-compare.py: GENE_SUPPORT + `optitype` key placeholder
  - [x] pgx-report.py: CPIC_DB entry (diplotype_check on rs2231142, rosuvastatin Level A)
- [x] HLA-A / HLA-B (CPIC Level A) added via OptiType:
  - [x] `docker/pgx-hla.sh` — new script: samtools MHC extraction → OptiType Apptainer SIF
  - [x] pgx-run.sh: DO_OPTITYPE flag, run_hla() function, HLA special-case bypass
  - [x] pgx-all-genes.sh: HLA-A + HLA-B in GENES list
  - [x] pgx-compare.py: parse_optitype() reading A1/A2 or B1/B2 from result TSV
  - [x] pgx-report.py: TOOLS now includes "OptiType"; HLA-A/HLA-B CPIC_DB with risk_alleles, min_tools=1
  - [x] Dockerfile: pgx-hla.sh installed + symlinked to /usr/local/bin/
  - [x] OptiType SIF pulled separately: `apptainer pull --name optitype.sif docker://quay.io/biocontainers/optitype:1.3.5--hdfd78af_1`
- [x] Key Clinical Findings cards now clickable → link to gene detail HTML page
- [x] CPIC coverage: 17/19 Level A genes (up from 15/19); gaps: CACNA1S (natural MH partner to RYR1), MT-RNR1 (mitochondrial)

---

## Phase 5: HTML Report Enhancements ✅ COMPLETE

- [x] CPIC Clinical Reference section on **all 29 gene detail pages** (previously absent for CYP1A1, CYP1A2, CYP2A6, CYP2C8, CYP2E1, CYP3A4, GSTM1, GSTT1, IFNL3, NAT1, POR)
  - [x] Drug–gene pair table with CPIC Level badge and recommendation text
  - [x] Patient-specific finding note derived from consensus diplotype / phenotype
  - [x] PharmVar + cpicpgx.org links
- [x] **In-page navigation bar** at top of every gene detail page → anchored links to Tool Results, Variant Evidence, CPIC Reference sections
- [x] **Chromosomal locus** (e.g. `chr22:42,116,498-42,155,810 (GRCh38)`) displayed under gene name on every detail page — monospace, from `GENE_LOCI` dict
- [x] **HLA-A and HLA-B depth** added to gene cards and detail pages:
  - [x] `pgx-bamstats.sh`: HLA-A (`chr6:29,905,000-29,920,000`) and HLA-B (`chr6:31,316,000-31,330,000`) added to mosdepth BED — now 27 primary-assembly genes
  - [x] `bam_stats.json` `gene_depth` populated for HLA-A and HLA-B
  - [x] Gene cards show `Depth: {mean}X | ≥30X: {pct}%` for all 29 genes
- [x] **OptiType parse_optitype() double-prefix bug fixed** — `HLA-A*A*01:01` → `HLA-A*01:01`; `HLA-B*B*57:01` → `HLA-B*57:01`
- [x] Validated end-to-end on T7_NA24385 (29/29 genes) and HG003 (29/29 genes)

---

## Known Limitations / Notes

- `StellarPGx` requires `--privileged` Docker flag (Apptainer inside Docker)
- `pypgx-bundle` (500 MB, 1KGP VCFs) is volume-mounted, not baked into image
- Stargazer and Aldy are **non-commercial academic use only** — do not publish image to public registries
- Stargazer's GDF-based CNV calling requires hg38-aligned depth data (not yet integrated)
- StellarPGx `nextflow.config` `$PWD` references assume pipeline is run from the StellarPGx repo directory — the volume mount at `/pgx/stellarpgx` handles this
- mosdepth is I/O bound at ~490 MB/s on HDD; ~550 sec on 167 GB WGS is a hardware ceiling — no threading or tool choice overcomes it. Running bamstats in background (overlapping gene calling) is the practical mitigation.
- sambamba v1.0.1 is installed in the image but not used by default (benchmarked at same speed as samtools flagstat — both I/O bound, not CPU bound on full BAM scan)
- Future speedup option: pre-generate a persistent PGx-region subset BAM upstream; mosdepth on the subset takes ~4 sec vs ~550 sec on the full BAM. Not implemented — the subset extraction itself takes ~560 sec (BAM re-compression overhead).
