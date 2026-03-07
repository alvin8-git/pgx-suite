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
- [ ] Dockerfile updated to install pgx-run.sh + pgx-compare.py → rebuild image needed

### 2e. Integration test with real BAM

- [ ] Obtain GeT-RM NA12878 CYP2D6 BAM (chr22 region) for end-to-end test
  - [ ] Alternatively: user-supplied BAM from clinical/research pipeline
- [ ] Run `pgx-run.sh CYP2D6 /pgx/data/NA12878.bam` end to end
- [ ] Verify all 4 tools return concordant diplotype
- [ ] Benchmark runtime for single-gene call

---

## Phase 3: Multi-gene batch mode (future)

- [ ] `pgx-run.sh --all-genes <BAM>` — run the full supported gene panel
- [ ] Parallel execution per gene (GNU parallel or Nextflow wrapper)
- [ ] HTML report generation

---

## Known Limitations / Notes

- `StellarPGx` requires `--privileged` Docker flag (Apptainer inside Docker)
- `pypgx-bundle` (500 MB, 1KGP VCFs) is volume-mounted, not baked into image
- Stargazer and Aldy are **non-commercial academic use only** — do not publish image to public registries
- Stargazer's GDF-based CNV calling requires hg38-aligned depth data (not yet integrated)
- StellarPGx `nextflow.config` `$PWD` references assume pipeline is run from the StellarPGx repo directory — the volume mount at `/pgx/stellarpgx` handles this
