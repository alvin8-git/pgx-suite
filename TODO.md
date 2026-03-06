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

- [ ] `docker/pgx-run.sh` — main entry point
  - [ ] Accept `<GENE> <BAM_FILE> [--ref /pgx/ref/hg38.fa] [--output /pgx/results]`
  - [ ] Validate: BAM file exists and is indexed (`.bai` present)
  - [ ] Validate: gene is supported by at least one tool
  - [ ] Check reference FASTA is mounted at `/pgx/ref/`

### 2b. Shared preprocessing

- [ ] bcftools mpileup → gene-region VCF (shared input for PyPGx and Stargazer)
  - [ ] Extract gene chromosomal coordinates from a lookup table
  - [ ] `bcftools mpileup -r chr{N}:{start}-{end} -f /pgx/ref/hg38.fa <bam>`
  - [ ] `bcftools call -m -v -o results/<gene>.vcf.gz`
  - [ ] `tabix -p vcf results/<gene>.vcf.gz`

### 2c. Per-tool callers

- [ ] PyPGx caller block
  - [ ] `pypgx run-ngs-pipeline --gene <GENE> --assembly GRCh38`
  - [ ] `--variants results/<gene>.vcf.gz --output results/pypgx/`
  - [ ] Parse output: `results/pypgx/*/results.zip` → `data.tsv`

- [ ] Stargazer caller block (genes: CYP2D6, CYP2C19, CYP2C9, CYP2B6, CYP2C8, CYP3A4/5, CYP4F2, SLCO1B1, G6PD, VKORC1, CYP1A1/2, CYP2A6, CYP2E1)
  - [ ] `stargazer -t <gene> -d wgs -a grc38 -i results/<gene>.vcf.gz -o results/stargazer/`
  - [ ] Optionally add `-c <control_gene> -g <gdf_file>` for CNV analysis when GDF available
  - [ ] Parse output: `results/stargazer/<gene>.genotype.txt`

- [ ] Aldy caller block
  - [ ] `aldy genotype -g <GENE> -p hg38 <bam> -o results/aldy/`
  - [ ] Parse output: `results/aldy/<gene>.aldy`

- [ ] StellarPGx caller block (genes: CYP2D6, CYP2C19, CYP2C9, CYP2B6, CYP2C8, CYP3A4/5, CYP4F2, SLCO1B1, NUDT15, TPMT, UGT1A1, NAT1, NAT2, GSTM1, GSTT1, POR/CYPOR, CYP1A1/2, CYP2A6, CYP2E1)
  - [ ] `nextflow run /pgx/stellarpgx/main.nf --gene <gene> --in_bam "<bam>*{bam,bai}" --ref_file /pgx/ref/hg38.fa --out_dir results/stellarpgx/`
  - [ ] Parse output: `results/stellarpgx/star_allele_calls/`

### 2d. Result aggregation

- [ ] `docker/pgx-compare.py` — parse all 4 tool outputs
  - [ ] Gene support matrix (which tools run for each gene)
  - [ ] Extract: diplotype, activity score, phenotype per tool
  - [ ] Calculate concordance (how many tools agree on top diplotype)
  - [ ] Output format:
    ```
    ===== PGx Star Allele Results =====
    Gene:   CYP2D6
    Sample: NA12878
    Build:  GRCh38
    ────────────────────────────────────
    Tool         Diplotype  Score  Phenotype
    ────────────────────────────────────
    PyPGx        *1/*4      1.0    Normal Metabolizer
    Stargazer    *1/*4      1.0    Normal Metabolizer
    Aldy         *2/*4      1.0    Normal Metabolizer
    StellarPGx   *1/*4      -      -
    ────────────────────────────────────
    Concordance: 3/4 tools agree (*1/*4)
    ```
  - [ ] Save `results/comparison.tsv`

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
