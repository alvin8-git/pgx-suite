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
- [x] CPIC coverage: 17/19 Level A genes (up from 15/19); gaps closed in Phase 7

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

## Phase 6: Report Quality & Coverage Flags ✅ COMPLETE

- [x] **Depth coverage colour flags** on gene cards and gene detail panels
  - Thresholds based on CAP/AMP 2021 PGx guidance and CPIC tool recommendations
  - ≥80% at ≥30× → green (`depth-ok`) — reliable calling
  - 50–79%        → amber (`depth-caution`) — SNV OK, SV uncertain
  - 10–49%        → orange (`depth-poor`) — results may be unreliable
  - <10%          → red (`depth-critical`) — insufficient for calling
  - Applied in: gene card `.gene-depth` text, detail panel header, per-gene depth table
- [x] `_depth_css()` helper added to `pgx-report.py`

---

## Investigation Required — Low ≥30× Coverage Genes

Three genes consistently show near-zero ≥30× coverage in test samples (HG003, T7_NA24385).
These are flagged red in the HTML report. Root causes:

### GSTT1 — BED coordinates point to wrong region
- The GSTT1 gene body is on `chr22_KI270879v1_alt` (alternate contig), NOT on the primary
  chr22 assembly. Standard WGS alignment pipelines (BWA-MEM2 with GRCh38 primary assembly)
  do not align reads to alternate contigs.
- **Our mosdepth BED uses `chr22` coordinates** from the primary assembly where GSTT1 is
  represented only as a gap/placeholder. This region has essentially zero real read coverage.
- Fix options:
  - [ ] Use alt-aware alignment with GRCh38 full + alt contigs; add `chr22_KI270879v1_alt`
        BED entry to pgx-bamstats.sh for true GSTT1 depth
  - [ ] Or accept that GSTT1 depth is unmeasurable on primary assembly BAMs and display
        "Alt contig — depth not measurable" in the report (already done for the null entry;
        mosdepth BED entry should be removed to avoid false red flag)

### GSTM1 — Whole-gene deletion (null allele)
- GSTM1 has a common whole-gene deletion polymorphism (*0 null allele). Population
  frequency of homozygous deletion (GSTM1*0/*0) is ~40–50% in European populations.
- If the sample is GSTM1*0/*0 (null/null), there are literally no reads from this locus;
  mean depth ≈ 0× and ≥30× ≈ 0%.
- This is NOT a QC failure — it is the true biological finding (homozygous deletion).
- Fix options:
  - [ ] Annotate the depth card/table: if GSTM1 diplotype is *0/*0 or *0/*, display
        "Deletion detected — depth reflects null allele" instead of a red flag
  - [ ] Suppress depth flag for GSTM1 when SV caller reports deletion

### G6PD — X-linked hemizygous coverage in males
- G6PD is on chrX. For XY (male) samples, coverage of chrX is ~half that of autosomes
  (hemizygous, one copy vs two copies).
- At typical WGS depth of ~30× autosomal, chrX is ~15× → ≥30× fraction near 0%.
- This is expected and NOT a QC failure for male samples.
- Fix options:
  - [ ] Use inferred sex from `bam_stats.json` (`inferred_sex` field) to skip or adjust
        the ≥30× threshold for X-linked genes in male samples
  - [ ] Display "X-linked — hemizygous coverage expected (~15X in males)" annotation
        when inferred sex = MALE and gene = G6PD

### Recommended fixes (priority order)
1. [ ] Remove GSTT1 from mosdepth primary BED; show null entry as before (already partially done)
2. [ ] Use `inferred_sex` from `bam_stats.json` to annotate G6PD depth appropriately
3. [ ] Cross-reference GSTM1 diplotype call to explain zero depth as deletion finding

---

## Phase 7: CPIC Level A Gap Closure — CACNA1S and MT-RNR1 ✅ COMPLETE

Pipeline now covers **19/19 CPIC Level A genes** across **31 genes total** with **6 callers**.

---

### Phase 7a: CACNA1S ✅ COMPLETE

**Gene:** `CACNA1S` · chr1:201,006,956-201,083,927 (GRCh38) · ~77 kb · CPIC Level A
**Drug:** Volatile anaesthetics + succinylcholine (malignant hyperthermia susceptibility)
**Callers:** PyPGx + StellarPGx (`1 0 0 1`)

- [x] `pgx-run.sh`: GENE_COORDS + GENE_SUPPORT (`1 0 0 1`)
- [x] `pgx-all-genes.sh`: added to GENES list
- [x] `pgx-bamstats.sh`: BED entry `chr1 201006955 201083927 CACNA1S`
- [x] `pgx-compare.py`: GENE_SUPPORT entry (pypgx=True, stellarpgx=True)
- [x] `pgx-report.py`: GENE_LOCI + CPIC_DB (volatile anaesthetics, succinylcholine, Level A)
- [ ] Validate on HG002 (expected: MHS-normal / Reference/Reference — no MHS variants in NA24385)

---

### Phase 7b: MT-RNR1 ✅ COMPLETE

**Gene:** `MT-RNR1` (12S rRNA) · chrM:648-1601 · ~954 bp · CPIC Level A
**Drug:** Aminoglycosides — m.1555A>G (~1/500 Europeans) / m.1494C>T → irreversible hearing loss
**Caller:** mutserve v2.0.0 (standalone JAR; AF-based; Java 21 already in image)

- [x] `docker/pgx-mt.sh`: chrM read extraction → `mutserve.jar call` → JSON result
  - `--level 0.01` (1% AF detection threshold; CPIC carrier threshold ≥2%)
  - Output: `mt-rnr1/<sample>_mtrna1_result.json` with diplotype, AF, classification, notes
- [x] `Dockerfile`: mutserve.jar downloaded + pgx-mt.sh installed
- [x] `pgx-run.sh`: DO_MUTSERVE flag, `run_mt()` function, MT-RNR1 special-case bypass
- [x] `pgx-all-genes.sh`: added to GENES list
- [x] `pgx-bamstats.sh`: BED entry `chrM 647 1601 MT-RNR1`
- [x] `pgx-compare.py`: GENE_SUPPORT (mutserve=True) + `parse_mutserve()` reading JSON
- [x] `pgx-report.py`: GENE_LOCI + CPIC_DB (aminoglycosides, Level A, min_tools=1)
- [x] `ToolsDocumentation.md`: Section 6 — mutserve overview, workflow, limitations
- [ ] Validate on HG002 (expected: Reference — no m.1555A>G or m.1494C>T in NA24385)

---

## Phase 8: Documentation & Validation

- [ ] `PGxDocumentation.md` — scientific gene summary for all 31 genes
  - Overview, pharmacogenomic significance, key variants, population frequencies, CPIC level, tool support
  - Special notes: IFNL3/IFNL4 causal variant distinction (ss469415590 in IFNL4), GSTT1 alt-contig caveat,
    GSTM1 whole-gene deletion, G6PD X-linked hemizygous males, MT-RNR1 AF-based heteroplasmy model
- [ ] End-to-end validation of CACNA1S and MT-RNR1 on HG002/NA24385 and a second sample
- [ ] `test.sh` smoke tests: add mutserve JAR version check + CACNA1S support entry check

---

## Phase 9: Hypersensitivity Panel Expansion

The quick-win annotations below leverage OptiType's existing Class I genotyping output.
Class II alleles require a dedicated tool (HLA-HD or arcasHLA) and are not yet implemented.

### 9a. Annotation-only additions (OptiType already calls these — DONE ✅)

- [x] `pgx-report.py` CPIC_DB — HLA-B entry: add **HLA-B\*13:01** (dapsone → DRESS/HSS, CPIC Level B, Asian ancestry)
- [x] `pgx-report.py` CPIC_DB — HLA-B entry: add **HLA-B\*57:03** (flucloxacillin → DILI, CPIC Level B, European ancestry)

### 9b. New tool required — HLA Class II typing (HLA-HD or arcasHLA)

OptiType covers only **MHC Class I** (HLA-A, -B, -C). Class II genes need a separate tool.
Candidate tools: **HLA-HD** (IMGT-based, high resolution) or **arcasHLA** (RNA-seq + WGS).

| Gene | Allele | Drug | Reaction | Evidence |
|------|--------|------|----------|----------|
| **HLA-DRB1** | *07:01 | Abacavir (additional locus), lapatinib | DILI | PharmGKB Level 2A |
| **HLA-DQA1** | *02:01 | Lapatinib | DILI | CPIC annotation |
| **HLA-DQB1** | *06:02 | Clozapine (agranulocytosis protective) | — | Observational |

**Implementation tasks:**
- [ ] Evaluate HLA-HD vs arcasHLA for WGS (non-RNA) input compatibility with GRCh38
- [ ] Add Class II tool to `Dockerfile` (Apptainer SIF or conda install)
- [ ] Write `docker/pgx-hla2.sh`: MHC read extraction (chr6:28510020-33480577) → HLA Class II calling
- [ ] Parse Class II results in `pgx-compare.py` and `pgx-report.py`
- [ ] Add HLA-DRB1, HLA-DQA1, HLA-DQB1 to gene support matrix and CPIC_DB

### 9c. CYB5R3 (dapsone methaemoglobinaemia)

CYB5R3 encodes NADH-cytochrome b5 reductase; deficiency causes Type I congenital
methaemoglobinaemia and exacerbates dapsone-induced methaemoglobinaemia.

| Variant | Clinical impact | CPIC Level |
|---------|----------------|------------|
| *2 (rs1799931 → G72E) | Reduced enzyme activity | B |

- [ ] Assess gene size (chr22:42,511,583-42,556,659, ~45 kb) — bcftools lookup is feasible
- [ ] Determine PyPGx / Aldy support (likely no star-allele support; custom bcftools genotyper)
- [ ] Implement variant-level lookup for key CYB5R3 variants (rs1799931)
- [ ] Add to CPIC_DB with cross-reference to HLA-B\*13:01 + dapsone context

### 9d. HLA-C clinical reporting (OptiType already genotypes internally)

OptiType outputs HLA-C alleles but `pgx-report.py` currently ignores them.

| Allele | Drug | Reaction | Evidence |
|--------|------|----------|----------|
| C*04:01 | Clozapine | Agranulocytosis | Observational (CPIC annotation) |
| C*08:01 | Abacavir (additive) | Hypersensitivity | PharmGKB |

- [ ] Add HLA-C to `GENE_LOCI` (chr6:31,268,749-31,272,092)
- [ ] Parse HLA-C allele from OptiType output in `pgx-hla.sh` + `pgx-compare.py`
- [ ] Add CPIC_DB entry for HLA-C with C*04:01 (clozapine) annotation

---

## Phase 10: Beyond CPIC Level A — Population PGx Expansion

After 19/19 CPIC Level A coverage, the next highest-yield additions for a population-based
genome project are:

### Priority 1 — Strong evidence, transporter / drug interaction scope

| Gene | CPIC Level | Clinical relevance | Notes |
|------|-----------|-------------------|-------|
| **ABCB1** (MDR1/P-glycoprotein) | B | Efflux transporter; wide drug interaction scope: tacrolimus, digoxin, antiretrovirals, statins, loperamide | chr7:87,133,500-87,342,000; PyPGx supported |
| **F5** (Factor V Leiden) | — | OCP / HRT / surgery thrombotic risk; rs1799963 G20210A | Not a classic PGx gene but population-critical for prescribing safety |
| **F2** (Prothrombin) | — | Thrombophilia; partner to F5 for VTE risk assessment | Relevant for warfarin initiation and OCP safety |

### Priority 2 — Expanded metaboliser genes

| Gene | CPIC Level | Notes |
|------|-----------|-------|
| **UGT2B15** / **UGT2B17** | B | Androgen metabolism; tamoxifen, lorazepam; copy-number deletions common in Asian populations |
| **CYP1B1** | — | Cancer pharmacogenomics; tamoxifen oestrogen metabolism; widely studied |
| **MTHFR** | — | Methotrexate, folate; controversial evidence but widely ordered clinically; C677T + A1298C |

### Priority 3 — Psychiatric / neurology focus

| Gene | Notes |
|------|-------|
| **HTR2A** | SSRI, clozapine, antipsychotic response (CPIC B) |
| **COMT** | Pain, dopamine, antipsychotic response (Val158Met, rs4680) |
| **ANKK1/DRD2** | Tardive dyskinesia risk with antipsychotics |

### Southeast Asian population priorities (Singapore / TTSH context)

| Gene | Rationale |
|------|-----------|
| HLA-B*15:02 | Already covered via HLA-B ✅ |
| HLA-B*58:01 | Already covered via HLA-B ✅ |
| CYP2C19 poor metabolizer | High frequency in East Asian populations ✅ |
| G6PD A- / Mediterranean | Already covered ✅ |
| CYP2B6 *6 | ~50% allele frequency in African populations → efavirenz/HIV ✅ |
| **SLCO1B1 *15** | Higher frequency in Asian populations → already covered ✅ |
| **UGT2B17 deletion** | ~70% deletion frequency in East Asians; tamoxifen, steroid metabolism |

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
