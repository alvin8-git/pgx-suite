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

## Phase 7: CPIC Level A Gap Closure — CACNA1S and MT-RNR1

Two CPIC Level A genes remain unsupported. Coverage would bring the pipeline to 19/19 Level A genes.

---

### CACNA1S — Malignant Hyperthermia Susceptibility (MHS) partner to RYR1

**Gene:** `CACNA1S` · chr1:201,006,956-201,083,927 (GRCh38) · ~77 kb · CPIC Level A
**Drug:** Volatile anaesthetics (halothane, sevoflurane, desflurane, isoflurane) and succinylcholine
**Clinical significance:** Loss-of-function variants → malignant hyperthermia susceptibility (MHS). Key actionable
variants: c.520C>T (p.Arg174Trp), c.3257G>A (p.Arg1086His). Diplotype-based recommendation: MHS-susceptible vs.
MHS-normal phenotype.

**Tool support (current):**
| Tool       | CACNA1S support | Notes |
|------------|----------------|-------|
| PyPGx      | ✅ Yes          | Supported since ≥ v0.19; SV panel available; requires `--panel wgs` |
| Aldy       | ❌ No           | Not in Aldy 4.x gene database |
| Stargazer  | ❌ No           | Not in Stargazer 2.x gene list |
| StellarPGx | ✅ Yes          | Supported in StellarPGx ≥ 1.2; part of default gene set |

**What is needed:**
- [ ] Add `CACNA1S` to `GENE_COORDS` and `GENE_SUPPORT` (`1 0 0 1`) in `pgx-run.sh`
- [ ] Add `CACNA1S` to `GENES` list in `pgx-all-genes.sh`
- [ ] Add BED entry `chr1:201006956-201083927` to `pgx-bamstats.sh` mosdepth BED
- [ ] Extend `pgx-compare.py` `GENE_SUPPORT` dict and verify PyPGx `parse_pypgx()` handles CACNA1S
  (PyPGx outputs `MHS-susceptible` / `MHS-normal` in the Phenotype column — not a standard metabolizer label)
- [ ] Add CPIC_DB entry to `pgx-report.py`:
  - drugs: volatile anaesthetics + succinylcholine
  - phenotypes: MHS-susceptible → avoid triggering agents; MHS-normal → standard care
  - risk_alleles: `*2` (Arg174Trp), `*3` (Arg1086His)
- [ ] Add `GENE_LOCI` entry for CACNA1S in `pgx-report.py`
- [ ] Validate on HG002 (expected: no MHS variants → MHS-normal)

**No new tools required** — PyPGx and StellarPGx already handle CACNA1S with standard WGS BAM input. This
is a configuration-only addition, similar to how ABCG2 was added in Phase 4.

---

### MT-RNR1 — Aminoglycoside-Induced Hearing Loss

**Gene:** `MT-RNR1` (12S rRNA) · chrM:648-1601 (GRCh38, NC_012920.1) · ~954 bp · CPIC Level A
**Drug:** Aminoglycosides (gentamicin, tobramycin, amikacin, streptomycin, neomycin)
**Clinical significance:** m.1555A>G (most common, ~1/500 Europeans) and m.1494C>T → ribosomal hypersensitivity
to aminoglycosides → irreversible sensorineural hearing loss. Even a single dose can be ototoxic in carriers.

**Tool support (current):**
| Tool       | MT-RNR1 support | Notes |
|------------|----------------|-------|
| PyPGx      | ❌ No           | No mitochondrial gene support |
| Aldy       | ❌ No           | No mitochondrial gene support |
| Stargazer  | ❌ No           | No mitochondrial gene support |
| StellarPGx | ❌ No           | No mitochondrial gene support |

**Why none of the current tools work:**
- All four tools are diplotype callers designed for nuclear (biallelic) genes.
- Mitochondrial DNA is polyploid (hundreds to thousands of copies per cell), maternally inherited,
  and haploid (no second allele). Standard genotype callers report the wrong ploidy and fail.
- Heteroplasmy: a variant may be present in 2–98% of mt copies. A threshold of ≥5% allele
  fraction is typical for clinical reporting.
- NUMTs (nuclear mitochondrial DNA segments) can create false-positive calls at chrM positions
  if not properly filtered.

**Dedicated mitochondrial variant callers required:**

| Tool | Approach | Availability | Notes |
|------|----------|--------------|-------|
| **mtDNA-Server 2 + mutserve** | Nextflow pipeline; mutserve for SNPs, Mutect2 for indels, or "fusion" mode for both | Free; Nextflow + Docker; runs via `nextflow run genepi/mtdna-server-2` | Outputs `variants.annotated.txt` (39 columns incl. AF, gene locus, pathogenicity, population freq) + haplogroup + HTML report. Detection limit 0.02 (2%) by default. NUMT annotation column from 1000GP Phase 1 data. Requires ≥50× mean chrM coverage. **Simplest integration path.** |
| **mutserve standalone** | SNP-only AF caller for chrM | Free; JAR/CLI; `--contig-name chrM` for WGS BAMs | Outputs VCF with AF field or tab-delimited text. `--level 0.01` default (1%). Handles homoplasmic + heteroplasmic. Does NOT call indels. Used as SNP engine inside mtDNA-Server 2. |
| **GATK Mutect2 `--mitochondria-mode`** | Somatic-style AF caller adapted for chrM; SNP + indel | Free, GATK 4.2+; ~1.5 GB image | GATK best-practices mtDNA workflow; strand-bias + NUMT contamination filters via `FilterMutectCalls`. More complex setup. Outputs VCF with AF field. |
| **haplocheck** | Contamination detection from mtDNA VCF | Free, CLI | Post-calling step; not a variant caller |
| **HaploGrep3 / Yleaf** | Haplogroup assignment | Free | Context only; not actionable for m.1555A>G |

**Recommended approach — mtDNA-Server 2 (fusion mode):**

mtDNA-Server 2 is the most practical integration path: it is a Nextflow pipeline requiring only
Docker (no GATK installation), handles a WGS BAM directly, and its "fusion" mode combines
mutserve (SNPs, including m.1555A>G / m.1494C>T) with Mutect2 (indels). It annotates variants
with 1000GP NUMT flags, pathogenicity predictions, and population frequencies — all in one run.

```bash
# Minimal run (single BAM, fusion mode, detection limit 2%)
nextflow run genepi/mtdna-server-2 -r v2.1.16 \
    -profile docker \
    --input  /pgx/data/${SAMPLE}.bam \
    --output /pgx/results/MT-RNR1/ \
    --mode fusion \
    --detection_limit 0.02 \
    --contig-name chrM
```

Key output: `variants.annotated.txt` — parse for positions 1555 (A→G) and 1494 (C→T).
- AF ≥ 0.02 at m.1555A>G or m.1494C>T → carrier (flag for aminoglycoside avoidance)
- No variant detected → non-carrier (standard care)
- Report AF value alongside classification (e.g. homoplasmic 0.99 vs heteroplasmic 0.15)

**What is needed:**
- [ ] New script `docker/pgx-mt.sh`:
  - Run mtDNA-Server 2 Nextflow pipeline in Docker-in-Docker (or Apptainer equivalent)
  - Parse `variants.annotated.txt` for chrM:1555 and chrM:1494
  - Output: `{gene}/{sample}_mtrna1_result.json` — variant, AF, classification (carrier/non-carrier)
  - Minimum coverage check: warn if mean chrM depth < 50× (mtDNA-Server 2 requirement)
- [ ] Add **Nextflow** to `Dockerfile` (already present — used by StellarPGx) ✅
- [ ] Add **Docker-in-Docker** capability or convert mtDNA-Server 2 to **Apptainer SIF**:
  - Preferred: `apptainer pull mtdna-server-2.sif docker://ghcr.io/genepi/mtdna-server-2:v2.1.16`
  - Then run via `apptainer exec mtdna-server-2.sif nextflow run ...` (consistent with OptiType pattern)
  - The Nextflow pipeline itself handles internal tool execution — no separate GATK install needed
- [ ] Add `MT-RNR1` to `pgx-run.sh` as a special case (like HLA-A/B calling pgx-hla.sh)
  - `run_mt()` function that calls `pgx-mt.sh` instead of standard mpileup → tool callers chain
- [ ] Add `MT-RNR1` to `pgx-all-genes.sh` GENES list
- [ ] Extend `pgx-compare.py` with `parse_mtrna1()`:
  - Read `_mtrna1_result.json`; emit diplotype-style row: e.g. `m.1555A>G(AF=0.98)` or `Reference`
  - Phenotype: `Aminoglycoside-ototoxicity risk` vs. `Standard risk`
- [ ] Add CPIC_DB entry to `pgx-report.py`:
  - drugs: gentamicin, tobramycin, amikacin, streptomycin, neomycin
  - risk_alleles: `m.1555A>G`, `m.1494C>T`
  - min_tools: 1
- [ ] Add mosdepth BED entry `chrM:448-1801` (MT-RNR1 ± 200 bp flank) to `pgx-bamstats.sh`
  - Note: chrM depth is typically very high (500–5000×) — flag if < 100× as insufficient for AF calling
- [ ] Add `GENE_LOCI` entry for MT-RNR1: `chrM:648-1,601 (GRCh38)`

**Complexity assessment:** High — requires a new caller paradigm (AF-based, not diplotype-based),
GATK in Docker, and new parsing logic. Suggest implementing as a standalone Phase 7b after CACNA1S (Phase 7a).

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
