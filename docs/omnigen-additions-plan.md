# OmniGen additions — implementation plan (pgx-suite)

> **Status (2026-07-11): SCOPE CHANGED — features 3 (HLA class II) and 4 (ABO) IMPLEMENTED
> in pgx-suite; features 1 (mtDNA) and 2 (Y haplogroup) MOVED TO OmniGen.** After validating
> that haplogrep3 reproduces both GIAB trios' maternal lineages exactly directly from the 30X
> WGS VCFs, the user decided both haplogroups live in OmniGen (classified downstream), not
> pgx-suite. The mtDNA/Y sections below are retained for historical/design context only — the
> corresponding pgx-suite code (`pgx-mito.sh`, `pgx-ychr.sh`, `ychr_isogg_markers.bed`, the
> `mito_haplogroup`/`ychr_haplogroup` rules, haplogrep3/yhaplo Dockerfile lines, the Ancestry
> & Lineage report card) has been removed.
>
> **Landed in pgx-suite:** `docker/omnigen_addons.py` (HLA-II + ABO logic), `pgx-hla2.sh`,
> `parse_abo_vcf`/`parse_arcashla` in `pgx-compare.py`, `genes.tsv` `arcashla` column + 4 rows
> (`HLA-DQA1/DQB1/DRB1`, `ABO`), the report **HLA Class II & Blood Group** section, and
> `test_omnigen_addons.py`. See `docs/CHANGES.md` (2026-07-11) for the shipped summary.
>
> **Update (2026-07-12): Feature 3 arcasHLA typing VALIDATED on GIAB HG002.** arcasHLA was
> provisioned and run end-to-end on the HG002 30–40X MHC reads (all 7 loci) against a
> reference built from IMGT/HLA 3.63.0. Class-II concordance vs the Chin et al. MHC-benchmark
> truth (Suppl. Table 4, clinical trio-phased): **DRB1 2/2, DQB1 2/2** at 2-field;
> **DQA1 1/2** (`03:01` exact, `01:05` vs truth `01:01` — same allele group). Class-I
> cross-check: A 2/2, B/C allele-group-correct but 2-field-discordant (WGS depth is thin for a
> tool tuned to RNA-seq). See Feature 3 §(f) below. ABO still requires its reference/coordinate
> + known-control validation.

Planning doc for four new features that pgx-suite will **call/type from the aligned
GRCh38 BAM** and emit as small stable contract files for OmniGen to render. No new
germline caller is required — every feature works from the same BAM the suite already
consumes. This is a plan only; nothing here has been implemented.

## 0. Architecture recap (verified against the repo)

The pipeline is a **Snakemake** DAG (`docker/Snakefile`) driven by a single config
table `docker/genes.tsv`. Key facts that constrain the design:

- **Gene routing.** `genes.tsv` columns: `gene region pypgx stargazer aldy stellarpgx
  optitype mutserve vcf_check pypgx_sv stargazer_sv stargazer_control cyrius pharmcat`.
  `TARGET = config["genes"].split() or ALLGENES`. Adding a row makes a "gene" flow
  automatically into the per-gene compare + report path (see `docs/howto-add-a-gene.md`).
- **Caller rules are non-fatal.** Each caller rule writes a *status sentinel*
  `STAT(gene, tool)` = `genes/<gene>/logs/<tool>.status` holding the exit code; the rule
  always succeeds, so one caller failing never aborts the batch. `rule compare` reads the
  sentinels and globs whatever real outputs exist.
- **Special callers already exist as shell scripts** invoked by thin rules:
  `rule optitype` → `pgx-hla.sh` (HLA-A/HLA-B), `rule mutserve` → `pgx-mt.sh` (MT-RNR1).
  Both live in `docker/` and write into the gene dir `G(gene)` = `<OUT>/genes/<gene>`.
- **Per-gene contract.** `rule compare` runs `docker/pgx-compare.py` and writes:
  `genes/<gene>/<gene>_<sample>_comparison.tsv` (cols:
  `Gene Sample Build Tool Diplotype ActivityScore Phenotype Status SVMode`) and
  `genes/<gene>/<gene>_<sample>_detail.json` (carries the authoritative `verdict`, per-tool
  dict, `coverage`, `sv_mode`). OmniGen already reads the HLA `*_comparison.tsv`.
- **Report.** `rule report` → `docker/pgx-report.py --sample --output --genes-dir --bam
  --bam-stats --provenance` → `results/<S>/<S>_pgx_report.html`. It loops `TARGET`,
  loading `genes/<gene>/<gene>_<sample>_detail.json` (glob fallback). Clinical text comes
  from a CPIC DB plus `RISK_ALLELE_PATTERNS` (line ~137) and `GENE_REGION_DISPLAY`
  (line ~98/111) dicts, and `build_clinical_findings_section` / `build_landing`.
- **Containers.** Main image `pgx-suite:latest` (`Dockerfile`, multi-stage): Java 21 JRE,
  samtools, bcftools, tabix, mosdepth, sambamba, snakemake, pypgx, aldy, stargazer,
  nextflow, and `mutserve.jar` baked at `/usr/local/bin/mutserve.jar`. Heavy tools ship as
  **Apptainer SIFs** in `StellarPGx/containers/` (`optitype.sif`, `pharmcat.sif`,
  `stellarpgx-dev.sif`), run with `--privileged`. Host launcher `run_pgx_suite.sh` mounts
  resources and calls the Snakefile inside the container.
- **VCF-Check star-allele machinery** (the model for ABO): `parse_ugt1a1_vcf` /
  `parse_g6pd_vcf` / `parse_cacna1s_vcf` / `parse_ryr1_vcf` in `pgx-compare.py`. Pattern:
  extract region → `bcftools mpileup|call` indel-aware → `<gene>.vcf.gz` → parse against a
  per-gene `_<GENE>_SNP_ALLELES` def table → assign alleles. Registered in
  `_VCF_CHECK_PARSERS`; authority via the `AUTHORITATIVE` dict.

**Two output styles.** (1) CPIC/diplotype-shaped features (HLA class II, ABO) fit the
existing `genes.tsv` → compare → detail.json → report card path *and also* emit their own
top-level contract file. (2) Lineage features (mtDNA, Y haplogroup) are **not** CPIC
diplotypes and do not fit the concordance/verdict model — they get a dedicated small
contract file plus a new "Ancestry / Lineage" section in the report, not a gene card.

---

## Feature 1 — mtDNA haplogroup (maternal lineage)

### (a) Pipeline home
Extend the existing mutserve path. `pgx-mt.sh` already extracts **all** chrM reads
(`samtools view -b BAM chrM`) and the whole chrM reference (`samtools faidx REF chrM`) and
runs `java -jar mutserve.jar call --input chrM_sorted.bam --output <raw>.txt --reference
chrM_ref.fa --level 0.01`. mutserve therefore already calls the **whole 16.5 kb
mitochondrion**; only the *parsing* step is narrowed to MT-RNR1 pos 1494/1555. So the
whole-genome call already exists — we add haplogroup **classification** on top.

Two clean options:
- **Preferred: new script `docker/pgx-mito.sh` + new Snakefile rule `mito_haplogroup`.**
  Keep `pgx-mt.sh` (MT-RNR1 CPIC) untouched so the pharmacogenetic contract is unchanged;
  the new rule depends on the mutserve raw/VCF output (or re-runs the extract, cheaply).
  Add a `mito_haplogroup` column to `genes.tsv`'s MT-RNR1 row is **not** needed — instead
  wire the rule at the sample level (like `bamstats`/`provenance`) and add its output to
  `rule report`'s inputs.
- Alternative: fold classification into `pgx-mt.sh` and have it emit both files. Rejected —
  keeps CPIC and lineage concerns separable and independently testable.

Sample-level rule sketch (mirrors `rule bamstats`):
```python
rule mito_haplogroup:
    input:  STAT("MT-RNR1", "mutserve")        # ensures mutserve ran
    output: os.path.join(OUT, "mito", "haplogroup.tsv")
    shell:  "pgx-mito.sh {BAM} {SAMPLE} {OUT} --ref {REF}"
```
and add `os.path.join(OUT,"mito","haplogroup.tsv")` to `rule report.input`.

### (b) Tool + container
- **mutserve v2** (already baked) emits a VCF: run `mutserve.jar call ... --output
  sample.vcf.gz` (mutserve writes a bgzipped multi-sample VCF alongside the txt) — this is
  the native haplogrep input.
- **haplogrep** for classification. It is a single Java jar (Java 21 already present), so
  **bake `haplogrep.jar` into the Docker image exactly like `mutserve.jar`** (add a
  download+verify step to `Dockerfile` §4, fail loudly if the asset is not a real jar).
  Use **haplogrep 2.4.x** (`haplogrep classify --in sample.vcf.gz --format vcf --out
  hg.txt`) or **haplogrep3** (`haplogrep3 classify --tree phylotree-rcrs@17.2 --in ...`).
  haplogrep3 is preferred (current Phylotree 17.x, richer quality output) and is also a
  self-contained Java distribution. Prereq: mutserve/haplogrep expect **rCRS chrM**;
  GRCh38 `chrM` **is** rCRS (NC_012920), so no lift is needed — assert the contig name is
  `chrM` and length 16569 as a guard.

### (c) Contract output
`results/<S>/mito/haplogroup.tsv` — single data row, tab-separated:

| Column | Example | Notes |
|--------|---------|-------|
| `sample` | `HG002` | |
| `haplogroup` | `H2a2a1` | haplogrep top hit |
| `quality` | `0.9721` | haplogrep quality score (0–1) |
| `defining_variants` | `263G;750G;1438G;...` | `;`-joined polymorphisms vs rCRS |

`pgx-mito.sh` runs haplogrep, then a small embedded Python block (same style as
`pgx-mt.sh` step 4) reads haplogrep's output and writes the TSV. Also emit a sidecar
`mito/haplogroup.json` (`sample, haplogroup, quality, rank, defining_variants[], n_variants,
overall_mt_coverage`) for richer OmniGen rendering — TSV stays the stable minimal contract.

### (d) Report integration
mtDNA haplogroup is **not** a CPIC diplotype → add a new **"Ancestry & Lineage"** section
to `pgx-report.py` `build_landing`. Add a helper `mito_lineage_card(out_dir)` that reads
`mito/haplogroup.tsv` and renders haplogroup + quality + a short "maternal lineage /
haplogroup macro-clade" caption. Gate on file presence so older runs still render. No
`RISK_ALLELE_PATTERNS`/CPIC entry needed. Keep it clearly informational (ancestry, not
clinical).

### (e) Testing
- **In-session (no full run):** unit-test the haplogrep-output → TSV parser with a captured
  fixture (feed a canned haplogrep `.txt`/`.json`); assert column schema and that a known
  GIAB VCF classifies to its published haplogroup (HG002 ≈ H2a2a1, HG003/HG004 known).
  Add `docker/test_mito.py` alongside `test_parsers.py`. Validate the contig guard logic
  standalone.
- **Needs full pipeline run:** end-to-end mutserve→haplogrep on a real BAM (container +
  Java + jars). Validate against a sample with a known mtDNA haplogroup (GIAB trio have
  published haplogroups) and confirm `results/<S>/mito/haplogroup.tsv` appears and the
  report renders the lineage card.

---

## Feature 2 — Y-chromosome haplogroup (paternal lineage)

### (a) Pipeline home & feasibility
**Feasibility flag: MODERATE fit, feasible with caveats.** pgx-suite is a *targeted*
region-extraction architecture, and Y-haplogroup calling needs genotypes at thousands of
ISOGG/yhaplo SNP positions spread across chrY — larger than a typical gene region but still
a bounded, BAM-only region extraction, so it fits the "extract region → call → classify"
mold. Caveats to document prominently:
- **Females (X46,XX) have ~no chrY reads** → the feature must gate on chrY coverage and
  emit `NA (no Y chromosome / insufficient chrY coverage)` rather than a spurious call.
- GRCh38 chrY has **PAR1/PAR2 + XTR** pseudoautosomal/homologous regions that attract
  X-derived reads; restrict genotyping to the **Y-unique (MSY)** ISOGG marker positions to
  avoid X contamination.
- Resolution scales with depth; 30× WGS gives good sub-clade resolution, low-coverage
  samples resolve only to macro-haplogroup.

**If judged a poor fit for a targeted PGx suite**, note the alternative: **poly-suite**
(the sibling polygenic/ancestry suite) is a more natural home for genome-wide lineage, and
OmniGen could consume the Y contract from there instead. But since pgx-suite is the
assigned home, the plan below implements the pgx-suite route.

Home: new `docker/pgx-ychr.sh` + sample-level rule `ychr_haplogroup` (same wiring style as
`mito_haplogroup`), output added to `rule report.input`.

### (b) Tool + container
- **yhaplo** (Poznik) — pip-installable pure Python; classifies from a VCF of chrY SNP
  genotypes at its bundled ISOGG marker set. **Add to the venv in the Dockerfile builder
  stage** (`pip install yhaplo`) — no separate SIF needed, keeping with the "small Java/py
  tools baked in, heavy tools as SIF" convention. Alternatives if yhaplo packaging is
  troublesome: **snappy** or **Y-LineageTracker** (both heavier; would go in a SIF).
- Genotyping step (in `pgx-ychr.sh`): ship a `docker/ychr_isogg_markers.bed` (yhaplo's
  marker positions lifted to GRCh38, MSY only), then
  `bcftools mpileup -R ychr_isogg_markers.bed -f REF -r chrY BAM | bcftools call -m` →
  `ychr/ychr.vcf.gz` (indel-aware not required; these are SNP markers). Then
  `yhaplo --input ychr.vcf.gz --output-dir ychr/`.
- **Coverage gate:** `samtools coverage -r chrY:2781480-56887902 BAM` (MSY span, GRCh38);
  if mean depth < threshold (e.g. < 1×) → write the `NA` contract and skip yhaplo.

### (c) Contract output
`results/<S>/ychr/haplogroup.tsv` — single data row:

| Column | Example | Notes |
|--------|---------|-------|
| `sample` | `HG002` | |
| `haplogroup` | `R1b1a1b1a1a2` | yhaplo long-form label |
| `haplogroup_short` | `R-M269` | terminal-SNP short name |
| `quality` | `0.94` | fraction of derived-consistent markers / yhaplo score |
| `n_markers_used` | `812` | informative markers genotyped |
| `chrY_mean_depth` | `28.3` | for QC / female detection |
| `status` | `ok` \| `no_call_female` \| `low_coverage` | |

Plus sidecar `ychr/haplogroup.json` with the derived-marker path for audit.

### (d) Report integration
Render inside the same new **"Ancestry & Lineage"** section next to mtDNA. Add
`ychr_lineage_card(out_dir)` reading `ychr/haplogroup.tsv`; when `status != ok` show an
explicit "No Y chromosome detected / insufficient coverage" note (correctly handles female
and low-coverage samples). Informational only.

### (e) Testing
- **In-session:** unit-test the yhaplo-output parser and the coverage-gate/female-branch
  logic with fixtures (canned yhaplo output + a synthetic zero-chrY-coverage case).
  `docker/test_ychr.py`. Confirm `ychr_isogg_markers.bed` parses and lifts cleanly.
- **Needs full run:** genotyping + yhaplo on a real male BAM (HG002 → R-M269 class per
  published data) and on a female BAM (→ `no_call_female`). Verify contract + report card.

---

## Feature 3 — HLA class-II typing (HLA-DQA1 / DQB1 / DRB1)

### (a) Pipeline home
OptiType (`pgx-hla.sh`) is **class I only**. Add a class-II-capable typer as a new special
caller alongside it. Because a single typer call resolves all HLA genes at once, run it
**once per sample** (like PharmCAT) rather than per-gene:
- New `docker/pgx-hla2.sh` + a **sample-level** rule `hla2` that types DQA1/DQB1/DRB1 (and
  optionally class I for cross-check) in one pass, writing `hla2/<sample>.genotype.json`.
- Add three rows to `genes.tsv`: `HLA-DQA1`, `HLA-DQB1`, `HLA-DRB1`, each with a **new
  support column `hla2=1`** (region = MHC `chr6:28510020-33480577`, all other caller flags
  `0`). Update the column list + `docs/howto-add-a-gene.md` table + `test_genes_config.py`
  frozen dicts. The per-gene `compare` for each of these genes reads the shared
  `hla2/<sample>.genotype.json` (sample-level dependency, exactly like the PharmCAT wiring
  in `compare_inputs`).

### (b) Tool + container
- **Preferred: arcasHLA.** Types class I **and** class II directly from a WGS BAM, is
  lighter than HLA-LA, and outputs a simple per-gene genotype JSON. Ship as an Apptainer
  SIF beside OptiType:
  `apptainer pull --name arcashla.sif docker://quay.io/biocontainers/arcas-hla:<ver>`,
  placed in `StellarPGx/containers/arcashla.sif`, mounted read-only like `optitype.sif`.
  arcasHLA needs its IMGT/HLA reference DB built once (`arcasHLA reference --version …`);
  bake the built DB into the SIF or mount a prepared DB dir. Flow in `pgx-hla2.sh`:
  extract MHC reads (reuse the `pgx-hla.sh` extraction, `chr6:28,510,020-33,480,577`) →
  `arcasHLA extract` → `arcasHLA genotype -g A,B,C,DQA1,DQB1,DRB1` → JSON.
- **Alternative: HLA-LA** — higher accuracy on class II, but ~30–40 GB RAM and much longer
  runtime; would be its own SIF and a heavier resource block. Note it as the fallback if
  arcasHLA accuracy is insufficient for DRB1.
- Runs `--privileged` like the other SIFs.

### (c) Contract output
**Extend the existing `genes/HLA-*/*_comparison.tsv` pattern** (the schema OmniGen already
reads). A new parser `parse_arcashla(output_dir, gene, sample)` in `pgx-compare.py`
(registered so `HLA-DQA1/DQB1/DRB1` route to it, analogous to `parse_optitype`) reads
`hla2/<sample>.genotype.json` and emits, per gene, the standard row:
```
Gene            Sample  Build   Tool      Diplotype                       ActivityScore Phenotype Status SVMode
HLA-DQB1        HG002   GRCh38  arcasHLA  HLA-DQB1*06:02/HLA-DQB1*03:01    -             -         ok     no SVs expected
```
and the matching `genes/HLA-DQB1/HLA-DQB1_HG002_detail.json`. No new top-level contract
directory is required — this deliberately reuses the class-I contract so OmniGen needs zero
new reader logic; it just sees three more `HLA-*` genes.

### (d) Report integration
- Add `HLA-DQA1/DQB1/DRB1` to `GENE_REGION_DISPLAY` and give them clinical text.
- Extend `RISK_ALLELE_PATTERNS` with class-II risk lambdas so the **Key Clinical Findings**
  tiering fires:
  - **Celiac disease** — DQ2.5 (`DQA1*05:01` + `DQB1*02:01`) and DQ8 (`DQA1*03` +
    `DQB1*03:02`); note these are *haplotype* combinations across DQA1+DQB1, so the lambda
    must inspect both gene calls (add a small cross-gene evaluator in
    `build_clinical_findings_section`, which already receives all `genes_data`).
  - **Narcolepsy type 1** — `DQB1*06:02`.
  - **Type 1 diabetes** — DR3-DQ2 / DR4-DQ8 risk haplotypes (DRB1*03:01, DRB1*04 +
    DQB1*03:02); DQB1*06:02 is protective — encode both directions.
- Emit these as **informational/risk** findings (HLA disease associations, not CPIC drug
  guidance) with clear "susceptibility, not diagnostic" language.

### (e) Testing
- **In-session:** unit-test `parse_arcashla` against a captured `hla2/<sample>.genotype.json`
  fixture; assert comparison.tsv/detail.json schema for each of the three genes. Unit-test
  the celiac/narcolepsy/T1D cross-gene risk logic with synthetic diplotype combos
  (DQ2.5 present/absent, DQ8, DQB1*06:02). Extend `test_parsers.py` / add `test_hla2.py`.
  `test_genes_config.py` after adding the rows + `hla2` column.
- **Needs full run:** build/pull `arcashla.sif` + DB and type a real BAM; validate DQ calls
  against a sample with known class-II type (e.g. a GIAB/1000G sample with published
  IMGT-typed HLA, or the OptiType class-I concordance as a sanity anchor).

### (f) Validation — GIAB HG002 (2026-07-12) ✅ DONE

arcasHLA (RabadanLab, wrapper v0.2.0) provisioned locally: kallisto 0.44.0 + git-lfs +
python-3.10/numpy/scipy/pandas/biopython in dedicated conda envs; samtools/bedtools/pigz from
host. Reference built with `arcasHLA reference --rebuild` from **IMGT/HLA 3.63.0** `hla.dat`
(fetched from the EBI IPD-IMGT/HLA FTP because the ANHIG GitHub LFS budget was exhausted).
Input = the HG002 30–40X MHC reads extracted at `chr6:28,510,020-33,480,577`
(`HG002.hla.{1,2}.fq.gz`). Command:
`arcasHLA genotype -g A,B,C,DPB1,DQA1,DQB1,DRB1 R1 R2 -o hla2/` (temp under `/data/alvin/tmp`).

Truth = **Chin et al. 2020, "A Diploid Assembly-based Benchmark for Variants in the MHC",
Suppl. Table 4** (HG002 clinical trio-phased typing). Concordance at 2-field
(allele-group:protein), unordered:

| Locus | arcasHLA (2-field) | Truth | 2-field match | Reads | Class |
|-------|--------------------|-------|---------------|-------|-------|
| HLA-DRB1 | 10:01 / 04:02 | 10:01 / 04:02 | **2/2 ✅** | 134 | II (target) |
| HLA-DQB1 | 03:02 / 05:01 | 03:02 / 05:01 | **2/2 ✅** | 137 | II (target) |
| HLA-DQA1 | 03:01 / 01:05 | 01:01 / 03:01 | 1/2 (⚠ 01:05 vs 01:01, same group) | 167 | II (target) |
| HLA-DPB1 | 04:01 / 677:01 | — (not in truth table) | n/a (04:01 plausible; 677:01 low-conf) | 123 | II |
| HLA-A | 01:01 / 26:01 | 01:01 / 26:01 | 2/2 ✅ | 86 | I (cross-check) |
| HLA-B | 35:15 / 38:02 | 35:08 / 38:01 | 0/2 (groups 35/38 right, protein wrong) | 77 | I (cross-check) |
| HLA-C | 05:01 / 12:233 | 04:01 / 12:03 | 0/2 (C\*12 group right on one) | 54 | I (cross-check) |

**Verdict:** the class-II validation target passes — DRB1 and DQB1 are exact at 2-field and
DQA1 is allele-group-correct on both alleles. DQA1's `01:05` vs `01:01` and the class-I B/C
protein-field errors reflect thin WGS coverage over the MHC (54–86 reads/locus) for a tool
tuned to high-depth RNA-seq; **OptiType remains the primary class-I typer** (its HG002 A call
already agrees with truth). arcasHLA's own confidence signals corroborate: the class-II loci
have the clearest heterozygosity separation (DQA1 0.99, DQB1 0.96) and the most reads.
Raw output: `results/HG002/hla2/HG002.genotype.json` (+ `.genes.json`, `.genotype.log`).
Contract rows: `results/HG002/genes/HLA-{DQA1,DQB1,DRB1}/*_comparison.tsv` + `*_detail.json`.

---

## Feature 4 — ABO blood-type typer

### (a) Pipeline home
ABO is **SNP + indel typing, not a CNV** — model it directly on the existing VCF-Check
star-allele machinery (`parse_ugt1a1_vcf` is the closest template because it already mixes
an indel call with SNP-def-table calls). Two coordinated pieces:
1. **Reuse `genes.tsv` + VCF-Check** for the per-gene detail/comparison + report card: add
   an `ABO` row (`region chr9:133,250,000-133,280,000`, `vcf_check=1`, everything else `0`),
   a new parser `parse_abo_vcf` in `pgx-compare.py`, register it in `_VCF_CHECK_PARSERS`,
   and (to make it authoritative) add `"ABO": "VCF-Check"` to `AUTHORITATIVE`. The generic
   `vcf_check` rule already produces `ABO.vcf.gz` via indel-aware `bcftools mpileup|call`.
2. **Emit the dedicated contract** `results/<S>/abo/abo_type.tsv`. Simplest: have
   `parse_abo_vcf` (or a tiny sample-level `pgx-abo.sh` post-step) write the extra file in
   addition to the standard detail/comparison outputs.

### (b) Tool + container
No new tool — **bcftools/samtools already in the base image**, indel-aware calling is the
exact path VCF-Check genes already use. Only new data artifacts:
- `docker/abo_allele_defs.json` — the ABO allele-definition table (positions, ref/alt,
  rsIDs) — companion to `vcf_check_variants.json`, guarded by a new
  `test_abo.py` (mirrors `test_vcf_check.py`).

Core defining variants (GRCh38, ABO on chr9q34.2):
- **rs8176719** — `c.261delG` single-base **deletion** → **O** (frameshift, non-functional).
- **rs8176746** (`c.796C>A`, Leu266Met) and **rs8176747** (`c.803G>C`, Gly268Ala) — **A vs
  B** discrimination.
- (Extendable: rs8176743 for A2/weak subgroups later; keep v1 to O/A/B.)

Genotype logic: per haplotype, `c.261delG` present ⇒ O; else A vs B from rs8176746/747;
combine the two haplotype calls (unphased — report as compound genotype like the UGT1A1
`*x(het)+*y(hom)` style) into `A/O`, `B/O`, `A/B`, `O/O`, `A/A`, `B/B`.

### (c) Reference-allele question (must be resolved before coding, and validated)
State explicitly which ABO allele the **GRCh38 reference encodes**, because every call is
relative to it:
- At **rs8176719** the reference genome **carries the G** (the deletion is defined
  *relative to* ref), i.e. the reference is the **non-deleted / functional** state → ref is
  **not O**.
- At **rs8176746 / rs8176747** the reference bases are the **A-consistent** alleles
  (`c.796C`, `c.803G`).
- Therefore GRCh38 reference ABO is expected to type as a **functional A-like allele
  (≈ ABO\*A1.01), i.e. non-O**. **This must be empirically confirmed**, not assumed: as the
  first implementation step, `samtools faidx GRCh38 chr9:<pos>` the three positions and
  record the actual reference bases in `abo_allele_defs.json` provenance. A homozygous-ref
  sample must then type as **A/A**, and a true O/O sample must show the deletion at both
  haplotypes — this is the primary correctness check (a common failure mode is inverting
  the O logic because "reference = O" is a frequent but wrong assumption).

### (d) Contract output
`results/<S>/abo/abo_type.tsv` — single data row:

| Column | Example | Notes |
|--------|---------|-------|
| `sample` | `HG002` | |
| `ABO_type` | `O` | serological group: A / B / AB / O |
| `alleles` | `ABO*O.01.01/ABO*O.01.01` | genotype (two haplotypes, unphased) |
| `confidence` | `high` | `high`/`medium`/`low` from depth + variant clarity at the 3 sites |

Sidecar `abo/abo_type.json` (`sample, ABO_type, alleles, per_site:[{rsid,pos,gt,depth,ref,alt}],
reference_allele_note, confidence`). The `genes/ABO/ABO_<sample>_comparison.tsv` +
`detail.json` still emit via the normal path so ABO also appears as a standard report card.

### (e) Report integration
- Add `ABO` to `GENE_REGION_DISPLAY` (`chr9:133,255,...`) so its region shows.
- Because `ABO` is `vcf_check`/`AUTHORITATIVE`, it renders as a normal **gene diplotype
  card** automatically (verdict from detail.json). Add a short blood-group interpretation
  string (A/B/AB/O + Rh-not-covered caveat) to the CPIC/clinical text block, and optionally
  surface ABO type in **Key Clinical Findings** as informational (transfusion / blood-group
  context — explicitly *not* a substitute for serological cross-match).

### (f) Testing
- **In-session:** unit-test `parse_abo_vcf` against synthetic `ABO.vcf.gz` fixtures covering
  O/O (del/del), A/O, A/B, A/A, homozygous-ref → must be A/A; assert both the standard
  comparison/detail rows and the `abo/abo_type.tsv` contract. Add `test_abo.py`. Confirm the
  three reference bases via `samtools faidx` (in-session, needs the GRCh38 FASTA which is a
  repo symlink — cheap, no pipeline).
- **Needs full run + control:** run the `vcf_check` path on a real BAM and **validate
  against a known-ABO-type control** (a sample with serologically/orthogonally known ABO —
  e.g. a GIAB or 1000G sample with a published/curated ABO type). Do not ship until the
  control types correctly and the reference-allele assumption in (c) is confirmed on real
  data.

---

## Cross-cutting work & sequencing

1. **genes.tsv column additions.** New `hla2` column (Feature 3) and new rows
   (`HLA-DQA1/DQB1/DRB1`, `ABO`). Update: `docker/genes.tsv`, the column-order docstring in
   `Snakefile`, `docs/howto-add-a-gene.md` table, and the frozen `ORIG_SUPPORT`/`ORIG_COORDS`
   in `docker/test_genes_config.py`.
2. **Dockerfile.** Bake `haplogrep.jar` (F1) next to `mutserve.jar` with download+verify;
   `pip install yhaplo` into the builder venv (F2). Pull `arcashla.sif` into
   `StellarPGx/containers/` and build/stage its IMGT DB (F3). No image change for ABO (F4).
3. **New scripts / rules.** `pgx-mito.sh` (+ `mito_haplogroup` rule), `pgx-ychr.sh` (+
   `ychr_haplogroup` rule, `ychr_isogg_markers.bed`), `pgx-hla2.sh` (+ `hla2` rule). Wire
   all lineage/sample-level outputs into `rule report.input`.
4. **pgx-compare.py.** Add `parse_arcashla` (route the 3 HLA-II genes), `parse_abo_vcf`
   (+ `_VCF_CHECK_PARSERS["ABO"]`, `AUTHORITATIVE["ABO"]`), `abo_allele_defs.json`.
5. **pgx-report.py.** New "Ancestry & Lineage" section (mtDNA + Y cards); extend
   `GENE_REGION_DISPLAY`, `RISK_ALLELE_PATTERNS` (+ cross-gene HLA-II celiac/T1D evaluator),
   clinical text for the new genes.
6. **Tests.** `test_mito.py`, `test_ychr.py`, `test_hla2.py`, `test_abo.py`; extend
   `test_parsers.py`, `test_genes_config.py`, `test_vcf_check.py` sibling.
7. **run_pgx_suite.sh / provenance.** Add the new SIF mounts (arcasHLA) and tool-version
   lines (haplogrep, yhaplo, arcasHLA) to `rule provenance`.

### What is testable in-session vs needs a full pipeline run
- **In-session (parsers, schema, logic, small CLI checks):** every `*_comparison.tsv` /
  contract-file writer and every parser can be unit-tested against captured tool-output
  fixtures without running the callers. The ABO reference-base check
  (`samtools faidx`) and genes.tsv config regression run in-session.
- **Needs a full pipeline run (container + heavy tools + real BAM):** mutserve→haplogrep,
  chrY genotyping→yhaplo, arcasHLA typing, and the bcftools ABO call — plus the end-to-end
  report render and all control-sample validations (GIAB haplogroups, known HLA-II type,
  known ABO type).

## OmniGen-facing contract summary
| Feature | New contract file(s) | Reuses existing pattern |
|---------|----------------------|-------------------------|
| mtDNA haplogroup | `results/<S>/mito/haplogroup.tsv` (+ `.json`) | new lineage file |
| Y haplogroup | `results/<S>/ychr/haplogroup.tsv` (+ `.json`) | new lineage file |
| HLA class II | `genes/HLA-{DQA1,DQB1,DRB1}/*_comparison.tsv` (+ `detail.json`) | **same as class-I HLA** |
| ABO | `results/<S>/abo/abo_type.tsv` (+ `.json`) and `genes/ABO/*_comparison.tsv` | VCF-Check star-allele |
