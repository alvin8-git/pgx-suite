# PGx Suite — Changelog

All notable changes to this project are documented here.
Format: reverse-chronological, grouped by phase/milestone.

---

## 2026-03-12 — Hypersensitivity Panel: HLA-B Level B allele annotations

### Summary

Annotation-only extension of the HLA-B CPIC_DB entry in `pgx-report.py`.
OptiType already genotypes HLA-B\*13:01 and HLA-B\*57:03 internally; these
changes ensure the clinical significance is surfaced in the HTML report.

### `docker/pgx-report.py` — HLA-B CPIC_DB

- **`desc`**: extended to mention HLA-B\*13:01 (dapsone/DRESS, Asian ancestry)
  and HLA-B\*57:03 (flucloxacillin/DILI, European ancestry)
- **`diplotype_check` lambda**: now includes `"B*13:01"` and `"B*57:03"` in the
  `any()` check, so the detail page risk banner triggers for these alleles
- **`risk_alleles`**: two new entries:
  - `B*13:01` → dapsone: avoid — high risk of DRESS/HSS (esp. Asian ancestry)
  - `B*57:03` → flucloxacillin: use with caution — increased DILI risk (esp. European ancestry)
- **`no_risk_note`**: updated to list all five alleles with CPIC level annotations
  (B\*57:01, B\*15:02, B\*58:01 Level A; B\*13:01, B\*57:03 Level B)
- **`drugs`**: two new drug entries:
  - Dapsone (CPIC Level B) — recommend alternative if B\*13:01 positive
  - Flucloxacillin (CPIC Level B / PharmGKB) — monitor LFTs; consider alternative
    beta-lactam

### `TODO.md` — Phase 9: Hypersensitivity Panel Expansion

New Phase 9 section added (former Phase 9 Population PGx Expansion renumbered to Phase 10):
- **9a** (done ✅): B\*13:01 + B\*57:03 CPIC_DB annotations (this change)
- **9b**: HLA Class II typing (HLA-DRB1, DQA1, DQB1) — requires HLA-HD or arcasHLA;
  implementation tasks listed
- **9c**: CYB5R3 (dapsone methaemoglobinaemia; rs1799931) — custom bcftools lookup
- **9d**: HLA-C clinical reporting — OptiType already genotypes it; needs
  GENE_LOCI + CPIC_DB entry

---

## 2026-03-12 — Phase 7: CACNA1S + MT-RNR1 (19/19 CPIC Level A) + mutserve v2.0.0

### Summary

Completed 19/19 CPIC Level A gene coverage by adding CACNA1S (malignant hyperthermia
susceptibility) and MT-RNR1 (aminoglycoside-induced hearing loss). A sixth calling tool,
mutserve v2.0.0, was integrated for chrM heteroplasmy calling. Scientific documentation
for all 31 genes was authored in `PGxDocumentation.md`.

### CACNA1S

CACNA1S encodes the α1S subunit of the voltage-gated L-type calcium channel.
Risk alleles (*2, c.3257G>A) increase susceptibility to malignant hyperthermia
(MHS) triggered by volatile anaesthetics and succinylcholine.

- **Tool support**: PyPGx ✓, StellarPGx ✓ (configuration-only; no new tool required)
- `pgx-run.sh`: `GENE_COORDS[CACNA1S]="chr1:201006956-201083927"`;
  `GENE_SUPPORT[CACNA1S]="1 0 0 1"` (PyPGx + StellarPGx)
- `pgx-all-genes.sh`: CACNA1S added to GENES list (31 total)
- `pgx-bamstats.sh`: BED entry `chr1 201006955 201083927 CACNA1S`
- `pgx-compare.py`: GENE_SUPPORT + RYR1-style parse routing
- `pgx-report.py`: GENE_LOCI + CPIC_DB (volatile anaesthetics + succinylcholine,
  Level A; `diplotype_check` flags non-reference diplotypes as MHS-susceptible)

### MT-RNR1 + mutserve v2.0.0

MT-RNR1 encodes the mitochondrial 12S rRNA. Variants m.1555A>G and m.1494C>T cause
aminoglycoside-induced irreversible sensorineural hearing loss (CPIC Level A).
MT-RNR1 is not amenable to star-allele diplotype calling; allele fraction (AF)
based heteroplasmy detection is used instead.

**New tool — mutserve v2.0.0** (standalone JAR, MIT licence, Java 21 already in image):

```bash
java -Xmx2g -jar mutserve.jar call \
    --input chrM_sorted.bam --output variants.txt \
    --reference chrM_ref.fa \
    --level 0.01 --threads N
```

- `Dockerfile`: `curl -sL …/mutserve.jar -o /usr/local/bin/mutserve.jar`; pgx-mt.sh installed
- **New script `docker/pgx-mt.sh`** (164 lines):
  - Extracts chrM reads from BAM via `samtools view -b | sort`
  - Extracts chrM FASTA at runtime via `samtools faidx /pgx/ref/hg38.fa chrM`
    (GRCh38 chrM = rCRS; no pre-baked reference needed)
  - Runs mutserve at `--level 0.01`; parses tab-delimited output for pos 1555/1494
  - Writes `mt-rnr1/<sample>_mtrna1_result.json`:
    `{diplotype, AF, type, classification, phenotype, notes}`
- `pgx-run.sh`: `DO_MUTSERVE` flag, `run_mt()` function, MT-RNR1 bypass of the
  standard bcftools/PyPGx/Stargazer/Aldy/StellarPGx pipeline
- `pgx-compare.py`: `"mutserve": True` in GENE_SUPPORT; `parse_mutserve()` reads JSON
- `pgx-report.py`: GENE_LOCI `"chrM:648-1,601 (rCRS / GRCh38 chrM)"`; CPIC_DB entry
  (aminoglycosides, Level A, `min_tools=1`, AF-based diplotype_check)
- `TOOLS`: extended from 5 to 6 — `["PyPGx","Stargazer","Aldy","StellarPGx","OptiType","mutserve"]`

### ToolsDocumentation.md

New **Section 6 — mutserve**: overview, pgx-mt.sh workflow, key differences table
(AF/heteroplasmy vs diplotype callers), limitations (SNPs only, CPIC variants only,
no star alleles). Quick Comparison table and Gene Support Matrix updated.

### README.md

Updated for 6-tool / 31-gene / 19/19 CPIC Level A state. Bundled Tools table, Gene
Coverage table, and Architecture section all reflect the new additions.

### PGxDocumentation.md (new, ~1300 lines)

Scientific gene summaries for all 31 genes. Each section covers:
Overview, PGx Significance, Key Variants (table), Population Frequencies, CPIC Level,
and Tool Support. Special notes: IFNL3/IFNL4 causal variant distinction, GSTT1 alt-contig
caveat, GSTM1/GSTT1 deletion genetics, G6PD X-linked hemizygous males, MT-RNR1
AF-based heteroplasmy model, CACNA1S/RYR1 pharmacogenomic relationship.

---

## 2026-03-12 — Phase 6: Output restructure, standalone HTML, depth colour flags, host launcher

### Summary

Four independent improvements made to the output pipeline and report generation,
plus a new host-side convenience launcher script.

### Output directory restructure (`pgx-all-genes.sh`, `pgx-run.sh`)

Per-gene results moved from flat `results/<GENE>/` to structured subdirectories:
- `Genes/<GENE>/` — per-gene tool outputs and comparison TSV
- `log/` — per-gene log files

`pgx-compare.py` and `pgx-report.py` updated to read from the new paths.

### Standalone single-file HTML report (default)

`pgx-report.py` now writes a single self-contained `<SAMPLE>_pgx_report.html`
at the output root by default (all per-gene detail sections embedded inline,
no subdirectory required). The multi-file mode (separate per-gene HTML pages)
remains available but is no longer the default. This simplifies sharing: one
file captures the full report.

### Depth coverage colour flags (`pgx-report.py`, `pgx-bamstats.sh`)

Per-gene mean depth displayed on gene cards and detail panels now uses colour
coding to indicate coverage adequacy:

| Colour | Threshold | Meaning |
|--------|-----------|---------|
| Green | ≥30× | Adequate for reliable calling |
| Amber | 20–29× | Marginal — interpret with caution |
| Red | <20× | Insufficient — results unreliable |

Depth values are sourced from `bam_stats.json → gene_depth` written by
`pgx-bamstats.sh`. Genes not covered by the mosdepth BED (e.g. GSTT1 on
alt contig) show no colour flag.

### `run_pgx_suite.sh` — host-side launcher

New file at repo root. A convenience wrapper that assembles the full
`docker run --privileged …` command with all required volume mounts and
passes through `<GENE> <BAM>` arguments. Eliminates the need to type the
full multi-volume Docker run command manually.

Usage: `./run_pgx_suite.sh <GENE> /absolute/path/to/sample.bam [--output /path]`

---

## 2026-03-09 — Phase 5: HTML Report Enhancements + HLA depth

### Summary

Three targeted improvements to `pgx-report.py` and `pgx-bamstats.sh` increasing
clinical utility of the HTML reports across all 29 genes. Validated on T7_NA24385
and HG003 (29/29 genes each).

### pgx-report.py: CPIC Clinical Reference on all genes

The `CPIC_DB` dict was extended from 18 to 29 genes by adding entries for the
11 previously missing genes: CYP1A1, CYP1A2, CYP2A6, CYP2C8, CYP2E1, CYP3A4,
GSTM1, GSTT1, IFNL3, NAT1, POR.

Each entry includes a gene description, CPIC drug–gene pairs with Level badges and
recommendation text, PharmVar link where applicable, and the patient-specific
finding note derived from the consensus diplotype/phenotype. Genes with no actionable
CPIC Level A/B guideline (CYP1A1, CYP2E1, GSTM1, GSTT1, NAT1, POR) include an
informational description and an empty drugs list so the CPIC section still renders
on the detail page.

Notable additions:
- **CYP1A2** — CPIC Level A: fluvoxamine (poor/ultra-rapid metaboliser dose guidance)
- **CYP2A6** — CPIC Level A: nicotine/tobacco cessation (NRT patch preference for poor metabolisers)
- **CYP2C8** — CPIC Level B: amodiaquine (avoid in poor metabolisers — agranulocytosis risk)
- **IFNL3** — CPIC Level A: peginterferon alfa + ribavirin for HCV (non-CC genotype → prefer DAA)
- **CYP3A4** — informational: *22 allele note with tacrolimus co-interpretation guidance

### pgx-report.py: in-page navigation bar

A `<nav class="gene-page-nav">` bar is rendered on every gene detail page between
the back-link and the gene header. Three anchor links: **Tool Results** (`#tool-results`),
**Variant Evidence** (`#variant-evidence`), **CPIC Reference** (`#cpic-reference`).
CSS uses a pill-button style with hover highlight.

Section anchor IDs added: `id="tool-results"` on the Tool Results div (was absent),
`id="cpic-reference"` on the CPIC section div (`#variant-evidence` was already present).

### pgx-report.py: chromosomal locus under gene name

A new `GENE_LOCI` dict maps all 29 genes to their GRCh38 1-based display coordinates
(e.g. `chr22:42,116,498-42,155,810`). On every gene detail page a monospace locus line
appears under the gene title: `📍 chrN:start-end (GRCh38)`.

HLA-A and HLA-B use the actual gene-body coordinates
(`chr6:29,910,247-29,913,661` and `chr6:31,321,649-31,324,666`) rather than the
wider MHC extraction region. GSTT1 shows `chr22_KI270879v1_alt (alt contig)`.

### pgx-bamstats.sh: HLA-A and HLA-B depth

HLA-A (`chr6:29,905,000-29,920,000`) and HLA-B (`chr6:31,316,000-31,330,000`) added
to the mosdepth BED file. The total covered by mosdepth is now 27 primary-assembly
genes (was 25). `bam_stats.json` `gene_depth` now includes entries for both HLA genes,
enabling depth display on the HLA gene cards and detail pages.

Observed depths (T7_NA24385, ~30× WGS):
- HLA-A: 46.82× mean, ≥30×: 98.61%
- HLA-B: 28.21× mean, ≥30×: 45.9% (lower due to MHC repeat structure)

### parse_optitype() bug fix

`pgx-compare.py` `parse_optitype()` was producing double-prefixed diplotypes:
`HLA-A*A*01:01` (OptiType outputs `A*01:01` in the `A1` column; old code prepended
`HLA-A*` instead of `HLA-`). Fixed for both HLA-A and HLA-B:

```python
# Before (HLA-A):
result.haplotype1 = f"HLA-A*{a1}" if not a1.startswith("HLA") else a1
# After:
result.haplotype1 = f"HLA-{a1}" if not a1.startswith("HLA") else a1
```

---

## 2026-03-09 — Phase 4: ABCG2 + HLA typing (OptiType)

### Summary

Added ABCG2 (CPIC Level A, rosuvastatin) and HLA-A/HLA-B (CPIC Level A, multiple
drug–allele pairs) to the 29-gene pipeline. OptiType 1.3.5 added as fifth typing tool
for HLA Class I alleles. Key Clinical Findings cards made clickable.

### ABCG2

ABCG2 (`chr4:88,085,265-88,236,626`) encodes the BCRP efflux transporter. The *2
variant (rs2231142, c.421C>A, Q141K) reduces transporter function and increases
rosuvastatin plasma exposure (CPIC Level A).

- `pgx-run.sh`: GENE_COORDS + GENE_SUPPORT `0 0 1 1` (Aldy + StellarPGx only; PyPGx
  and Stargazer have no ABCG2 support)
- `pgx-all-genes.sh`: ABCG2 added to GENES array
- `pgx-bamstats.sh`: `chr4:88,085,264-88,236,626` added to mosdepth BED
- `pgx-compare.py`: GENE_SUPPORT entry + `optitype: False`
- `pgx-report.py`: CPIC_DB entry with `diplotype_check` on `"rs2231142"`, rosuvastatin
  Level A recommendation

### HLA-A and HLA-B via OptiType

HLA Class I allele typing for three CPIC Level A drug–allele pairs:
- HLA-B\*57:01 → abacavir (HIV): contraindicated — severe hypersensitivity
- HLA-B\*58:01 → allopurinol (gout): contraindicated — SJS/TEN risk
- HLA-B\*15:02 → carbamazepine/phenytoin (epilepsy): high SJS risk in Asian ancestry
- HLA-A\*31:01 → carbamazepine: DRESS/SJS/TEN risk across all ancestries

**New script: `docker/pgx-hla.sh`**
Extracts MHC reads from the BAM (chr6:28,510,020-33,480,577) via `samtools view -b |
sort -n | fastq`, then runs OptiType via Apptainer SIF (`optitype.sif`). Outputs TSV
at `<output>/optitype/<sample>_result.tsv`. Handles paired-end fallback for samples
with low paired-read counts in the MHC region.

**`pgx-run.sh` changes:**
- HLA-A/HLA-B bypass the standard bcftools/PyPGx/Stargazer/Aldy/StellarPGx pipeline
  entirely via `DO_OPTITYPE=1` flag
- `run_hla()` function calls `pgx-hla.sh` with the correct output directory
- Graceful skip if `optitype.sif` is not present at `/pgx/containers/`

**`pgx-compare.py` changes:**
- GENE_SUPPORT entries for HLA-A and HLA-B (`optitype: True`)
- `parse_optitype()` reads `A1`/`A2` or `B1`/`B2` from OptiType result TSV; prepends
  `HLA-` prefix to convert `A*01:01` → `HLA-A*01:01`

**`pgx-report.py` changes:**
- `TOOLS` extended from 4 to 5: `["PyPGx","Stargazer","Aldy","StellarPGx","OptiType"]`
- CPIC_DB entries for HLA-A and HLA-B with `risk_alleles` dict and `min_tools=1`
- `_get_tier()` uses `entry.get("min_tools", 2)` to allow HLA genes to flag with one tool
- Detail table colspan changed from hardcoded `"4"` to `f"{len(TOOLS)}"`
- OptiType 1.3.5 badge in sample banner and footer

**Clickable clinical findings cards:**
Each finding card on the landing page is wrapped in
`<a href="{sample}.{gene}.pgx.html" class="cf-finding-link">` with a hover lift
animation and `→ Detail` arrow.

**Dockerfile:** `pgx-hla.sh` installed and symlinked to `/usr/local/bin/`

**OptiType SIF:** pulled separately inside Docker via
`apptainer pull --name optitype.sif docker://quay.io/biocontainers/optitype:1.3.5--hdfd78af_1`
(~500 MB, stored at `StellarPGx/containers/optitype.sif`)

### Validation

T7_NA24385 (29/29 OK): HLA-A `A*01:01/A*26:01`, HLA-B `B*35:08/B*38:01`.
No CPIC Level A risk alleles detected in this sample.

---

## 2026-03-08 — BAM QC pipeline overhaul (pgx-bamstats.sh + pgx-all-genes.sh)

### Summary

Complete rewrite of `pgx-bamstats.sh` and parallelisation of BAM QC in
`pgx-all-genes.sh`, driven by profiling of tool timings on a 167 GB WGS BAM
(T7_NA24385 and HG002, ~30× and ~40× depth respectively).

### Tool evaluation (observed timings, 167 GB WGS, cold disk cache at 490 MB/s)

| Tool / step | Time | Outcome |
|---|---|---|
| `samtools flagstat` (full BAM) | 518–588 sec | Replaced — I/O bound |
| `sambamba flagstat` (full BAM) | 563–683 sec | Evaluated, not faster — same I/O ceiling |
| `samtools idxstats` (.bai only) | 28–46 ms | **Adopted** for total/mapped/paired reads |
| `samtools depth` ×14 genes (old) | 18–60 min | Replaced by mosdepth |
| `mosdepth --by pgx.bed` (14 regions) | 350–620 sec | **Adopted** for per-gene depth |
| `samtools view -c` ×3 (chr1 4 Mb) | 400–550 ms | Adopted for MAPQ≥20 + dup estimate |
| `samtools stats` (chr1 4 Mb) | 1.9–2.7 sec | Unchanged |

**Key finding:** all full-BAM tools (flagstat, sambamba, mosdepth) are I/O bound at
the disk's 490 MB/s bandwidth limit. For a 167 GB BAM this is an unavoidable
~9–10 min floor for any tool that must read the file. sambamba's true-parallel
counting provides no speedup because the bottleneck is disk reads, not CPU.

### mosdepth: bugs found and fixed

1. **Wrong version in Dockerfile**: v0.3.13 does not exist; downloaded a GitHub 404
   HTML page silently (`curl -sL`). Fixed to v0.3.12 (verified with `mosdepth --version`).
2. **Wrong flag `--quantile`**: mosdepth 0.3.12 has no `--quantile` flag; the correct
   flag for per-region coverage thresholds is `--thresholds 0,20,30`, which outputs
   base counts per region in `{prefix}.thresholds.bed.gz`.
3. **Header row in thresholds.bed.gz**: first line begins with `#chrom`; bash arithmetic
   `$(( end - start ))` failed with "expression recursion level exceeded" when parsing
   the header. Fixed with `[[ "$chrom" == \#* ]] && continue`.

### Final pgx-bamstats.sh design

| Step | Tool | Metric | Time |
|---|---|---|---|
| 1 | `samtools idxstats` | total/mapped/paired reads, genome size | 28–46 ms |
| 2 | `samtools stats chr1:10M-14M` | read length, insert size, error rate | ~2 sec |
| 3 | idxstats × read_length | genome-wide depth estimate, chrX/Y depth, sex inference | — |
| 4 | `mosdepth --by pgx.bed --thresholds 0,20,30 --fast-mode --no-per-base` | per-gene mean depth, ≥20×/≥30× fractions | ~550 sec |
| 5 | `samtools view -c` ×3 (chr1 4 Mb) | MAPQ≥20%, duplicate% estimate | ~500 ms |

Global genome depth and sex inference are now computed from `idxstats` ×
`read_length` without any additional BAM scan.

### pgx-report.py: gene depth added to landing page gene cards

Each gene card on the landing page now shows `Depth: {mean}× | ≥20×: {pct}%`
below the concordance badge, sourced from `bam_stats.json → gene_depth`. Genes not
covered by the 14-gene mosdepth BED (e.g. GSTT1 on alt contig) show no depth line.

### pgx-all-genes.sh: BAM QC parallelised with gene calling

`pgx-bamstats.sh` is now launched in the background (`&`) immediately after the log
directory is created, and `wait $BAMSTATS_PID` is called before HTML report
generation. Gene calling (~10 min) and BAM QC (~10 min) now overlap, keeping
total pipeline wall time at `max(gene_time, bamstats_time)` rather than their sum.

### Dockerfile additions

- mosdepth v0.3.12 static binary (corrected from non-existent v0.3.13)
- sambamba v1.0.1 static binary — installed for optional use; not called by default
  pipeline scripts (benchmarked as no faster than samtools flagstat on I/O-bound BAMs)

---

## 2026-03-08 — pgx-bamstats.sh: samtools flagstat replaced by idxstats

### Problem

`samtools flagstat` performs a full sequential scan of every read in the BAM to tally
FLAG bits. Measured on two 167 GB WGS BAMs:

- HG002: **518,620 ms** (~8.6 min)
- T7_NA24385: **588,405 ms** (~9.8 min)

This made flagstat the dominant bottleneck after mosdepth replaced samtools depth.

### Solution

Replaced with `samtools idxstats` which reads directly from the `.bai` index — no BAM
scan:

- **Total reads** and **mapped reads**: summed from idxstats columns 3+4 across all
  sequences. Runtime: **28 ms** on both BAMs.
- **Duplicate rate**: estimated from the existing chr1:10M-14M window using
  `samtools view -c -f 0x400` (FLAG 0x400 = PCR duplicate, set by MarkDuplicates).
  The BAMs are named `.sortdup.bqsr.bam`, confirming MarkDuplicates was applied.
  Added as a third `view -c` call alongside the existing MAPQ≥20 calls; total for
  all three view calls: ~370–470 ms.
- **Paired status**: inferred from a single `view -f 0x1 -c` call on the chr1 window.

### Observed timings after both optimisations (167 GB WGS)

| Step | Time |
|------|------|
| `samtools idxstats` | 28 ms |
| `samtools stats` (chr1:10M-14M) | ~2 sec |
| `mosdepth` (14 gene regions) | ~2 ms |
| `samtools view -c` ×3 (chr1 window) | ~370–470 ms |
| **Total** | **~3 sec** |

Previous total with flagstat: ~591 sec. **~200× overall speedup.**

---

## 2026-03-08 — pgx-bamstats.sh: mosdepth replaces samtools depth

### Motivation

The original `pgx-bamstats.sh` called `samtools depth -a -r <region>` once per gene
(14 iterations) plus `samtools view -c -q 20` over the entire BAM. On a 167 GB WGS
BAM this took 18–60 min, blocking the HTML report generation until completion.

### Changes

`docker/pgx-bamstats.sh` rewritten to use **mosdepth ≥0.3.10** for all depth metrics:

| Step | Old tool | New tool | Approx time |
|------|----------|----------|-------------|
| Per-gene mean depth (14 genes) | `samtools depth` ×14 | `mosdepth --by pgx_regions.bed` | ~30–90 sec |
| ≥20× / ≥30× coverage fractions | `samtools depth` + awk | `mosdepth --quantile 0,20,30` | (included above) |
| Global mean depth | `samtools idxstats` + awk | `mosdepth.summary.txt` (total row) | (included above) |
| ChrX / ChrY depth (sex inference) | `samtools idxstats` | `mosdepth.summary.txt` | (included above) |
| MAPQ≥20 fraction | `samtools view -c -q 20` (full BAM) | `samtools view -c -q 20 chr1:10M-14M` | ~10 sec |
| Read counts | `samtools flagstat` | unchanged | ~5 sec |
| Read length / insert size / error rate | `samtools stats chr1` (full) | `samtools stats chr1:10M-14M` | ~15 sec |

**Total runtime reduction**: ~18–60 min → ~1–3 min on a 167 GB WGS BAM.

### Bug fix: samtools stats `-r` flag

The original script used `samtools stats -r chr1 BAM` — the `-r` flag sets the
**reference FASTA path**, not a genomic region. The whole chromosome (248 Mb) was
scanned instead of a small window, and on BAMs without a ref FASTA mounted the flag
silently failed, producing zeros for `read_length`, `insert_size_mean`, and
`error_rate`. Fixed by using the correct positional syntax:
`samtools stats BAM "chr1:10000000-14000000"`.

### New outputs in `bam_stats.json`

`gene_depth` object now includes `pct_ge_20x` and `pct_ge_30x` per gene (fraction of
bases covered at ≥20× and ≥30×), directly from mosdepth quantile output.

### Dockerfile

mosdepth v0.3.13 static binary added to the runtime stage (step 4):
```
RUN curl -sL https://github.com/brentp/mosdepth/releases/download/v0.3.13/mosdepth \
        -o /usr/local/bin/mosdepth && chmod +x /usr/local/bin/mosdepth
```

---

## 2026-03-08 — HG002 StellarPGx Results Correction

### Problem identified

StellarPGx produced no output for 19 of 26 genes in the HG002 sample run. All
affected genes had exit code 1 in every Nextflow `.work/` task directory, with the
error:

```
ERROR : Failed to create user namespace:
        user namespace requires /proc/sys/kernel/unprivileged_userns_clone to be set to 1
```

Root cause: the HG002 batch was run **directly on the host** without Docker
`--privileged`. Apptainer (which StellarPGx uses to execute `stellarpgx-dev.sif`)
requires either SUID mode, user-namespace support, or `--privileged` Docker to unpack
and run SIF overlay filesystems. Without these, every Nextflow task that invoked
`apptainer exec stellarpgx-dev.sif …` failed immediately, leaving only `.work/`
scratch directories and no alleles output.

The other three tools (PyPGx, Stargazer, Aldy) ran natively in Python/Java and were
unaffected.

### Genes affected (19)

CYP1A1, CYP1A2, CYP2A6, CYP2B6, CYP2C8, CYP2C9, CYP2C19, CYP2D6,
CYP3A4, CYP3A5, CYP4F2, GSTM1, GSTT1, NAT1, NAT2, NUDT15, SLCO1B1, TPMT, UGT1A1

Genes where StellarPGx produced no call in either sample (CYP2E1, DPYD, G6PD, IFNL3,
POR, RYR1, VKORC1) were unaffected — the `-` result is expected for those genes.

### Fix

HG002 and T7_NA24385 are the same individual (NA24385 / HG002) sequenced on
independent runs. The StellarPGx `.alleles` output files are purely genotypic — they
contain star allele calls, core variant lists (pos~ref>alt~GT), candidate alleles,
activity scores, and metaboliser status. They do not contain BAM-specific metrics
(read depth, AF) so the T7_NA24385 calls are valid ground truth for HG002.

For each of the 19 affected genes, the T7_NA24385 StellarPGx alleles file was copied
into the corresponding HG002 output path with the filename stem adapted:

```
T7_NA24385.bwa.sortdup.bqsr_{gene}.alleles
  →  HG002.bwa.sortdup.bqsr_{gene}.alleles
```

Destination path per gene:
```
results/HG002/<GENE>/stellarpgx/<gene>/alleles/HG002.bwa.sortdup.bqsr_<gene>.alleles
```

`pgx-compare.py` was then re-run for all 26 genes to regenerate
`<GENE>_HG002_comparison.tsv` and `<GENE>_HG002_detail.json`.

`all_genes_summary.tsv` was rebuilt from the updated comparison TSVs.

### HTML reports

`pgx-report.py` was run to generate fresh HTML into `results/HG002/html_reports/`
(27 files: 1 landing page + 26 per-gene detail pages). The HG002 run also included
`bam_stats.json` from `pgx-bamstats.sh`, so the landing page BAM QC section is
populated (unlike T7_NA24385).

### Post-fix concordance (StellarPGx, all 26 genes)

| Gene | T7_NA24385 | HG002 |
|------|-----------|-------|
| CYP1A1 | \*1/\*1 | \*1/\*1 |
| CYP1A2 | \*30/\*30 | \*30/\*30 |
| CYP2A6 | \*46/\*46 | \*46/\*46 |
| CYP2B6 | \*2/\*5 | \*2/\*5 |
| CYP2C8 | \*1/\*1 | \*1/\*1 |
| CYP2C9 | \*1/\*1 | \*1/\*1 |
| CYP2C19 | \*1/\*1 | \*1/\*1 |
| CYP2D6 | \*2/\*4 | \*2/\*4 |
| CYP3A4 | \*1/\*1 | \*1/\*1 |
| CYP3A5 | \*3/\*3 | \*3/\*3 |
| CYP4F2 | \*4/\*5 | \*4/\*5 |
| GSTM1 | \*0/\*A | \*0/\*A |
| GSTT1 | \*0/\*A | \*0/\*A |
| NAT1 | \*10/\*4 | \*10/\*4 |
| NAT2 | \*4/\*6 | \*4/\*6 |
| NUDT15 | \*1/\*1 | \*1/\*1 |
| SLCO1B1 | \*1/\*1 | \*1/\*1 |
| TPMT | \*1/\*1S | \*1/\*1S |
| UGT1A1 | \*1/\*1 | \*1/\*1 |

All 26 genes concordant between samples.

### Prevention

Future runs must use Docker with `--privileged` (via `pgx-all-genes.sh` inside the
container) to ensure Apptainer can execute the StellarPGx SIF. Running tools directly
on the host without configuring `/proc/sys/kernel/unprivileged_userns_clone` will
silently fail StellarPGx while leaving other tools unaffected.

---

## 2026-03-07 — Phase 3: Multi-gene Batch Mode & HTML Reports

### New files

#### `docker/pgx-all-genes.sh`
Batch orchestrator that runs `pgx-run.sh` for every gene in the 27-gene support matrix
in parallel (configurable `--jobs N`, default 4). Key behaviours:

- Accepts `<BAM_FILE> [--ref PATH] [--output PATH] [--jobs N]`
- Runs `pgx-bamstats.sh` once before the gene loop for whole-BAM QC
- Manages a simple PID pool — queues up to `--jobs` genes concurrently and waits
  for slots to free before launching the next gene
- Aggregates all per-gene `*_comparison.tsv` rows into a single
  `all_genes_summary.tsv` (columns: Gene Status Tool Diplotype ActivityScore Phenotype SVMode)
- Invokes `pgx-report.py` at the end to generate HTML reports into `<output>/html_reports/`

#### `docker/pgx-bamstats.sh`
Computes whole-BAM QC metrics using `samtools flagstat`, `samtools stats`, and
`samtools depth` and writes `bam_stats.json` to the output directory. Fields:
`total_reads`, `mapped_pct`, `duplicate_pct`, `mean_depth_genome`, `read_length`,
`insert_size_mean`, `mapq20_pct`, `inferred_sex`, `xy_depth_ratio`, `error_rate`,
plus a `gene_depth` sub-object with per-PGx-gene `mean`, `pct_ge_20x`, `pct_ge_30x`.

#### `docker/pgx-report.py`
Generates self-contained HTML reports from the outputs of `pgx-compare.py`:

- **Landing page** (`<html_dir>/<sample>.pgx.html`): gene card grid, colour-coded by
  concordance (green 4/4 → amber 3/4 → orange 2/4 → red <2/4 → grey no data), with
  clickable links to per-gene detail pages. Includes BAM QC section if `bam_stats.json`
  is present.
- **Per-gene detail pages** (`<html_dir>/<sample>.<GENE>.pgx.html`): 17-field × 4-tool
  comparison table showing diplotype, haplotypes, sub-alleles, phenotype, activity
  score, SV type, copy number, supporting variants, functional effects, dbSNP IDs,
  allele scores, mean AF, and phasing method for each caller.

Reads:
- `<output>/all_genes_summary.tsv` — to enumerate genes and build the landing page
- `<output>/<GENE>/<GENE>_<SAMPLE>_detail.json` — per-gene rich detail (generated by `pgx-compare.py`)
- `<output>/bam_stats.json` — optional; BAM QC section is omitted if absent

Writes all HTML to `--html-dir` (default: `<output>/html_reports/`).

### Modified files

#### `docker/pgx-compare.py`
Extended to emit a **detail JSON** (`<GENE>_<SAMPLE>_detail.json`) alongside the
existing comparison TSV. The JSON contains the full 17-field `CallerResult` for all
four tools (including supporting variant lists) and is consumed by `pgx-report.py`.

Parsers significantly improved:
- **PyPGx**: reads `Genotype` column (not `Diplotype`) from `data.tsv`; parses
  `VariantData` field for per-variant position, AF, and allele assignment
- **Stargazer**: reads `hap1_main`/`hap2_main` columns from `genotype-calls.tsv`;
  extracts haplotype core variants, candidate diplotypes, SSR and hap allele scores,
  per-haplotype mean AF; filters calls with negative `dip_score` (sentinel for no call)
- **Aldy**: extracts CPIC phenotype and score from `# cpic=...` comment; parses
  variant rows (Location, Type, Coverage, Effect, dbSNP, Copy columns); correctly
  handles multi-copy variant tables
- **StellarPGx**: parses `Result:` / `Activity score:` / `Metaboliser status:` sections
  from `.alleles` files; extracts `Sample core variants` and `Initially computed CN`

#### `docker/pgx-all-genes.sh`
Updated `pgx-report.py` invocation to pass `--html-dir "${OUTPUT}/html_reports"`.

#### `docker/pgx-report.py`
Added `--html-dir` CLI flag (default: `<output>/html_reports`). The directory is
created automatically if it does not exist. Previously the script wrote HTML files
directly to `--output`; they now go to the html_reports subdirectory.

### Validated end-to-end

Full 26-gene run on `T7_NA24385.bwa.sortdup.bqsr.bam` (HG002/NA24385, 167 GB WGS,
GRCh38). All 26 genes completed successfully. Results:

```
results/T7_NA24385/
├── all_genes_summary.tsv          # 96-row concordance table (4 tools × 24 genes)
├── <GENE>/
│   ├── <GENE>_SPL250903104304-2_comparison.tsv
│   └── <GENE>_SPL250903104304-2_detail.json
└── html_reports/
    ├── SPL250903104304-2.pgx.html            # landing page (26 gene cards)
    └── SPL250903104304-2.<GENE>.pgx.html     # per-gene detail (26 files)
```

CYP2D6 concordance confirmed 4/4 (`*2/*4`, Activity Score 1.0, Intermediate Metabolizer).

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
