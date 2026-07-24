# WES (Whole-Exome Capture) Limitations

**Scope:** whether PGx Suite produces trustworthy results on a **whole-exome
capture (WES) BAM/CRAM** rather than the WGS BAM it is designed for. This is an
assessment document — it describes behaviour, it does not change pipeline logic.

## TL;DR

PGx Suite is the **most WES-viable** of the OmniGen suites. Its architecture is
*fail-safe by construction*: every gene either produces a call backed by covered
reads, or it reports **`NO_CALL`** — it never manufactures a wild-type
(`*1/*1`, "normal metabolizer") diplotype from missing data.

- **VALID on WES** — the SNV/genotype-based callers, *for any gene whose
  star-allele-defining exonic positions are on the capture kit*: PharmCAT,
  VCF-Check, and the SNV component of PyPGx / Stargazer / Aldy (the per-gene
  `bcftools` VCF feeds them).
- **NO_CALL on WES (fail-safe, never wrong)** — the copy-number / structural
  methods: Cyrius (CYP2D6 CN/hybrid), PyPGx-SV and Stargazer-SV CN genes,
  StellarPGx (WGS-only by design). Exome's non-uniform capture depth breaks the
  depth-ratio CN model; the coverage floor and ALT-aware gate convert that into
  `NO_CALL`, not a false deletion.
- **Capture-dependent** — HLA (OptiType class I, arcasHLA class II) works *iff*
  the MHC exons are captured (most exome kits do); mutserve MT-RNR1 needs chrM
  reads (most exome kits **do not** capture mtDNA) → typically `NO_CALL`.
- **The one caveat:** a gene **absent from the capture kit** falls below the
  depth floor and is reported `NO_CALL`. `NO_CALL` means *"no data — do not
  prescribe on this"*. It must **never** be read as "normal metabolizer". This
  is the single most important WES misread to guard against, and the pipeline's
  coverage floor is exactly what enforces it (see below).

## Why it fails safe — the two gates (evidence)

### 1. Coverage floor → `NO_CALL`, never wild-type

`docker/pgx-compare.py:44`
```
DEFAULT_MIN_DEPTH = 10.0
```

`docker/pgx-compare.py:1664-1669` — the gate:
```python
# Gate the call when (a) coverage is too low, or (b) alt_block is set …
low_cov = region_depth is not None and region_depth < min_depth
no_call = low_cov or bool(alt_block)
verdict = compute_verdict(results, no_call, gene)
```

`docker/pgx-compare.py:1556-1558` — `compute_verdict` when gated returns a hard
`NO_CALL`, agreement counts zeroed:
```python
if no_call:
    return {"consensus_diplotype": "NO_CALL", "consensus_phenotype": "-",
            "n_agree": 0, "n_called": 0, "status": "no_call", "reconciled": False}
```

`docker/pgx-compare.py:1719-1724` — **every tool's diplotype is suppressed** in
the TSV so neither the per-gene output nor the batch headline can surface an
unreliable call; the reason recorded is `Insufficient coverage (Nx<10x)`.

Region mean depth is measured directly from the BAM per gene
(`docker/Snakefile:383`, `samtools coverage -r {region} …`), using the region
column of `docker/genes.tsv`. **Consequence for WES:** a gene not on the kit has
near-zero depth → `NO_CALL` (correct, fail-safe). Note also a *conservative
false-negative*: because the floor uses the **whole-gene** mean (exon + intron),
a captured gene whose intronic zeros drag the mean under 10× can be `NO_CALL`
even though its exons are well covered — safe, but it lowers WES yield. Override
per run with `--min-depth` (`docker/pgx-compare.py:1770-1772`) if assay
validation supports it.

### 2. ALT-aware / paralog-SV gate → `NO_CALL`

`docker/pgx_altcheck.py:4-7` and `docker/pgx-compare.py:1818-1824`: a
non-ALT-aware alignment (or a paralog/SV gene from data that cannot resolve the
paralog) sets `alt_block`, which forces the same `NO_CALL` path — e.g. it
prevents a false CYP2D6 `*5/*5` deletion. This is independent of WES but
compounds it: exome capture of paralogous CYP genes is not a reliable CN
substrate.

## Per-caller / per-gene-class verdict

| Caller / method | Model | WES verdict | Evidence |
|---|---|---|---|
| **bcftools** (per-gene VCF) | SNV/indel genotyping over the gene region | **VALID** for captured exons — supplies the VCF the star-allele callers vote on | `docker/Snakefile` standard-gene rules; `region()` from `genes.tsv` |
| **PharmCAT** | CPIC reference named-allele matcher; **authority** for UGT1A1, CYP2B6, CYP4F2 | **VALID** if the defining positions are captured. Caveat: UGT1A1 `*28` is a TA-repeat *promoter* indel that some exome kits miss | `AUTHORITATIVE` `docker/pgx-compare.py:1466-1474` |
| **VCF-Check** (in-house) | Direct genotype lookup vs CPIC table; **authority** for RYR1, CACNA1S, G6PD, VKORC1, IFNL3, CYP2C-CLUSTER, ABO | **VALID** if the single position is captured — the most WES-robust method (one SNP, no depth model) | `AUTHORITATIVE` `docker/pgx-compare.py:1466-1471`; `docker/vcf_check_variants.json` |
| **PyPGx** (SNV part) | Beagle-phased star alleles | **VALID** for captured SNV-defined alleles | `docs/ToolsDocumentation.md:34` (WES = Yes) |
| **Stargazer** (SNV part) | Beagle-phased star alleles | **VALID** for captured SNV-defined alleles | `docs/ToolsDocumentation.md:34` (WES = Yes) |
| **Aldy** (SNV part) | ILP over genotypes | **VALID** for SNVs; **CN unreliable** on exome | `docs/ToolsDocumentation.md:34,271-272` |
| **StellarPGx** | Genome-graph (graphtyper) | **NO_CALL / breaks** — WGS-only by design; on exome yields no/failed output (non-fatal, does not abort the gene) | `docs/ToolsDocumentation.md:34,345-346` |
| **Cyrius** | Targeted CYP2D6 CN/paralog/hybrid, depth-based; **authority** for CYP2D6 | **NO_CALL** on exome — non-uniform capture depth breaks CN; declines → consensus fallback → coverage/ALT gate → `NO_CALL`, never a false `*5` | `AUTHORITATIVE` `docker/pgx-compare.py:1470`; ALT gate `docker/pgx_altcheck.py:4-7` |
| **PyPGx-SV / Stargazer-SV** (CN genes) | Depth-of-coverage + control-region CN normalisation | **NO_CALL** for the CN component on exome. Affected genes below | `docker/pgx-compare.py:1480-1483` (`_sv_note`) |
| **OptiType** (HLA class I A/B/C) | ILP on MHC-filtered reads | **Capture-dependent** — valid iff HLA exons captured (most exome kits do) | `docs/ToolsDocumentation.md:34` ("Yes if MHC captured") |
| **arcasHLA** (HLA class II DQA1/DQB1/DRB1) | Read-partitioned HLA typing | **Capture-dependent** — valid iff HLA exons captured | `README.md` Bundled Tools; `genes.tsv` HLA-DQA1/DQB1/DRB1 |
| **mutserve** (MT-RNR1) | chrM allele-fraction pileup | **NO_CALL** typically — most exome kits do not capture mtDNA | `docs/ToolsDocumentation.md:34` ("only if chrM captured") |

### CN/SV gene lists (the genes whose *structural* call is unreliable on WES)

From `docker/genes.tsv` (`pypgx_sv` / `stargazer_sv` flags, and the `ALT_CONTIG`
region):

- **`pypgx_sv=1`** — CYP2A6, CYP2B6, CYP2D6, CYP2E1, CYP4F2, G6PD, GSTM1, GSTT1,
  SLC22A2, UGT2B15
- **`stargazer_sv=1`** — CYP2A6, CYP2B6, CYP2D6
- **`ALT_CONTIG`** — GSTT1 (whole-gene deletion detection needs ALT-contig
  reads; unreliable on exome)

For these genes the **SNV-defined** sub-alleles may still be called from captured
exons; only the **copy-number / gene-deletion** determination degrades to
`NO_CALL`.

## Does PGx need a target BED to run on WES?

**No — a target BED is not required to run.** The per-gene coordinates come from
the `region` column of `docker/genes.tsv` (`region()` in `docker/Snakefile:110,363`),
and depth is measured with `samtools coverage -r {region}` (`docker/Snakefile:383`).

`docker/pgx_cram_regions.bed` is a **CRAM-extraction convenience only**: it slices
a CRAM down to the PGx loci into a region BAM
(`samtools view -b -T ref.fa -L pgx_cram_regions.bed input.cram`), documented at
`docker/Snakefile:20-22`. With a **BAM** input, `depth_bam` defaults to the input
BAM and no BED is consulted (`docker/Snakefile` config: `DEPTH_BAM = config.get("depth_bam", BAM)`).
It is *not* the WES capture kit's BED and is not used to define what "captured"
means — WES capture is an intrinsic property of the input reads, which the
coverage floor detects per gene.

## Bottom line

- **Trustworthy on WES (captured genes):** the CPIC-critical SNV-driven genes via
  their authorities — **PharmCAT** (UGT1A1*, CYP2B6, CYP4F2) and **VCF-Check**
  (RYR1, CACNA1S, G6PD, VKORC1, IFNL3, CYP2C-CLUSTER), plus SNV star alleles from
  PyPGx / Stargazer / Aldy for genes like CYP2C9, CYP2C19, TPMT, NUDT15, SLCO1B1.
- **`NO_CALL` on WES (correctly, never wrong):** CYP2D6 (Cyrius CN), the CN/SV
  genes above, StellarPGx contributions, and mtDNA MT-RNR1 if chrM isn't captured.
- **Capture-dependent:** HLA class I/II (usually captured).
- **The caveat, restated:** a gene **not on the capture kit** → depth floor →
  **`NO_CALL`**. `NO_CALL` is *absence of evidence*, not evidence of a normal
  metabolizer. Any downstream report or clinician consuming a WES run must treat
  `NO_CALL` as "not assessed", never as `*1/*1`.

*(`\*` UGT1A1 `*28`/`*80` is a promoter TA-repeat indel; capture varies by kit —
treat a UGT1A1 call on exome as kit-dependent.)*
