# pgx-suite — Snakemake Migration Plan

**Date:** 2026-06-17
**Decision:** migrate orchestration from bash (`pgx-run.sh` + `pgx-all-genes.sh`) to Snakemake. Chosen over "harden bash" to structurally resolve config duplication, resource governance, Nextflow-per-gene overhead, resume, and provenance in one move.
**Companion docs:** `2026-06-17-pgx-suite-eng-review.md` (findings), `2026-06-17-pgx-suite-architecture-critique.md`.

---

## Guiding principles

- **Strangler fig, not big bang** (Fowler). The bash pipeline stays runnable until Snakemake reproduces the validated HG002 result exactly at each stage. No stage is "done" until its output matches the current pipeline on the golden sample.
- **Make the change easy, then make the easy change** (Beck). Land the single-source `genes.tsv` FIRST — it's the input both the bash fixes and Snakemake need. That refactor is not throwaway; it's the shared foundation.
- **Don't let a clinical bug wait for a multi-week migration.** The P0 fixes (especially A1, the HTML gate-bypass) ship on the *current* pipeline in Stage 0.
- **Smallest blast radius:** Snakemake runs *inside the existing Docker image*, replacing only the two bash scripts. Callers, `pgx-compare.py`, `pgx-report.py`, `pgx-bamstats.sh`, the image contents, and the `--privileged`/Apptainer/StellarPGx reality are all unchanged. (StellarPGx's `--privileged` dependency is a separate decision — out of scope here.)

## What changes vs what stays

```
REPLACED                          KEPT (wrapped as Snakemake rules)
─────────                         ─────────────────────────────────
pgx-run.sh        ──►  Snakefile   bcftools / Aldy / StellarPGx / PyPGx / Stargazer commands
pgx-all-genes.sh  ──►  rules +     pgx-compare.py   (parsers + unified verdict)
 (job pool)            config       pgx-report.py    (HTML)
hardcoded arrays  ──►  genes.tsv    pgx-bamstats.sh  (becomes a rule)
                       + config.yaml CRAM extraction  (becomes a rule)
```

## Target DAG (per sample)

```
                         config.yaml + genes.tsv  (single source of truth)
                                     │
              ┌──────────────────────┼───────────────────────┐
              ▼                                               ▼
   rule extract_regions (CRAM→pgx_input.bam, once)   rule bamstats (on _ORIG_BAM)
              │                                               │
   ── per {gene} wildcard ──────────────────────────         │
   ▼            ▼            ▼            ▼                    │
 bcftools     aldy      stellarpgx   pypgx_prep (SV genes)    │
   │            │            │            │                   │
   └──► VCF     │            │       pypgx_pipeline           │
        │       │            │       stargazer_genotype       │
        └───────┴────────────┴────► compare {gene} ──► detail.json / comparison.tsv
                                                  │            │
   special rules: hla (OptiType) · mt (mutserve) · gstt1 (alt-contig)
                                                  │            │
                                                  └──► report ◄┘  (reads coverage.verdict)
```

Special cases become **conditional rules / input functions** keyed off `genes.tsv` columns (tool-support flags, SV flag, region), so the bash bypass branches disappear into declarative config.

## Stages (each independently shippable + validated against current pipeline)

### Stage 0 — Foundation on the CURRENT pipeline (do first, no Snakemake yet)
This is shared groundwork + the urgent fixes. Ships value immediately and de-risks everything after.
- **A1 [P0, clinical]** — `pgx-report.py` honors `detail.json → coverage.verdict == NO_CALL`; short-circuit gene card/detail before `compute_concordance`.
- **A2 [P0]** — make `pgx-compare.py` write the final verdict (consensus / DISCORDANT / NO_CALL, n_agree/n_called) into `detail.json`; `pgx-report.py` and the batch summary **read** it instead of recomputing. Kills the duplicate tie logic.
- **C1 [P0]** — create `genes.tsv` (gene, region, per-tool support, SV flag, CRAM region). Refactor `pgx-run.sh`, `pgx-compare.py`, and a BED generator to read it. This file becomes Snakemake's config.
- **T1 [P0]** — `test_parsers.py` + one real output fixture per tool under `docker/testdata/`.
- Validate: full bash pipeline on HG002 still yields the known 4/4 CYP2D6 `*2/*4`, and an induced uncovered gene now shows NO_CALL in the HTML.

### Stage 1 — Snakefile for one gene, one sample
- Wrap existing per-tool commands as rules; wildcards `{sample}`, `{gene}`; read `genes.tsv`.
- Scope: CYP2D6 only, parallel within the gene.
- Validate: Snakemake CYP2D6 output == Stage 0 bash CYP2D6 output, byte-for-byte on the calls.

### Stage 2 — All genes via Snakemake scheduler
- Expand wildcards to all genes; delete the hand-rolled job pool. Use `--cores`/`resources: mem_mb` so peak load is governed (fixes P1).
- Validate: all-genes summary == current `all_genes_summary.tsv` on HG002.

### Stage 3 — Special cases + QC + report as rules
- CRAM extraction, HLA (OptiType), MT-RNR1 (mutserve), alt-contig GSTT1, bamstats, and report all become rules with proper input/output deps.
- Validate: end-to-end HTML report == current, on both a BAM and a CRAM sample.
- **Decommission `pgx-run.sh` + `pgx-all-genes.sh`** only after this passes.

### Stage 4 — Free wins the migration unlocks
- Provenance: per-rule tool-version capture + reference build + pipeline git SHA into `detail.json` and Snakemake `--report` (fixes the provenance gap).
- Resume: rely on Snakemake output-existence + `--rerun-incomplete` (fixes P3). Remove Stargazer's `rmtree`-on-start now that rules own their outputs.
- StellarPGx: evaluate one Nextflow invocation across genes, or replace with a direct graphtyper rule, to kill the 25× cold start (P2).

## Risks & mitigations

- **Behavior drift.** Mitigation: golden-output validation gate at every stage; the bash pipeline remains the oracle until Stage 3.
- **Snakemake learning curve = innovation token.** Mitigation: it's Python-native (matches the codebase) and the rules are thin wrappers around commands that already work.
- **Special-case rules (HLA/MT/GSTT1) are where DAG bugs hide.** Mitigation: they land last (Stage 3), each validated independently.
- **Scope creep into the StellarPGx/`--privileged` problem.** Explicitly out of scope; tracked separately.

## Effort (rough)

- Stage 0: ~1.5 days human / ~2 hr CC. Ships the clinical fix + foundation.
- Stages 1-3: ~3-5 days human / ~half day CC, dominated by validation, not code.
- Stage 4: ~1-2 days human / ~1-2 hr CC.

## NOT in scope

- StellarPGx `--privileged`/Apptainer removal (separate architectural decision).
- Any change to the callers, parser logic (beyond A2's verdict write), or clinical content of the report.
- Multi-sample batch UX beyond what wildcards give for free.
