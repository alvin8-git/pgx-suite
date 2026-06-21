# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A single Docker image that runs **multiple pharmacogenomics (PGx) star-allele callers** over one WGS BAM/CRAM and reconciles their calls into a clinical HTML report. **GRCh38 only.** All caller orchestration lives in `docker/*.sh` and `docker/*.py`; the rest of the repo is vendored caller source + reference data.

Callers: **PyPGx** 0.26.0, **Stargazer** 2.0.3, **Aldy** 4.8.3, **StellarPGx** 1.2.7 (Nextflow+Apptainer), plus **OptiType** (HLA-A/B, Apptainer SIF) and **mutserve** (MT-RNR1, baked-in jar). Plus three **verdict-authority** methods: **Cyrius** (CYP2D6 SV/CNV; vendored at `/opt/cyrius`, PolyForm Strict — image not redistributable), **PharmCAT** 3.2.0 (UGT1A1/CYP2B6/CYP4F2; official `pgkb/pharmcat` Apptainer SIF, runs once per sample), and in-house **VCF-Check** (RYR1/VKORC1/IFNL3/G6PD/CACNA1S; CPIC variant table in `docker/vcf_check_variants.json`).

## Build & run

```bash
# Build the image (multi-stage; ~15-25 min first time)
docker build -t pgx-suite:latest .

# Smoke tests — no BAM/volumes needed (PyPGx/Stargazer/Aldy):
docker run --rm pgx-suite:latest bash /opt/pgx/test.sh
# Full smoke tests incl. StellarPGx (needs volume mounts):
./docker/docker-run.sh -- bash /opt/pgx/test.sh

# Full pipeline, all genes, one sample — host-side launcher (handles all mounts):
./run_pgx_suite.sh /path/to/sample.bam [--ref REF.fa] [--output DIR] [--jobs N]
```

Inside the container the orchestrator is **Snakemake** (`/opt/pgx/Snakefile`):
- All genes: `snakemake -s /opt/pgx/Snakefile --cores N --config bam=/pgx/data/<sample>.bam outdir=/pgx/results`
- One gene (debug): add `genes="CYP2D6"` to `--config` (space-separated for several).
- CRAM input: also pass `depth_bam=<extracted-region-bam>` (see the `depth_bam` note in the Snakefile).

Snakemake replaced the former bash orchestrators `pgx-run.sh` / `pgx-all-genes.sh` (removed 2026-06-17 after full-set 31-gene equivalence on HG002). Unit checks live in `docker/test_*.py` (`python3 docker/test_parsers.py`, `test_verdict.py`, `test_coverage_gate.py`, `test_genes_config.py`); whole-pipeline validation is by running real samples (see `docs/AxiomValidation.md`).

## Architecture & data flow

`run_pgx_suite.sh` (host) → docker run → `snakemake -s /opt/pgx/Snakefile` (container) → per-gene caller rules → `pgx-compare.py` per gene → `pgx-bamstats.sh` + `pgx-report.py` (final HTML).

**`docker/Snakefile` is the orchestration core** — a wildcard `{gene}` DAG. Everything is driven by **`docker/genes.tsv`**, the single source of truth for per-gene region, tool support (pypgx/stargazer/aldy/stellarpgx/optitype/mutserve/vcf_check/cyrius/pharmcat), SV flags, and Stargazer control. **To add/change gene support, edit `genes.tsv`** — both the Snakefile and `pgx-compare.py` read it (guarded by `test_genes_config.py`). Each caller is a rule whose declared output is a `logs/<tool>.status` sentinel holding the real exit code, so a failing caller is **non-fatal** (the rule still succeeds); the `compare` rule reads the sentinels to pass `--failed-tools`. `--cores` budgets the genes×tools scheduling. HLA-A/B route to the `optitype` rule, MT-RNR1 to `mutserve`, and alt-contig GSTT1 runs pypgx(SV, no VCF)+stellarpgx with no bcftools — all selected by `genes.tsv` flags.

**`pgx-compare.py`** has one `parse_<tool>()` per caller (each writes a different format/location). It is the **single verdict authority** and reconciles in ordered tiers: (1) coverage gate (NO_CALL below `--min-depth`); (2) **variant-tier reconciliation** — `reconcile.py` + `allele_synonyms.json` synonym-collapse diplotypes before counting agreement (DPYD `*9A`=`*S10`=`c.85T>C`=`rs1801265`); (3) **authoritative callers** — the `AUTHORITATIVE` dict makes Cyrius the verdict for CYP2D6, PharmCAT for UGT1A1/CYP2B6/CYP4F2, VCF-Check for RYR1/VKORC1/IFNL3/G6PD/CACNA1S when that method called; (4) **phenotype tier** — `_phenotype_consensus()` upgrades a diplotype-string discordant/majority to concordant when ≥2 phenotype-emitting callers agree on the CPIC phenotype (`Indeterminate`=abstain), tagged `authority:"Phenotype"`. Emits `<GENE>_<SAMPLE>_comparison.tsv` + `<GENE>_<SAMPLE>_detail.json` — the latter's `verdict` block (status + `authority` + `phenotype_tier`) is **read** by `pgx-report.py` and the summary, never recomputed. The report colours cards/ticker by verdict **status** (not raw tool fraction) and shows an authority badge.

**`pgx-report.py`** reads the per-gene `*_detail.json` files + `bam_stats.json` (from `pgx-bamstats.sh`) and produces `<SAMPLE>_pgx_report.html` (landing page + per-gene 17-field detail tables). It is the largest file; CAP/CLIA clinical-report fields live here.

**`pgx-bamstats.sh`** is the `bamstats` rule (I/O-bound, ~6 min on 167 GB WGS); the `report` rule depends on its `bam_stats.json`. For CRAM inputs it takes a separate depth BAM (the Snakefile's `depth_bam` config) because mosdepth on full CRAM is too slow.

## Submodules & vendored code

- `pypgx/` and `StellarPGx/` are **git submodules**. `stargazer-grc38-2.0.3/` is vendored source (non-commercial UW licence — do **not** push the image to a public registry).
- **pypgx has local patches** (pandas 2.x fixes) committed inside the submodule. To change pypgx: commit inside `pypgx/`, then bump the submodule pointer in the parent repo.

## Conventions

- The Dockerfile is multi-stage: stage 1 compiles the Python venv, stage 2 is a lean runtime. scikit-learn is intentionally omitted (only pypgx's unused `train_cnv_caller` needs it).
- StellarPGx requires `--res_init/--db_init/--caller_init` flags because its `nextflow.config` uses `$PWD`-relative paths; it uses gene id `cypor` not `por`.
- Aldy profile flag is `-p illumina` (not `hg38`).
- Stargazer GDF must be written **outside** the stargazer output dir — Stargazer `rmtree`s its output dir on startup.

---

# context-mode — MANDATORY routing rules

You have context-mode MCP tools available. These rules are NOT optional — they protect your context window from flooding. A single unrouted command can dump 56 KB into context and waste the entire session.

## BLOCKED commands — do NOT attempt these

### curl / wget — BLOCKED
Any Bash command containing `curl` or `wget` is intercepted and replaced with an error message. Do NOT retry.
Instead use:
- `ctx_fetch_and_index(url, source)` to fetch and index web pages
- `ctx_execute(language: "javascript", code: "const r = await fetch(...)")` to run HTTP calls in sandbox

### Inline HTTP — BLOCKED
Any Bash command containing `fetch('http`, `requests.get(`, `requests.post(`, `http.get(`, or `http.request(` is intercepted and replaced with an error message. Do NOT retry with Bash.
Instead use:
- `ctx_execute(language, code)` to run HTTP calls in sandbox — only stdout enters context

### WebFetch — BLOCKED
WebFetch calls are denied entirely. The URL is extracted and you are told to use `ctx_fetch_and_index` instead.
Instead use:
- `ctx_fetch_and_index(url, source)` then `ctx_search(queries)` to query the indexed content

## REDIRECTED tools — use sandbox equivalents

### Bash (>20 lines output)
Bash is ONLY for: `git`, `mkdir`, `rm`, `mv`, `cd`, `ls`, `npm install`, `pip install`, and other short-output commands.
For everything else, use:
- `ctx_batch_execute(commands, queries)` — run multiple commands + search in ONE call
- `ctx_execute(language: "shell", code: "...")` — run in sandbox, only stdout enters context

### Read (for analysis)
If you are reading a file to **Edit** it → Read is correct (Edit needs content in context).
If you are reading to **analyze, explore, or summarize** → use `ctx_execute_file(path, language, code)` instead. Only your printed summary enters context. The raw file content stays in the sandbox.

### Grep (large results)
Grep results can flood context. Use `ctx_execute(language: "shell", code: "grep ...")` to run searches in sandbox. Only your printed summary enters context.

## Tool selection hierarchy

1. **GATHER**: `ctx_batch_execute(commands, queries)` — Primary tool. Runs all commands, auto-indexes output, returns search results. ONE call replaces 30+ individual calls.
2. **FOLLOW-UP**: `ctx_search(queries: ["q1", "q2", ...])` — Query indexed content. Pass ALL questions as array in ONE call.
3. **PROCESSING**: `ctx_execute(language, code)` | `ctx_execute_file(path, language, code)` — Sandbox execution. Only stdout enters context.
4. **WEB**: `ctx_fetch_and_index(url, source)` then `ctx_search(queries)` — Fetch, chunk, index, query. Raw HTML never enters context.
5. **INDEX**: `ctx_index(content, source)` — Store content in FTS5 knowledge base for later search.

## Subagent routing

When spawning subagents (Agent/Task tool), the routing block is automatically injected into their prompt. Bash-type subagents are upgraded to general-purpose so they have access to MCP tools. You do NOT need to manually instruct subagents about context-mode.

## Output constraints

- Keep responses under 500 words.
- Write artifacts (code, configs, PRDs) to FILES — never return them as inline text. Return only: file path + 1-line description.
- When indexing content, use descriptive source labels so others can `ctx_search(source: "label")` later.

## ctx commands

| Command | Action |
|---------|--------|
| `ctx stats` | Call the `ctx_stats` MCP tool and display the full output verbatim |
| `ctx doctor` | Call the `ctx_doctor` MCP tool, run the returned shell command, display as checklist |
| `ctx upgrade` | Call the `ctx_upgrade` MCP tool, run the returned shell command, display as checklist |
