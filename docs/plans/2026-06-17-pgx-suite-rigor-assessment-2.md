# pgx-suite — Second Rigor Assessment (post-migration)

**Date:** 2026-06-17
**Scope:** find design flaws *remaining* after the Snakemake migration, verdict/coverage
gate, single-source `genes.tsv`, and report hardening already shipped. Engineering-infra
review via plan-eng-review posture (boring-by-default, blast radius, systems-over-heroes,
reversibility). Each finding verified against the current code.

**Framing.** The pipeline now produces correct, validated calls on *good* input (31-gene
HG002 equivalence; TTSH ~80% vs Axiom). The remaining flaws are the gap between "correct on
good input" and "trustworthy in a regulated lab": it cannot prove *what produced a result*,
it cannot tell a clean all-normal sample from a broken run, and a rebuild months from now
isn't guaranteed to reproduce. None of these are caught by the existing tests.

---

## CRITICAL

### F1 — A broken run still exits 0 (silent whole-run failure)
Every caller is non-fatal by design (`logs/<tool>.status` sentinels), and the Snakefile's
final target is the report. So a wrong reference, an unreadable BAM, or a missing index makes
**every** caller fail → the report renders as all-`NO_CALL`/`failed`, and the pipeline exits
**0**. An operator reading the exit code (or a batch wrapper) cannot distinguish "this sample
is genuinely all wild-type" from "the whole run was misconfigured." This is the
systems-over-heroes / 3am failure mode.
**Fix:** after `compare`, a sanity gate (a Snakemake rule or a step in `pgx-report`/launcher):
if the fraction of genes with `status in {no_data, no_call}` or all-callers-failed exceeds a
threshold (e.g. >60%), still write the report but **exit non-zero** with a one-line reason.
Cheap; turns a silent misconfig into a loud one.

### F2 — No provenance in `detail.json` or the report
`detail.json` carries `verdict`, `coverage`, and per-tool calls — but **not** tool versions,
the reference build/checksum, the pipeline commit, the run timestamp, or the command lines.
For any CAP/CLIA-reportable result you must be able to answer "which PyPGx, which reference,
which pipeline version produced this call?" Today you cannot, from the artifact alone. If
PyPGx 0.26→0.27 changes a call, nothing records which version ran.
**Fix:** write a `provenance` block once per run — pipeline git SHA (or a baked `VERSION`),
`pypgx/aldy/stargazer/stellarpgx/samtools/bcftools --version`, reference path + size (md5 if
affordable), and an ISO timestamp — into `detail.json` and a report footer. Low effort, high
regulatory value.

---

## MAJOR

### F3 — Snakemake is installed unpinned (reproducibility regression)
`Dockerfile:57` is `RUN pip install --no-cache-dir snakemake` — **no version**. The
orchestrator that the whole pipeline now depends on floats to latest; a rebuild after a
Snakemake major bump can break the Snakefile silently. (Contrast: `aldy==4.8.3`, `ortools`
pinned.) Boring-by-default + reversibility both say pin it.
**Fix:** pin `snakemake==9.x` (whatever the validated build used). One-line.

### F4 — No per-rule resource declarations → oversubscription risk
No rule declares `threads:` or `resources: mem_mb:`. `--cores N` therefore schedules N jobs
each *accounted* as 1 core, while the tools multithread internally (Nextflow JVM, bcftools,
Aldy ILP, PyPGx sklearn). On a shared node or small cloud VM, N genes × several real threads
each oversubscribes CPU and can OOM (StellarPGx/graphtyper and PyPGx are the memory-heavy
rules). This is the unfinished half of the earlier "resource governor" finding — the job pool
is gone, but the scheduler still can't see real usage.
**Fix:** declare realistic `threads`/`resources: mem_mb` per rule and run with
`--resources mem_mb=<budget>`. Then `--cores` actually bounds the machine.

### F5 — No CI runs the test suite
`test_parsers.py`, `test_verdict.py`, `test_coverage_gate.py`, `test_genes_config.py` exist
and are pure-Python (no Docker/data) — but **nothing runs them automatically**. A parser
format-drift, a verdict-logic regression, or a `genes.tsv` typo ships silently. "Well-tested
code is non-negotiable" is undercut if the tests only run when someone remembers.
**Fix:** a GitHub Actions workflow (`.github/workflows/test.yml`): `setup-python`, then
`python3 docker/test_*.py`, on push + PR. Near-zero cost, catches the exact regressions the
tests were written for.

### F6 — No input-integrity gate
The pipeline trusts the BAM and its index. A truncated BAM or a stale `.bai` produces
empty/partial reads → now masked as `NO_CALL` by the coverage gate, i.e. a corrupt input
reads as "low coverage," not "bad file."
**Fix:** `samtools quickcheck` on the input + an index-freshness check at launch
(`run_pgx_suite.sh` or a Snakemake input function); fail loudly before any caller runs.

### F7 — CPIC clinical database hardcoded inside `pgx-report.py`
`pgx-report.py:117` defines `CPIC_DB` (drug recommendations, phenotype text, tiers) inline —
clinical content fused into the presentation layer. Adding or updating a gene's clinical
interpretation means editing a 2000-line Python file. This is the same single-source smell C1
fixed for `genes.tsv`, still present for the clinical knowledge base.
**Fix:** extract `CPIC_DB` → `docker/cpic.json` (or `.yaml`), load it in `pgx-report.py`, and
guard it with a small schema test. Lets clinical staff update dosing text without touching code.

---

## MINOR

### F8 — Nextflow `.work` dirs are never cleaned
Each StellarPGx rule writes `OUT/genes/<gene>/stellarpgx/.work` and nothing removes it — ~15
Nextflow work trees of intermediate files per sample, accumulating across runs.
**Fix:** `rm -rf` the `.work` dir on rule success, or run Nextflow with cleanup enabled.

### F9 — genes.tsv coordinates aren't validated against the reference
A typo'd region or a contig-naming mismatch (`chr1` vs `1`) silently yields empty coverage →
**NO_CALL** masks it as low depth, not "bad coordinates."
**Fix:** a startup check that each gene's contig/region exists in `ref.fai`.

### F10 — Sample-name derivation has two sources
The Snakefile derives the sample from the BAM `@RG SM` tag; `run_pgx_suite.sh` and naming
elsewhere may use the filename. If they disagree, the report filename and internal sample ID
diverge.
**Fix:** derive once (prefer `@RG SM`) and thread it through, or assert they match.

---

## Recommendation

The pipeline is functionally complete; these close the "trustworthy artifact" gap. Suggested
order, smallest-blast-radius first:

1. **F2 provenance + F1 fail-loud gate** — the two clinical-trust P1s. Medium effort, high value.
2. **F3 pin Snakemake + F5 CI** — reproducibility + regression safety. Tiny effort.
3. **F6 input integrity + F4 resources** — robustness on bad input / shared hardware.
4. **F7 CPIC → data file** — the bigger refactor; extends single-source to clinical content.
5. **F8/F9/F10** — housekeeping.

F2, F3, F5 are the highest value-per-effort and I'd do them first.
