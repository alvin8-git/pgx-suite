# pgx-suite — Engineering Architecture Review

**Date:** 2026-06-17
**Method:** plan-eng-review applied to the existing codebase (no plan doc; reviewing live architecture). Input: the office-hours critique (`2026-06-17-pgx-suite-architecture-critique.md`) + this session's coverage-gate/exit-code work. gstack telemetry/brain plumbing skipped per project context-discipline rules.
**Verdict up front:** the foundation is sound and validated (4/4 concordant on HG002) — this is hardening, not a rewrite. But there is one regression in today's fix, one serious DRY landmine, and a clinical-provenance gap.

---

## 0. What's actually good (don't touch)

- **Per-gene output isolation is correct.** `pgx-all-genes.sh` calls `pgx-run.sh ... --output ${GENES_DIR}/${gene}`, so parallel genes don't clobber each other's `pypgx/`, `stargazer/`, or `${GENE}.vcf.gz`. No collision bug.
- **CRAM single-extraction is a smart design.** Extracting all PGx regions + MHC + chrM + VDR into one small `pgx_input.bam` once per sample (then feeding it to all per-gene calls) solves the pysam-can't-decode-CRAM problem and avoids 25× CRAM region reads. Keep it.
- **Parsers are individually defensive.** Every `parse_*` in `pgx-compare.py` wraps its body in `try/except Exception`, so one caller's malformed output can't abort the whole report.

---

## 1. Architecture

### End-to-end flow

```
run_pgx_suite.sh (host)
        │  docker run --privileged + volume mounts
        ▼
pgx-all-genes.sh ─── CRAM? ──► extract all regions → pgx_input.bam (once)
        │                                    │
        │  job pool (--jobs N)               ├──► pgx-bamstats.sh  (bg, ~6 min)
        ▼                                    ▼
   for each of ~25 genes:  pgx-run.sh <gene> <bam> --output genes/<gene>
        │
        │  Phase 1 (parallel): bcftools | Aldy | StellarPGx(nextflow) | PyPGx-prep | Stargazer-gdf
        │  Phase 2 (parallel): PyPGx pipeline | Stargazer genotype   (wait on VCF)
        │  Phase 3:            pgx-compare.py  →  <gene>_<sample>_{comparison.tsv, detail.json}
        ▼
   aggregate ───► all_genes_summary.tsv ───► pgx-report.py ───► <sample>_pgx_report.html
```

### Finding A1 — [P0, CRITICAL] The coverage gate is bypassed by the HTML report

Today's fix made `pgx-compare.py` emit `NO_CALL` to the console table, the comparison TSV, and the batch headline. But the **HTML report has its own concordance path** that reads the raw per-tool diplotypes from `detail.json` (which I deliberately preserved for audit):

```
pgx-report.py:970   most_common_dip, n_agree = Counter(diplotypes).most_common(1)[0]
```

It never reads `detail.json → coverage.verdict`. So an uncovered gene that I "fixed" still renders as a confident `*1/*1` Normal Metabolizer in the clinical HTML — the one artifact a clinician actually reads. **The fix is incomplete.**
Remedy: `pgx-report.py` must honor `coverage.verdict == "NO_CALL"` and short-circuit the gene card / detail page to a NO_CALL state before computing concordance.

### Finding A2 — [P0] Concordance logic lives in three places (and two share a tie bug)

The "most tools agree" computation is implemented independently three times:

```
1. pgx-compare.py   print_table()        console + my gate
2. pgx-all-genes.sh awk 'sort|uniq -c|sort -rn|head -1'   batch headline
3. pgx-report.py    compute_concordance()  Counter.most_common(1)   HTML
```

Two of them (`#2`, `#3`) resolve a 2-vs-2 split by arbitrary order — a tie is silently presented as consensus. Any future change to the tie rule, the coverage gate, or normalization must be made in three languages/files or they diverge. This is the structural cause of A1.
Remedy: one shared concordance+verdict function. Practically, make `pgx-compare.py` the single authority — write the final verdict (consensus / DISCORDANT / NO_CALL + n_agree/n_called) into `detail.json`, and have `pgx-report.py` and the batch summary **read** it rather than recompute.

### Finding A3 — [P2, strategic] Bash orchestration is near its complexity ceiling

The orchestrator now carries: associative config arrays in two languages, three special-case bypasses (HLA, MT-RNR1, alt-contig GSTT1), SV-gene subsets, a hand-rolled `kill -0`/`sleep 5` job pool, per-tool status files, and CRAM pre-extraction. It works, but each of these is a thing a workflow engine gives for free. This is the one place to consider spending an innovation token — see Section 5's decision.

---

## 2. Code Quality (DRY)

### Finding C1 — [P0] Gene configuration is duplicated across 3-4 files in 2 languages

Single source of truth does not exist. The same gene facts are restated in:

```
                       coords   tool-support   SV-sets   CRAM-region
pgx-run.sh             ✔ :82    ✔ :120         ✔ :153,156    ·
pgx-compare.py         (note)   ✔ :88          ✔ :75,78      ·
pgx_cram_regions.bed   ✔                ·            ·         ✔
```

Adding or editing a gene means touching three files, in bash *and* python, in two different encodings (`"1 1 1 1"` vs `{"pypgx": True, ...}`), with **no check that they agree**. The `pgx_cram_regions.bed` header literally says "from pgx-run.sh GENE_COORDS" — a comment standing in for a constraint the code can't enforce. This is the highest-probability future-bug source in the repo.
Remedy: one `genes.tsv` (or `genes.json`) — gene, region, per-tool support, SV flag, CRAM-extract region — read by `pgx-run.sh` (via a tiny parse), `pgx-compare.py`, and a generator that emits the BED. ~1 day human / ~1 hr CC. Add a unit test asserting every gene resolves consistently.

### Finding C2 — [P2] Parse errors are stuffed into the diplotype field

`parse_pypgx:271` (and siblings) do `result.diplotype = f"parse error: {exc}"`. The error rides in the value field. It's at least flagged `status=failed`, but a parse error masquerading as a diplotype string is the kind of thing that leaks into a report cell. Move it to a dedicated `error` field; keep `diplotype = "-"`.

---

## 3. Tests

### Finding T1 — [P0] The parsers — the most fragile code — have no unit tests

`pgx-compare.py` has nine parsers (~600 lines) translating nine different tool output formats. Tool output formats are exactly what drifts between versions and silently breaks a parser into `status=failed` (or worse, a wrong-but-plausible call). Today there is **one** parser-adjacent test (the coverage gate I added). There are zero fixture-based tests asserting "this real Aldy/Stargazer/StellarPGx output parses to this diplotype."
Remedy: commit one tiny real output fixture per tool under `docker/testdata/<tool>/` and a `test_parsers.py` that runs each parser against its fixture and asserts the parsed `CallerResult`. This is the highest-value test work and it's cheap. Given your "rather too many tests than too few" stance, this is squarely in scope.

### Finding T2 — [P1] No end-to-end test on a small input

Everything is validated by hand on the 167 GB HG002 BAM. There's no checked-in tiny BAM (one gene region, a few hundred reads) that exercises `pgx-run.sh` → `pgx-compare.py` → a known diplotype in CI-time. `test.sh` is smoke-only (tools exist + run). A 1-gene mini-BAM would catch orchestration regressions (like A1) before a 53-second-per-gene full run.

---

## 4. Performance

### Finding P1 — [P1] Unbounded nested parallelism, no global resource budget

```
all-genes job pool:   N genes concurrent  (--jobs N)
   each gene:         ~5 tools in parallel
      each tool:      multithreads internally (bcftools, nextflow, aldy…)
```

Peak concurrency = `N × ~5 × tool_threads`, with no global CPU/RAM cap. On a shared node or small cloud instance this oversubscribes and can thrash or OOM (StellarPGx/graphtyper and PyPGx's sklearn step are the memory-heavy ones). The job pool also polls with `sleep 5`, a crude scheduler.
Remedy (cheap): a single `--threads` budget threaded down to each tool's thread flag, and/or cap `jobs × tools`. Remedy (better): `xargs -P`/GNU parallel instead of the hand-rolled pool.

### Finding P2 — [P2] StellarPGx pays Nextflow startup 25× per sample

`run_stellarpgx` calls `nextflow run main.nf` once **per gene**. Nextflow is built to schedule a whole batch as one DAG; invoking it as a single-gene subprocess pays JVM + Nextflow cold-start (~seconds each) ~25 times and gains none of its scheduling. This is the largest avoidable per-sample overhead after PyPGx's CNV step.

### Finding P3 — [P2] No resume / idempotency

A failed run reprocesses everything from scratch; Stargazer even `rmtree`s its output dir on start. For 167 GB inputs where extraction + bamstats alone is minutes, a gene-25-of-25 failure is expensive. Resume is free in any workflow manager.

---

## 5. Prioritized findings & the one decision

| # | Finding | Severity | Effort (human / CC) |
|---|---------|----------|---------------------|
| A1 | HTML report bypasses the coverage gate (shows confident call) | **P0 regression** | 1 hr / 10 min |
| C1 | Gene config duplicated across 3-4 files, 2 languages | **P0** | 1 day / 1 hr |
| A2 | Concordance logic in 3 places; 2 share a tie bug | **P0** | ½ day / 30 min |
| T1 | Parsers (~600 lines) have no unit tests | **P0** | 1 day / 1 hr |
| (done) | Coverage gate + real exit codes | shipped today | — |
| P1 | Unbounded nested parallelism, no resource cap | P1 | ½ day / 20 min |
| T2 | No small-input end-to-end test | P1 | ½ day / 30 min |
| C2 | Parse error in diplotype field | P2 | 1 hr / 10 min |
| P2 | Nextflow invoked per-gene (25× cold start) | P2 | strategic |
| P3 | No resume/idempotency | P2 | strategic |

### The decision: harden the bash pipeline, or migrate orchestration?

Most of the P1/P2 items (C1 config duplication, missing provenance, P1 resource budget, P2 Nextflow overhead, P3 resume) are symptoms of one root choice: hand-rolled bash orchestration. There are two coherent paths.

- **Approach A — Harden bash (recommended now).** Fix A1, fold the three concordance implementations into one authority in `detail.json`, single-source the gene config, add parser fixtures + a mini-BAM e2e test, add a `--threads` budget, and write provenance (tool versions, reference build, pipeline git SHA, timestamp, command lines) into `detail.json`. Keeps boring tech, right-sized, no rewrite. Clears every P0 and most P1.
- **Approach B — Migrate to a workflow manager (Snakemake) (strategic, defer).** Python-native (fits the codebase), gives DAG scheduling, resume, per-rule resource limits, and provenance for free — consolidating C1/P1/P2/P3 structurally. But it's an innovation token and a real migration; only worth it when the goal shifts from occasional single-sample runs to batch throughput across many samples.

Recommendation: **A now** (A1 is an active clinical-correctness regression), **B only when batch scale becomes the actual requirement.**
