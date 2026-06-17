# pgx-suite — Adversarial Architecture Review

**Date:** 2026-06-17
**Mode:** Engineering/infra critique via the office-hours posture (anti-sycophancy, premise challenge, forced alternatives). The YC demand questions were skipped — this is a pipeline architecture review, not a product-market diagnostic.
**Rule honored:** no code written or modified. This is a strategy brief only.

---

## 0. Premise correction (read before anything else)

Your Block 1 asks me to help you "trim the voter graph down to a single dominant engine" because you are "spinning up multiple heavy, parallel sequence-alignment algorithms and hyper-variable localized assembly engines, focusing on heavy HLA-calling sub-modules, targeting overlapping polymorphic regions."

Three of those claims are false against the actual code, and a sycophantic answer would have you cutting the wrong thing.

1. **None of your star-allele callers do sequence alignment.** PyPGx, Stargazer, Aldy and StellarPGx all consume an *already-aligned* BAM/CRAM. They are genotypers/diplotype callers, not aligners. The only shared compute is one `bcftools mpileup` per gene. Measured cost: ~53s wall for CYP2D6 on a 167 GB WGS (your own timing profile). That is not "heavy."

2. **You do not have multiple HLA engines.** HLA-A/B are called by exactly **one** engine — OptiType — and it *bypasses* the entire star-allele pipeline (`GENE_SUPPORT[HLA-*] = "0 0 0 0"`, handled by the `run_hla()` special case in `pgx-run.sh`). There is no overlap, no voter graph, no reconciliation on HLA. There is nothing to trim there.

3. **The callers do not target overlapping regions.** Each runs per-gene on disjoint coordinates from `GENE_COORDS`. The only local-assembly components are StellarPGx (graphtyper) and OptiType (ILP over reads) — and they operate on different loci entirely.

So the "many heavy concurrent consensus callers competing on HLA" you want to kill does not exist. Your real costs live somewhere else (Sections 1, 3, 4). Redirect the energy.

---

## 1. Block 1 — the multi-caller "voter graph": keep it, it's your product

The 4-caller concordance is not redundant compute you tolerate. For a clinical (CAP/CLIA) report it **is** the product — independent callers agreeing is your defensibility story. Cutting to a single engine doesn't streamline a clinical product; it removes the reason a lab would trust the output.

And it's cheap. The callers are I/O-light genotypers running in parallel; the gating cost is PyPGx's sklearn CNV step (~49s of the ~53s), not "too many engines."

**Where the real cost actually is:**

- **Operational weight, not CPU.** StellarPGx drags in Nextflow + a JVM + Apptainer + a SIF + a JAR that downloads on first run. That is the heaviest, most fragile dependency in the stack by a wide margin — and it's one caller, not the consensus model.
- **Image size.** 4 Python/Java toolchains + Nextflow + Apptainer + the StellarPGx SIF + the ~500 MB OptiType SIF + the ~500 MB Beagle panel. The bloat is bundling, not concurrency.
- **Reconciliation overhead** is real but already paid: `pgx-compare.py` carries a hand-written parser per tool because each writes a different format in a different place. That's a maintenance tax, not a runtime tax.

**Brutal version:** if you want a single "dominant, resilient engine," the honest pick is **Aldy** — pip-installable, no JVM, no Nextflow, no Apptainer, no Beagle panel, and its ILP solver handles CNV/SV natively. It covers the widest single-tool gene set in your matrix. The genes it misses (POR, GSTT1, CACNA1S) fall back to PyPGx. **But** picking it means deleting your concordance signal. That is a clinical-vs-research decision (Section 5), not a performance optimization.

## 2. Block 2 — the 1-star baseline

Stripped to the minimum that still emits diplotypes:

- **Keep:** one caller (Aldy is leanest; PyPGx if you need its wider gene table + CACNA1S/POR/GSTT1), `bcftools` for the few SNP-only genes, a thin TSV writer.
- **Delete:** StellarPGx + Nextflow + Apptainer + `--privileged`, OptiType (HLA), mutserve (MT-RNR1), Stargazer, Aldy *or* PyPGx (whichever you didn't pick), `pgx-compare.py` consensus logic, the HTML report, `pgx-bamstats.sh`.
- **Result:** a single unprivileged container, ~10× smaller image, deploys anywhere, runs in seconds.

What you lose is exactly what makes it a clinical product: multi-caller concordance, HLA risk alleles (abacavir/carbamazepine), mtDNA aminoglycoside risk, and the CAP/CLIA report. The 1-star is a **research-grade genotyper**, not a clinical suite. That's a fine product — it's just a *different* product.

## 3. Block 3 — persistence & inter-module blind spots (this is the actual fire)

This is your strongest question and where I'd spend the next sprint. The silent failure modes are real and concentrated:

1. **Every tool runner returns success unconditionally.** `run_aldy`, `run_stellarpgx`, `run_pypgx_pipeline`, `run_stargazer_genotype`, `run_hla`, `run_mt` all `return 0` even on failure (Aldy is explicitly commented "always return 0"). A crashed caller is indistinguishable downstream from a caller that ran and found nothing. The script is `set -uo pipefail` but **not** `-e`, so nothing propagates.

2. **Concordance ties resolve silently and arbitrarily.** `pgx-all-genes.sh` picks the "top" diplotype with `... | sort | uniq -c | sort -rn | head -1`. A 2-vs-2 split between callers doesn't flag as discordant — `sort` order picks a "winner." For a clinical report, a coin-flip presented as consensus is the most dangerous bug in the repo.

3. **No coverage floor before calling.** I see no per-gene gate that hard-fails when a region has no/low reads. `pgx-bamstats.sh` runs *alongside* calling for the report, not *before* it as a gate. A gene with zero coverage produces an empty/garbage diplotype rather than an explicit "insufficient coverage — no call." That is the classic clinical silent-failure.

4. **VCF failures degrade silently.** `bcftools` failure logs "WARN VCF may be missing" and continues; VCF-dependent callers then fail and (per #1) report success with no output.

5. **Known QC corruption already in your tree.** Your own TODO notes `bam_stats.json` has `read_length=0` for CRAM samples (mosdepth-on-CRAM path). QC silently wrong is worse than QC absent.

**Fix direction (no code here, just the shape):** give each runner a real exit code; treat "tool failed" and "tool reported no result" as distinct states in `pgx-compare.py`; emit an explicit `DISCORDANT` status when callers disagree instead of `head -1`; add a coverage precondition that produces a loud `NO_CALL / INSUFFICIENT_COVERAGE` row. None of this requires removing a single caller — it makes the callers you have trustworthy.

## 4. Block 4 — deployment & hardware boundaries

- **`--privileged` is your enterprise killer, and it's 100% StellarPGx.** Apptainer-in-Docker needs it. GKE Autopilot, ECS/Fargate, and most shared HPC and managed K8s **forbid** privileged containers. As long as StellarPGx-via-Apptainer is in the default path, "enterprise-ready cloud deployment" is blocked — not by performance, by the security model. This single dependency, not the caller count, is what gates productionization.
- **Memory spikes** come from (a) PyPGx's sklearn CNV step and (b) StellarPGx/graphtyper assembly. In all-genes mode, `--jobs N` runs N genes × ~5 tools concurrently with no global core/RAM budget — tools each multithread independently, so you can oversubscribe badly on a shared node.
- **I/O, not CPU, is the wall.** `pgx-bamstats.sh` is ~6 min on 167 GB WGS at 490 MB/s; mosdepth-on-CRAM is the random-seek bottleneck you already split out.
- **Storage retention:** Nextflow `.work/` dirs, per-gene outputs, baked/mounted SIFs and the Beagle panel. Nextflow work dirs in particular accumulate silently.

## 5. The governing fork (the decision only you can make)

Everything above collapses to one question, and your codebase contradicts your brief on it:

- The brief says: cut to a minimal single engine for deployment velocity.
- The code says: CAP/CLIA fields, a clinical HTML report, multi-caller concordance — i.e. a **clinical** product where concordance is the value.

You can't have both the 1-star single-engine baseline *and* the clinical-defensibility story. Pick the product:

### Alternatives

**APPROACH A — Harden the clinical suite (recommended if the goal is clinical).**
Keep all callers. Spend the sprint on Section 3: real exit codes, explicit DISCORDANT/NO_CALL states, coverage gate. Separately attack the *real* deployment blocker by making StellarPGx optional/removable so the default path drops `--privileged`.
Effort: M. Risk: Low. Reuses: everything; this is mostly correctness wiring.
Net: keeps the product that justifies the complexity, fixes the failures that actually threaten it.

**APPROACH B — Research-grade 1-star (recommended if the goal is velocity/scale).**
Single engine (Aldy), bcftools for SNP-only genes, TSV out. Delete StellarPGx/Nextflow/Apptainer/`--privileged`, HLA, mtDNA, consensus, HTML.
Effort: S–M (mostly deletion). Risk: Med (loses clinical claims). Reuses: one caller + the gene-coords table.
Net: deploys anywhere, ~10× smaller, seconds to run — but it's a different, non-clinical product.

**APPROACH C — Tiered build (the lateral option).**
One repo, two image targets: a lean `pgx-core` (Approach B, unprivileged, the default enterprise artifact) and a `pgx-clinical` superset (Approach A, privileged, for the lab pipeline). The split is a build arg + making StellarPGx/OptiType/mutserve opt-in layers.
Effort: M. Risk: Low–Med. Reuses: the whole tree, reorganized.
Net: you stop forcing one image to be both the fast research tool and the heavy clinical instrument — which is the actual root tension in your brief.

## The Assignment

Before any cutting: **run the suite on one sample with a deliberately uncovered gene** (e.g. point a gene's coordinates at an empty region, or use a targeted BAM missing a locus) and watch what the final report says. I predict it shows a confident diplotype rather than "no call." If I'm right, your problem was never too many engines — it's that the engines you have fail silently. Fix that before you delete anything.
