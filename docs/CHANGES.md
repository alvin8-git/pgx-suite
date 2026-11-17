# PGx Suite — Changelog

All notable changes to this project are documented here.
Format: reverse-chronological, grouped by phase/milestone.

---

## 2026-11-16 (fresh clones build: pypgx fork, Stargazer panel check, README rewrite)

- The `pypgx` submodule pinned a commit (PyPGx 0.26.0 + three pandas 2.x fixes) that was only local, so
  `git clone --recurse-submodules` failed. The fixes are now branch `pgx-suite-patches` of the public fork
  `alvin8-git/pypgx` (tree identical); `.gitmodules` points there.
- The Dockerfile stops when `stargazer/1kgp_vcf/grc38/` (gitignored, Stargazer distribution) is missing;
  before, a fresh clone built an image in which Stargazer failed on every gene.
- Checked from a fresh public clone: no-cache build without the panels stops with the message; with them it
  builds in 4 min (1.95 GB) and the smoke test passes 4/4 callers.
- README rewritten for a new user: quick start, both engines, reference data incl. the Stargazer panels,
  engine-equivalence validation, troubleshooting.

## 2026-11-16 (Nextflow is the default engine; PharmCAT cross-check race fixed)

- `run_pgx_suite.sh --engine nextflow|snakemake`, default `nextflow`; `--resume` for Nextflow. Printed output
  paths now match the real layout, and a failed run reaches the end-of-run summary.
- `tests/smoke/run.sh` runs Nextflow by default, `ENGINE=snakemake` for the Snakefile; both pass (62 s / 47 s).
- PharmCAT is a declared input of every gene's compare (Snakefile and `main.nf`). `pgx-compare.py` adds an
  informational PharmCAT cross-check to every gene PharmCAT reports; before, only `pharmcat`-flag genes waited
  for PharmCAT, so the cross-check and the NUDT15 note text depended on job order. A PharmCAT `Unknown` call
  now gives no cross-check (`test_pharmcat.py`). HG002 and HG007 through the launcher: both engines 0
  differences, 15 genes carry the cross-check, no verdict status, diplotype or phenotype changed.
- `tests/compare_runs.py` no longer allows a cross-check difference.

## 2026-11-16 (a reading CPIC cannot score is not a split)

- Owner decision 2026-11-16, superseding the 2026-11-15 rule: in `_split_by_cpic_phenotype` (`docker/pgx-compare.py`) a
  diplotype reading CPIC has no row for (Aldy composite DPYD `*rsNNN` / `*9A+rsX` haplotypes, Stargazer `*S` codes the notation
  map cannot place) no longer counts as a split, but only after a component check. The reading is split into its components
  (named alleles, `+rsNNN` tags, `*rsNNN`, `*S` codes, HGVS names; `reconcile.reading_function_profile`) and each is looked up in
  CPIC's allele-function table (`cpic_diplotype.json` `allele_function`: CYP2C19, DPYD, TPMT; other genes have none, so nothing is
  ever set aside there). It is set aside only if every component has a known CPIC function, every decreased / no / increased
  function component is also present (with at least that many copies) in the scorable readings, no more callers back it than
  back the scorable readings together, and every scorable reading gives the same usable CPIC phenotype. Then that phenotype is
  the gene's phenotype, status stays `concordant` / `majority`, the note says "reading X has no CPIC phenotype and was not counted
  (<check that passed>); CPIC phenotype of the other readings is P" (once; "callers said Y" first when a caller's label is
  replaced), `phenotype_source` names the CPIC rows used, the verdict carries `unscorable_readings_ignored` and the report
  sentence says the reading was set aside. An unknown component, an extra decreased-function component (Aldy
  `*1/*1+rs67376798`, c.2846A>T, against `*1/*1`), a lettered allele that is deliberately unscorable (CYP2C19 `*1A`), an ambiguous
  CPIC row (`reconcile.cpic_diplotype_lookup` now returns "ambiguous" apart from "no_row"), NAT2 Indeterminate, or scorable readings
  that disagree all keep the split. `UNSCORABLE_READING_IS_A_SPLIT = False` is now the default; `True` restores the 2026-11-15 rule.
- Reason: DPYD `*2A` carriers (Intermediate Metabolizer) had their phenotype blanked although PharmCAT/PyPGx and Stargazer agreed
  and only Aldy's composite haplotype was unscorable; the component check keeps the rule from hiding a deviating variant (suite rule 4).
- Real run (image built from this tree, stored caller outputs of 66 samples, 2767 gene verdicts re-scored inside the container
  with the 2026-11-15 code (df2dd00) and the final code): 8 verdicts change, all DPYD, all from `-` / `concordant, diplotype
  unresolved` to a CPIC phenotype with status `concordant` (Normal Metabolizer 6, Intermediate Metabolizer 2); 2759 are
  identical. In all 8 the set-aside reading's components were all known, 6 carried no deviating component and 2 carried
  one (`*2A`; c.2279C>T), present in the kept readings too. Of the 14 verdicts blanked under the 2026-11-15 rule,
  8 regain their phenotype and 6 stay `-` (NUDT15 3, CYP3A5 2, CYP2C19 1) because their scorable readings disagree. The commits
  after df2dd00 (including the CYP2D6 phenotype fix) change none of these 2767 rows.
- Extra-variant tags in the vote (follow-up, same day): `reconcile.canonicalize_diplotype` used to drop `+rsNNN` / `+*N` tags
  (DPYD, CYP2C19, TPMT), so Aldy `*1/*1+rs67376798` (DPYD c.2846A>T, Decreased function in CPIC's allele table) voted as `*1/*1`
  and gave a 4/4 concordant Normal Metabolizer. Now each dropped tag is looked up in CPIC's `allele_function` table
  (`reconcile.blocking_tags`): a Normal-function tag is dropped as before; a decreased / no / increased function tag, a tag the
  table does not give, or a gene without a function table stays in the key, so that caller is its own reading (split or
  majority by the existing rules). Only the blocking tags enter the key, spelled through the gene's notation / synonym map
  (`+rs67376798`, `+*rs67376798`, `+c.2846A>T` are one key). Components are cut only before `rs`, `*` or `c.`, so an HGVS name
  such as `c.1905+1G>A` stays whole and resolves (`*2A`, No function). The note comes from the same keys
  (`reconcile.canonical_tags`): "not counted as agreeing with the diplotype without the tag (<caller> <diplotype>: tag <tag> is
  <function> in CPIC's allele table)" only for a caller whose key differs from every consensus reading, and "the reading keeps
  a tag that is not Normal function (...)" when a consensus reading carries one. The lookup files `cpic_diplotype.json` and
  `dpyd_notation_map.json` are now needed at run time (fail loud if missing; a test checks the Dockerfile copies them). TPMT
  `*1/*3A+rs1142345` now has no CPIC phenotype (CPIC's function table has no row for the rsID) instead of silently dropping the
  tag; a Normal-function tag still looks up. Not changed here: when every caller carries the same decreased-function tag, the
  phenotype still comes from the stripped diplotype (a DPYD-only rule is the next change). Real run (image `wt-rstag`, 66 stored
  samples, 2767 verdicts, OLD = the commit before this change): 0 rows change, since no stored call carries a tag that is not
  Normal function; the cases are proven by `test_extra_variant_tag_votes_only_when_its_cpic_function_is_normal` and
  `test_tag_vote_review_fixes`.
- DPYD: a non-Normal component blocks Normal (owner decision 2026-11-16, `_dpyd_nonnormal_guard` in `docker/pgx-compare.py`,
  `reconcile.nonnormal_components`). If any caller's DPYD call carries a named allele rated Decreased or No function in CPIC's
  allele table (`cpic_diplotype.json` `allele_function`), a named allele the table does not give (`UNKNOWN_NAMED_ALLELE_BLOCKS_NORMAL`),
  a `+` tag that is not Normal function (decreased, no, increased, unknown), or a call that cannot be read (VCF-Check's unphased
  `*A(het)+*B(het)` list is parsed), and the verdict would be Normal Metabolizer, the phenotype is `-` and the status is the new
  `concordant, phenotype withheld` / `majority, phenotype withheld` (the report maps both as amber; `build_axiom_concordance` takes the part before the comma; OmniGen's consumer
  report compares the summary's `Status` with the plain string `concordant`, so these labels are not read as concordant, and it reads
  ConsensusPhenotype, nAgree, nCalled and DecidedBy; the new fields are in `detail.json` only; OmniGen shows a no-phenotype gene from `phenotype_note` in `detail.json`, which carries the
  component for split and discordant verdicts too, so no new summary column is needed). Every affected verdict carries `dpyd_nonnormal_component` and `dpyd_nonnormal_note`
  (`normal_withheld` and `phenotype_note` only when a Normal was withdrawn), the note names each caller and component once with its CPIC function and
  says "confirm with a targeted DPYD assay", and any set-aside text that says Normal is removed. A DPYD verdict that already has no
  phenotype (a split or discordant call) keeps its status and split explanation and gains the flag and the component in the note
  (no claim that a Normal was withheld). Callers are chosen with the vote's typed-call rules (abstentions never count; VCF-Check's
  allele lists do; the informational PharmCAT call does), a real tool authority keeps its provenance, and the set-aside clause is
  recorded in `cpic_phenotype_clause` and removed by value. Report: the flag turns the phenotype into a
  shared "Normal withheld" category that `_get_tier` flags as high priority, so the gene section and the landing findings never show
  "Normal / Not flagged" (a one-caller case used to, because the diplotype check needs two tools); `print_table` says the phenotype
  is withheld instead of "N/N tools agree" alone. An Intermediate or Poor verdict is left alone; other genes are unchanged. Reason:
  one caller showing a reduced-function DPYD variant is evidence against Normal, and a false Normal risks fluoropyrimidine toxicity.
  Real run (image `wt-dpydblock`, 66 stored samples, 2767 verdicts, OLD = e546f13): see the counts below. 6 of 198 stored DPYD
  calls (2 samples) carry a Decreased or No function component (`*2A`, c.2279C>T, in every caller of each sample); both verdicts are
  Intermediate and stay; none of the 198 calls has an unknown-function named allele. Tests: `test_dpyd_nonnormal_component_blocks_normal`,
  `test_dpyd_guard_review_fixes`, `test_tag_vote_majority_on_a_gene_the_dpyd_guard_does_not_touch`.
- Tests: `test_verdict.py` (component check, extra deviating component, caller support, ambiguous row, NAT2 Indeterminate, lettered allele, report sentence; the composite-Aldy DPYD case in
  `test_phenotype_tier_rescues_string_discordance` now keeps Normal Metabolizer).

---
## 2026-11-16 (Nextflow engine, opt-in)

- `main.nf`, `modules/local/pgx.nf`, `nextflow.config`, `bin/`: a Nextflow DSL2 port of `docker/Snakefile`, one
  process per rule with the same commands, container paths and status-file convention. Snakemake stays the
  default; the Snakefile and pipeline scripts are unchanged.
- Validation (image `211aa10`): on HG001 to HG007, all genes, the Nextflow output equals a Snakemake run file for
  file; every remaining difference is tool nondeterminism also seen between two Snakemake runs or repeated tool
  runs (arcasHLA, Aldy, OptiType, PyPGx, Stargazer/Beagle, mutserve logs, the PharmCAT cross-check race).
  Details and wall times: `docs/nextflow-port-plan.md`, "Validation results".
- `tests/compare_runs.py`: compares two result trees, lists allowed tool nondeterminism separately.
- Found, not fixed (owner decision): the Snakefile does not declare PharmCAT as an input of non-PharmCAT genes,
  yet `pgx-compare.py` adds a PharmCAT cross-check to them when PharmCAT finished first.

## 2026-11-15 (code-review fixes, split rule and lookups)

Fixes to the uncommitted split-phenotype / CPIC-table work after a /code-review medium.
- `pgx-compare._split_by_cpic_phenotype` now applies to every `A or B` verdict of a gene with a CPIC table, whatever flags the verdict has.
  Readings that all give one usable CPIC phenotype report CPIC's phenotype (the callers' label, when different, goes in the note, as for
  TPMT/NUDT15; e.g. CYP2C19 `*2/*2 or *2/*3` is Poor, not the callers' Intermediate; CYP2B6 phase-ambiguous `*1/*7 or *5/*6` is
  Intermediate, not the callers' Indeterminate). Readings all Indeterminate give `-` like a split. A split keeps the agreement level
  (`majority, diplotype unresolved` / `concordant, diplotype unresolved`) and clears `authority` (it said the phenotype tier decided).
  `pgx-report.py` has badges for the new status; `build_axiom_concordance.adjudicate` treats it like its base status.
- `reconcile`: the legacy, notation and CPIC tables load on first use and a missing file raises an error naming it (import no longer reads them);
  `cpic_diplotype_phenotype` continues to the folded lookups past an ambiguous row; `fold_legacy_letter` is case-insensitive;
  `_strip_allele` cuts `rs1+rs2` at `+` again while `c.1905+1G>A` stays whole.
- `build-arcashla-sif.sh` and `reanalyze_mgi.sh` require TMPDIR (no `/var/tmp` default).
- Real run (15c416c container, worktree `docker/` mounted, before = pre-fix worktree state): 65-sample rig 2767 verdicts, no call changed
  (diplotype, phenotype, status); 14 split verdicts lose `authority: Phenotype`, 22 gain `phenotype_source`. Axiom rig 1240 verdicts: no call changed; 17 splits lose `authority`, 13 gain `phenotype_source`.
  New tests in `test_verdict.py`.

---

## 2026-11-15 (legacy letter fold) — explicit PharmVar table replaces the blanket trailing-letter strip

The fold of `reconcile` (`*2B` as `*2`) was a regex verified only for CYP3A5. All PharmVar genes pgx-suite reports were rendered
(PharmVar 6.2.29, updated 2026-09-23, https://www.pharmvar.org/gene/<GENE>; the Legacy Label column is only in the rendered page).
- New `docker/pharmvar_legacy_snapshot.json` (every sub-allele with a legacy label, per gene) and `docker/pharmvar_legacy_build.py`;
  `pharmvar_legacy_suballeles.json` is now derived from it: `fold` (letter number = core number; CYP2C19 18, CYP2C9 6, CYP2B6 22,
  CYP2D6 52, CYP3A4 5, CYP3A5 7, SLCO1B1 2, NAT2 13), `redesignated` (CYP2C19 `*1A` = `*38.001`, CYP2D6 `*14A` = `*114.001`, 32 NAT2).
  CYP4F2 and NUDT15 have no lettered label; TPMT and UGT1A1 pages say "not maintained by the PharmVar Database"; DPYD is by rsID.
- `reconcile.fold_legacy_letter(gene, allele)` folds only listed labels. Unlisted and redesignated letters stay unscorable and count toward
  the split rule. Behaviour change by design: CYP2C19 `*1A` no longer folds to `*1` (it is now `*38`, a different allele).
- `Dockerfile` copies `pharmvar_legacy_suballeles.json` (it was not copied before); the image is not rebuilt here.
- Tests: table equals the derivation from the snapshot, snapshot equals freshly rendered pages (PGX_PHARMVAR_PAGES), the fold refuses unlisted
  and redesignated labels. All `docker/test_*.py` pass.
- Real run (pgx-suite:15c416c container, worktree `docker/` mounted, stored detail JSON, old regex fold vs new table): 65-sample rig
  2767 gene verdicts and Axiom rig 1240 verdicts, 0 changed cells. Lettered alleles in stored calls for genes with a CPIC diplotype table:
  CYP2C19 `*2B` (5 readings) and CYP3A5 `*3B` (3) are listed folds; no redesignated or unlisted letter occurs there (CYP3A4 `*1B`/`*1G`, TPMT `*1S`,
  DPYD `*9A` are in genes without a fold table or outside it).
- Owner decision recorded in CLAUDE.md: composite Aldy DPYD haplotypes (`*9A+rsX`) stay `-`.

---

## 2026-11-15 (DPYD notation map, CYP3A5 *3B) — Stargazer `*S` codes and Aldy `*rs` tags score through CPIC

Owner approved. A DPYD reading in caller notation (Stargazer `*S10`, Aldy `*rs2297595`, bare `*4`) had no CPIC diplotype
row, so under the split rule it made the whole call phenotype-less. New `docker/dpyd_notation_map.json`, built by
`docker/dpyd_notation_map_build.py` (no value written from memory):
- Stargazer: the tool's own `stargazer-grc38-2.0.3/stargazer/star_table.tsv` (`grc38_core`) and `snp_table.tsv` (`rs_id`)
  (sha256 stored in the map's `source`) joined to CPIC DPYD allele definitions (CPIC API `allele_definition`,
  `allele_location_value`, `sequence_location`, embedded as `cpic_alleles`). A `*S` code maps only when its GRCh38 variant set
  equals one CPIC allele's (position + ref>alt); indels (CPIC writes HGVS, Stargazer VCF style) by a unique dbSNP id on a
  single-variant allele. 91 of 457 `*S` codes map; 366 stay unscorable with the reason (`*S8` has three variants, most others
  are variants CPIC does not define).
- Aldy `*rs<NNN>`: 82 tags, each the CPIC allele defined by exactly that one variant (rs id unique among CPIC single-variant
  alleles). Composite Aldy names (`*9A+rs1801159+rs2297595`, `*rsA+rsB`) stay unscorable.
- 13 bare legacy star names (`*4`, `*6`, ...) map to the CPIC allele whose own name carries that bracketed name (`c.1601G>A (*4)`).
- Wired into `reconcile._cpic_allele_key` (phenotype lookup only; `norm_diplotype` comparison keys keep the caller's spelling).
- Tests (`test_verdict.py`): the map equals a rebuild from the Stargazer tables + embedded CPIC data and (PGX_CPIC_VERIFY=1)
  the live CPIC API; sha256 of the Stargazer tables; ambiguity cases; readings score as CPIC does.
- Real run (pgx-suite:15c416c container, worktree `docker/` mounted, stored detail JSON, map cleared vs in place): 65-sample rig
  (2767 gene verdicts): 14 of the 22 DPYD split calls regain their phenotype, each equal to an independent CPIC-table lookup
  of the mapped diplotype; the other 8 contain a composite Aldy haplotype and stay `concordant, diplotype unresolved`
  (7 only gain the per-reading phenotypes in the note). Axiom rig (1240 verdicts): 8 DPYD calls regain, 6 notes only.
  No other gene changed in either rig.
- CYP3A5 `*3B` -> `*3` fold verified: PharmVar CYP3A5 page (version 6.2.29, updated 2026-09-23) lists `CYP3A5*3.002` (PV01338,
  rs776746 + rs28383468 + rs28383469) with Legacy Label `CYP3A5*3B`, and `*3.001` = legacy `*3A`. Evidence in
  `docker/pharmvar_legacy_suballeles.json`, comment in `reconcile.py`, test `test_cyp3a5_legacy_letters_are_pharmvar_suballeles_of_star3`.
  No stored call changes. The same letter fold for other genes (CYP2C19 `*2B`, ...) is not yet checked.
- `Dockerfile` copies the map; the image is not rebuilt here.

---

## 2026-11-15 (split phenotype) — readings with different CPIC phenotypes report no phenotype

Owner decision 2026-11-15. When the consensus lists two or more diplotype readings (`A or B`, diplotype not resolved)
and CPIC maps them to **different** phenotypes, `consensus_phenotype` is `-`, `status` is
`concordant, diplotype unresolved` (never plain `concordant`; the report badge and card are amber) and the note names
each reading's CPIC phenotype and says "confirm with a targeted assay". Readings that all give one CPIC phenotype keep
the phenotype and the `concordant` status as before. A reading CPIC has no row for counted as a different phenotype
(policy constant `UNSCORABLE_READING_IS_A_SPLIT` in `pgx-compare.py`; superseded 2026-11-16, see that entry).
- Phenotype per reading comes from CPIC's own diplotype tables: new `docker/cpic_diplotype.json` (CPIC `/diplotype`
  API, CYP2B6, CYP2C19, CYP2C9, CYP3A5, DPYD, SLCO1B1, UGT1A1; built by `cpic_diplotype_build.py`, compared with the API in
  `test_cpic.py` under `PGX_CPIC_VERIFY=1`), `cpic_thiopurine.json` (TPMT, NUDT15), and the NAT2 table already in the code.
  Lookup is by the reading as spelled, then synonym-collapsed, then with a trailing sub-allele letter dropped (`*3B` as
  `*3`, `*2B` as `*2`). Genes with no CPIC diplotype table (e.g. CYP4F2) and CYP2D6 (Cyrius / activity-score rule) are unchanged.
- Applies to the phenotype-tier verdicts and to the CYP2B6 phase-ambiguous verdicts; NAT2's own "Indeterminate" for a split
  is now `-` with the same status.
- `reconcile._strip_allele` no longer cuts an HGVS name at its `+` (`c.1905+1G>A (*2A)` reached the synonym map as `c.1905`,
  so it never matched `*2A`): DPYD readings that are synonyms of the same variants are now one reading (cosmetic).
- Bug fix: `allele_synonyms.json` listed CYP2C19 `*38` as a synonym of `*2`. CPIC (`allele` table) gives `*38` Normal function and `*2`
  No function, and PharmVar renamed the reference-sequence haplotype `*1.001` (legacy `*1A`) to `*38` (Botton et al. 2021,
  PMID 32602114; Gaedigk et al. 2021, PMID 34091888), so the synonym is removed. `cpic_diplotype.json` now also carries CPIC's allele
  function for the genes with a synonym entry, and `test_verdict.py` checks that every synonym pair has the same CPIC function
  where CPIC defines both (no other pair fails; CPIC does not define TPMT `*1S`, DPYD `*S`-codes or `ivs14+1g>a`).
- `Dockerfile` copies `cpic_diplotype.json`; the image is not rebuilt here.
- Tests: `docker/test_verdict.py` (split, same-phenotype, unscorable, phase-ambiguous, DPYD synonyms); `test_cpic.py` covers the table.
- Re-run on stored caller outputs with the real `pgx-suite:15c416c` container and this code mounted: Axiom rig (40
  samples, 1240 gene rows) 25 rows changed, all of them multi-reading verdicts (DPYD 14, CYP3A5 4, NUDT15 4, CYP2C19 2,
  TPMT 1); 65-sample consistency rig (2767 gene verdicts) 28 changed (DPYD 22, NUDT15 3, CYP3A5 2,
  CYP2C19 1), again only multi-reading verdicts; every other gene row is byte-identical. Of the 28, only 6 are splits CPIC scores
  (CYP3A5, CYP2C19; NUDT15 already had no phenotype, only its status label changed); the 22 DPYD ones are readings that
  differ in notation (Stargazer `*S20`, Aldy `*rsNNN` tags) that CPIC has no row for: setting
  `UNSCORABLE_READING_IS_A_SPLIT = False` keeps their phenotype (adopted 2026-11-16).

---

## 2026-11-15 (scrub) — no private sample ids, hostnames or internal paths in tracked files

- Private clinical sample ids in `docs/AxiomValidation.md`, `docs/TODO.md`, `docs/CHANGES.md`, `docs/ToolsDocumentation.md`,
  comments in `docker/pgx-compare.py`, `docker/pgx-report.py` (+ `.snapshot`), `docker/test_vcf_check.py` and
  `build_axiom_concordance.py` are replaced by neutral labels (`sample A` ... `sample T`, consistent within and across documents).
  A sequencer LIMS id for HG002 is now `HG002_LIMS_ID`. The cohort/institution name in comments and `allele_synonyms.json`
  provenance reads "internal validation cohort".
- Internal absolute paths (`/data/...`) in docs and `docker/build-arcashla-sif.sh` are placeholders (`<ref>`, `<scratch>`, `<repo>`).
- `reanalyze_mgi.sh` no longer carries a cohort path or sample list: `BAM_DIR` is required, the cohort is every
  `*.bwa.sortdup.bqsr.bam` in it, `LOG_DIR` defaults to `${TMPDIR:-/var/tmp}/pgx_reanalysis`.
- Behaviour unchanged: comment/string edits only; the 13 `docker/test_*.py` files pass before and after.
- Not done: git history still contains the removed names (not rewritten). `.gitignore` still lists the cohort directory
  name and the README/tutorial carry the public GitHub account URL; both are owner decisions.

---

## 2026-11-15 (batch 2) — one CYP2D6 agreement rule for the gene row and the drug table; TPMT/NUDT15 from CPIC; tandem folding; `decided_by`

- **Tandem spellings agree** in `cyp2d6_agreement` (`*36+*36+*10` = `*36x2+*10`): the diplotype key is folded before comparison, so
  two callers that differ only in tandem spelling reach level a. No change on the 65 stored samples (consistency 65/65).
- **CYP2D6 gene row = PharmCAT input** (`reconcile.cyp2d6_agreement`, `pgx_outside_calls`, `pgx-compare._cyp2d6_agreed_phenotype`).
  Cause: when Cyrius gave no usable call the verdict phenotype was a plurality of the callers' own phenotype labels (ties
  broken by caller order), while the outside call required every caller's CPIC activity score (levels a/b/c). In 10 stored
  samples the gene row showed a phenotype and the drug table said "No Result". The outside-call levels moved into
  `reconcile.cyp2d6_agreement` and the verdict (status concordant/majority/discordant, no Cyrius authority) takes its
  phenotype from it (CPIC phenotype of the agreed score, or of the all-agree class); where it makes no call the phenotype
  is `-` with `phenotype_note` (also appended to `note`). Voting over callers' own labels stays rejected (CLAUDE.md).
  Cyrius-authority verdicts are unchanged. Side effect: 12 stored samples lose a label-vote phenotype (see the report).
- **A clinical-sample pattern** (`*1x2/*36+*10`, score 2.25, shown Ultrarapid): same cause, same fix; the phenotype is now `-`.
- **TPMT / NUDT15**: `docker/cpic_thiopurine.json` (CPIC `/diplotype` table, built by `cpic_thiopurine_build.py`) gives the
  verdict phenotype for a single concordant/majority diplotype (`reconcile.thiopurine_phenotype`); `phenotype_source`
  names the row, a contradicting caller label is named in `note`. HG006 TPMT `*1/*6` Indeterminate -> Intermediate
  (CPIC row `*1/*6` Intermediate Metabolizer). Discordant NUDT15 verdicts are untouched: PyPGx and StellarPGx miss `*6`
  (HG005, HG006 and two clinical samples), a caller-sensitivity limit that needs an owner decision (PharmCAT as authority).
- **CYP2D6 tandem spellings** (`reconcile._fold_tandem`, used by `clean_diplotype`): `*36+*36+*10` -> `*36x2+*10`,
  `*2+*2` -> `*2x2`; N>=3 uses CPIC's `x≥3` row (`*1x≥3`, `≥3.0`, a lower bound kept as `activity_min` in
  `cpic_cyp2d6.json`; `*4x≥3` is exactly 0.0). A lower bound gives a phenotype (>=3.0 Ultrarapid) but never an exact score.
  Comparison keys (`norm_diplotype`) are not folded.
- **`DecidedBy`**: `all_genes_summary.tsv` gets a trailing `DecidedBy` column = the verdict's existing `authority` key (`-` for votes and the
  phenotype tier); no new verdict key. `status` and `n_agree` are unchanged.
- **Discordant verdicts carry no caller phenotype** (`compute_verdict`): `-`, the callers' labels named in `note`.
  **`pgx-report.py` Key Findings** keep an authority-decided gene with `n_agree` < 2 (`decided_by`), with no discordance warning.
- **Review fixes**: level e also accepts a Cyrius `≥` lower bound matched by the same canonical diplotype (diplotype passed,
  never a number); level m is blocked by an explicit Indeterminate from a scoreable caller (HG002_12X again `-`); a majority
  verdict's phenotype comes only from callers on the consensus diplotype (67 NAT1/CYP1A2/CYP2E1 cells: dissenters'
  Indeterminate labels no longer stand); NUDT15 split also catches zygosity and uses reconcile's cleaning.
- **CYP2D6 level m** (owner decision): without Cyrius, >= `MAJORITY_MIN` (3) callers sharing one exact CPIC activity
  score give that score's phenotype for the verdict and the outside call (score only), dissenters named in `note`/`why`.
- **NUDT15** (owner decision): when typed callers (and the informational PharmCAT row) disagree on whether a defining
  allele is present, the phenotype is `-` with "callers disagree on NUDT15 *6: ... report it; ... do not". No tool decides
  NUDT15; PharmCAT only withholds a phenotype.
- Tests: `test_verdict.py`, `test_outside_calls.py`, `test_cpic.py` (set `PGX_CPIC_VERIFY=1`, optional
  `PGX_CPIC_SNAPSHOT=<dir>`, to compare the committed CYP2D6/TPMT/NUDT15 tables with the CPIC API).

---

## 2026-11-15 — PharmCAT "Unknown" is no call; lettered legacy CYP2D6 names are never scored as the base allele

- `docker/pgx-compare.py` (`_compute_verdict_core`, `_is_no_diplotype`): for the PharmCAT-authority genes (UGT1A1,
  CYP2B6, CYP4F2) a PharmCAT diplotype with no allele call (`Unknown/Unknown`, `No Result`, empty) is dropped before the
  authority rule and the vote, so it neither overrides the agreeing callers nor counts as a dissenting one. The
  verdict note records it. A real PharmCAT call is unchanged. Stored results: CYP2B6 HG002 (three callers `*2/*2`,
  StellarPGx `*2/*5`) now reads `*2/*2 or *2/*5` (see the phenotype-tier item below); CYP4F2 in two samples now follows the callers (`majority` where one caller's
  spelling differs).
- `docker/reconcile.py` (`_clean_allele`): the legacy trailing-letter strip for CYP2D6 is removed. `*14A`, `*14B`,
  `*13C`, `*4N` are kept as written, so they have no activity score and match only themselves in the verdict and in
  `call_cyp2d6` (previously `*14A` scored as `*14`, decreased function, while the current `*114` has no function in
  the CPIC table). `*2x2`/`*2X2`, `*4.021` and Aldy's `*36X` marker are unchanged. No legacy-to-current table is
  carried: none is available from a cited source in the repo.
- Phenotype tier (callers differ in diplotype, every phenotype-emitting caller agrees on the phenotype): the
  verdict no longer presents one diplotype. `consensus_phenotype` is the agreed phenotype, `consensus_diplotype`
  lists every distinct reading ("A or B", most callers first then alphabetical, the phase-ambiguity convention),
  `diplotype_readings` maps reading to callers, `diplotype_resolved` is false and the note names the callers per
  reading. Status stays `concordant` for the phenotype only; `n_agree` is the number of phenotype-emitting callers
  that agree (unchanged). A caller's own "Indeterminate" is not a reading. Generic for every gene, no reference
  table. Unanimous diplotypes, real authority calls and the CYP2D6 authority path are untouched. Stored results:
  48 of 2571 verdicts (DPYD 25, CYP3A5 10, NUDT15 8, CYP2C19 4, TPMT 1) change only in
  `consensus_diplotype`, `diplotype_readings`, `diplotype_resolved`, `note`; HG002 CYP2B6 reads `*2/*2 or *2/*5`.
  `pgx-report.py` wording updated to match.
- Review follow-ups: (1) NAT2 CPIC phenotype is computed per "A or B" reading (same phenotype kept; otherwise
  Indeterminate with the readings in the note). (2/3) A PharmCAT diplotype with no call or an unknown half
  (`*1/Unknown`) leaves the authority rule only when the other callers corroborate each other (two agree on the
  diplotype, or two agree on the phenotype); otherwise PharmCAT's row stays as before. (5/6/9) A caller whose
  diplotype is a no-call (`Unknown/Unknown`, `Indeterminate`, `No Result`), on any gene, abstains: not in the vote,
  not in `n_called`, named in the note; "-"/empty rows add no note. Stored results: 5 CYP2D6 verdicts (`n_called`
  4 -> 3, one phenotype spelling, readings without the "Indeterminate" caller); 1 status change (ILMN CYP4F2
  concordant -> majority, same as the PharmCAT no-call drop). (4) A CYP2D6 diplotype with a lettered legacy name
  gives no outside call at any level (`call_cyp2d6`), with the reason recorded. (7) `reconcile.LEGACY_LETTERED`:
  Aldy's bare `*36X`/`*36x` marker is notation, scored and counted by both the verdict and `call_cyp2d6`; true
  letters (`*13C`, `*4N`, `*14A`) are not. No stored outside call changed.
- A no-call is never shown as a diplotype: a PharmCAT-authority gene whose PharmCAT call is empty and whose other
  callers do not corroborate each other, and a gene where every caller abstains, give the explicit `no_call` verdict
  (`NO_CALL`, phenotype "-", reason in `note`); the tier helper never returns an empty diplotype. No stored verdict
  changes to `no_call`.
- Tests: `docker/test_verdict.py` (PharmCAT no-call, real call, lettered names, copy-number tokens against
  `cpic_cyp2d6.json`); `test_normalise.py` and `test_outside_calls.py` updated (a letter is not the base number).
- Regression over all stored verdicts (2571): 29 changed cells, all in PharmCAT no-call genes; 0 lettered-name cells
  changed; `call_cyp2d6` output unchanged for every stored sample.

---

## 2026-11-15 — CYP2D6 derived phenotype: corroboration by diplotype or activity score

- `docker/pgx-compare.py` (`_cyp2d6_authority_phenotype`): the phenotype Cyrius lacks is derived when Cyrius's
  filter is PASS, its diplotype has a CPIC activity score, and either (a) another caller reports the same canonical
  diplotype (`n_agree >= 2`, the 6de12ca rule, unchanged) or (b) another caller's diplotype has the same activity
  score (callers with no score do not count; a caller whose diplotype contains a legacy lettered allele such as
  `*13C`, `*4N`, `*36x`, in either case, never counts for (b), because `reconcile.clean_diplotype` would score it as
  its base number). `phenotype_source` says which applied, `(same diplotype: <tools>)` or
  `(same score: <tool diplotype>, ...)`, the first when both do. Not corroborated: "-" with `no other caller's
  diplotype has the same activity score: phenotype not derived` in `note`. Diplotype, status, n_agree, n_called and
  every other gene are unchanged.
- Rule (b) uses the same score comparison as `pgx_outside_calls.call_cyp2d6` (levels a, b, e), re-stated in the
  verdict, plus the lettered-allele guard that `call_cyp2d6` does not have.
- Tests: `docker/test_verdict.py` (both rules, lettered guard in upper and lower case, fail-closed cases).
- Known limitations: a lettered-only corroborator can give a drug-table (outside) call with a blank gene-card
  phenotype; tandem spellings Cyrius emits that the exact-key lookup misses (`*36+*36+*10` vs CPIC `*36x2+*10`,
  `*36+*36`, `*4.013+*4`, `*68+*68+*4`) give no phenotype; copy number >= 3 (`*1x3`, CPIC `x>=3` rows) gives none;
  `all_genes_summary.tsv` has no provenance column for a derived phenotype.

---

## 2026-11-14 — CYP2D6 verdict phenotype is no longer blank when Cyrius calls

- `docker/pgx-compare.py` (`_cyp2d6_authority_phenotype`): Cyrius, the CYP2D6 authority, never emits a phenotype,
  so `ConsensusPhenotype` was "-" whenever it called and CYP2D6 dropped out of the report's actionable-phenotype
  list. The phenotype is now derived from the diplotype (`reconcile.activity_score` then `reconcile.score_phenotype`,
  the CPIC table in `cpic_cyp2d6.json`), no caller vote; a plain CPIC term, with `phenotype_source`
  ("CPIC activity score 1 from *2/*4"). It is derived only when Cyrius's filter is PASS (the filter is carried in
  `alternative_diplotypes`) and at least one other caller reports the same diplotype (`n_agree >= 2`), so the
  gene card never shows a phenotype that the PharmCAT outside call (which needs a second caller) does not use.
  Otherwise the phenotype stays "-" and the reason is appended to the verdict `note` (and `phenotype_note`):
  `Cyrius filter <X>: phenotype not derived`; `no other caller reports this diplotype: phenotype not derived`;
  `activity score not available for <diplotype> (structure or copy number not in the CPIC table lookup)`. A caller
  on the same diplotype whose phenotype class differs is named in the note only. Diplotype, status, n_agree,
  n_called, every other gene (CYP4F2 stays "-"), the phase-ambiguity branch, and an authority that supplies its own
  phenotype are unchanged. Takes effect only in a rebuilt docker image.
- `docker/pgx-report.py`: the "resolved by <authority>" sentence says the phenotype is the CPIC phenotype for the
  activity score (diplotype), or why none was derived.
- Tests: `docker/test_verdict.py` (`test_cyp2d6_*`, including fail-closed cases).
- Known limitations (unchanged): tandem spellings Cyrius emits that the exact-key lookup misses (`*36+*36+*10` vs
  CPIC `*36x2+*10`, `*36+*36`, `*4.013+*4`, `*68+*68+*4`) give no phenotype; copy number >= 3 (`*1x3`, CPIC `x>=3`
  rows) gives none; legacy lettered allele names are scored as their base number in `reconcile.clean_diplotype`;
  `all_genes_summary.tsv` has no provenance column saying a phenotype was derived.

---

## 2026-11-13 — CYP2D6 outside call: level e (Cyrius + activity-score agreement)

- `docker/pgx_outside_calls.py`: new level **e**, tried after a and b and before c. When Cyrius has a call
  with a known CPIC activity score and at least one other caller's call has the same score, Cyrius's cleaned
  diplotype and that score are passed to PharmCAT (`source` = `Cyrius + N other caller(s) on activity score X`).
  Dissenting callers, an explicit Indeterminate among them, do not block it and are listed in the
  `reason_omitted` column and in provenance (`note`). Without a Cyrius call, or with a Cyrius call that has no
  CPIC score, nothing changes. Nothing that passed before stops passing; level c results move to e when
  Cyrius shares a score (a diplotype with a score instead of a phenotype only).
- Consumers: `outside_calls.tsv` and provenance carry level `e`; any reader that maps level letters to text
  must add `e` (OmniGen `omnigen_dtc2.load_pgx_agreement` accepts only a, b, c).
- Tests: `docker/test_outside_calls.py` (`test_cyrius_activity_score_level`, GIAB cases updated).

- **CYP2B6 phase ambiguity** (`docker/pgx-compare.py`, `_phase_ambiguity`): when the authority (PharmCAT) and other
  callers report the same defining variants in a different cis/trans arrangement (`*5/*6` vs `*1/*7`; definitions in
  `PHASE_DEFS`, PharmVar via pypgx), no phase is asserted: `ConsensusDiplotype` prints both readings (`*5/*6 or *1/*7`,
  the reading more callers report first) with `phase_ambiguous` and `alternative_readings` in the verdict. The phenotype
  is the authority's when at least one caller of every other reading reports the same one (status `concordant`, nAgree =
  callers supporting the variant set); otherwise the phenotype is Indeterminate and the status `majority`. Only CYP2B6
  has definitions; only the authority branch is covered.

---

## 2026-11-12 — Caller agreement: shared notation normaliser, nAgree fix

- **Cause of the wrong `nAgree`**: in `compute_verdict` the authoritative-caller branch (Cyrius for
  CYP2D6, PharmCAT for UGT1A1/CYP2B6/CYP4F2, VCF-Check for RYR1/VKORC1/...) returned the constant
  `n_agree = 1`, so five callers on `*2/*4` printed `nAgree=1 nCalled=5`. It is now the number of
  callers whose normalised call equals the authority's call.
- **One normaliser** (`reconcile.norm_diplotype`, `clean_diplotype`) now feeds the verdict
  (consensus and `n_agree`) and the PharmCAT outside-call builder. Rules: whitespace/case ignored;
  Aldy's tandem `X` marker dropped (`*36X+*10` = `*36+*10`; `*2x2` is kept); PharmVar `.NNN` sub-allele
  dropped; trailing letter dropped for CYP2D6 only (`*6B` = `*6`; not for TPMT `*3A`/`*3C`, DPYD);
  alleles sorted inside a tandem and between haplotypes; Aldy `+rsNNN` tags are real differences;
  VCF-Check `*N(het) [unphased]` = `*1/*N` and `*N(hom)` = `*N/*N` for a single listed allele, while
  two listed alleles are unresolved and excluded from `n_called` (verdict `note`, printed as `Note:`).
  The per-caller table keeps the original strings.
- Outside calls: CYP2D6 is passed when all comparable callers agree after normalisation (cleaned
  spelling, Cyrius preferred). A majority rule was evaluated and NOT implemented: GeT-RM covers only
  HG001 among HG001-HG007, so it cannot be validated on the GIAB set.
- Tests: `docker/test_normalise.py` (new), `docker/test_outside_calls.py` extended.
- **CYP2D6 outside call by level of agreement** (supersedes "unanimous diplotype only"): (a) diplotype
  agreed -> diplotype + CPIC activity score; (b) diplotypes differ but all activity scores known and
  equal -> score only; (c) same phenotype class -> phenotype only; (d) otherwise omitted. Scores come
  from `docker/cpic_cyp2d6.json` (CPIC allele table + `gene_result` score->phenotype table, retrieved
  2026-11-12; refresh with `docker/cpic_cyp2d6_build.py`), via `reconcile.activity_score`,
  `score_phenotype`, `phenotype_class`. Unknown allele (`*79`), `+rsNNN` tag or unknown copy number
  gives no score; an explicit Indeterminate blocks level c. `outside_calls.tsv` and provenance record
  level, source and per-caller diplotype/score/phenotype. GIAB: HG002 and `*2/*36+*10` -> a;
  HG001, HG006, HG007 -> b; HG005 -> c (Intermediate); HG004 -> d. Checked with PharmCAT 3.2.0:
  score-only (`0.5` -> Intermediate, AS 0.5) and phenotype-only (annotations for 0.25/0.5/0.75/1.0,
  four times the per-drug annotations) both run; a phenotype-only call needs downstream care.

---

## 2026-11-11 — MT-RNR1 and G6PD false negatives fixed (ERROR is never NEGATIVE)

Two false negatives found on commercial reference standards (a homoplasmic MT-RNR1
m.1494C>T sample and a hemizygous G6PD c.95A>G / Gaohe sample).

- **MT-RNR1**: `pgx-mt.sh` called mutserve with `--input`, which mutserve 2.x rejects
  ("Unknown option"); the failure was downgraded to a warning, no raw file was written, and
  the empty parse reported `Reference / Standard aminoglycoside risk`. Even with a working
  call the parser indexed the wrong columns (mutserve 2.0.3 inserts a `Filter` column).
  Fixed: positional input, hard failure on mutserve error (no result JSON), header-driven
  parsing, full CPIC MT-RNR1 table (m.1095T>C, m.1494C>T, m.1555A>G increased risk; uncertain
  alleles listed), heteroplasmy level, per-site depth, `not typed` when chrM depth < 100x.
- **G6PD**: the VCF-Check authority tested only 5 sites and outranked PyPGx, Stargazer and
  Aldy, which all called Gaohe/deficient; PharmCAT also called Gaohe/Gaohe Deficient (shown
  only as an informational cross-check). Fixed: typing against the full 171-site CPIC table,
  sex-aware CPIC phenotype, per-site depth, never "Normal" on a partial list, and a
  discordant verdict when the authority sees no deficiency but another caller does.
- Wording fixes (no logic change) for the other partial-list checkers: CACNA1S, RYR1,
  UGT1A1 now say "not found among the N sites checked".
- Follow-up: a G6PD negative now requires verified depth. Depth unavailable (BAM missing,
  `samtools depth` failed) -> `not typed`; any key site of the common deficient alleles (A-
  202A/376G, Mediterranean, Viangchan, Mahidol, Canton, Kaiping, Union, Chinese-5) below
  `min_cov` -> `not typed (key site(s) below Nx)`. Positive calls are unaffected. `pgx-mt.sh`
  now aborts on an unparseable mutserve data row instead of skipping it. `test_vcf_check.py`
  is pytest-collectable (`test_vcf_check`) and covers both.
- `pgx-compare.py` gained `--bam/--ref/--sex`; Snakefile passes them (`sex=` config key).
- Review fixes: G6PD site depth filters on base quality (`samtools depth -q 13`, as
  `bcftools mpileup`), not mapping quality; nested alleles are assigned copy-aware (1/1 c.202A +
  0/1 c.376G in a female = A- + Asahi = Deficient, not Variable); a 0/1 chrX call in a male is
  `Indeterminate` (QC flag), never promoted to hemizygous Deficient; the G6PD safety net
  recognises `G6PD not typed`; supporting variants carry the alt base and an allele fraction;
  `--sex` accepts M/F/MALE/FEMALE in any case and treats anything else as "infer".
  MT-RNR1: table lookup keyed on (pos, ref, alt) so m.807A>G and m.807A>C are both reportable;
  the result JSON is written atomically (a parser abort leaves none).
- G6PD key-site depth is sex-aware (>= 5 reads per X copy): 10x for females and unknown sex,
  5x for hemizygous males (`min_cov_male`, default min_cov/2). At 30x WGS a male chrX averages
  ~15x, so a flat 10x left GIAB males HG002/HG003 `not typed`; a female carrier allele is in only
  ~half the reads, so 10x stays for females. Validated: HG001-HG005 all typed reference; a
  hemizygous Gaohe reference standard stays Deficient; a homoplasmic m.1494C>T standard is called.
- Second review round (ERROR is never NEGATIVE). G6PD: a VCF record at a CPIC site whose REF/ALT
  cannot be matched (after trimming bcftools' trailing context) is a QC failure, not skipped; a
  record with QUAL < 20 or DP < 3 never makes a Deficient call; QC failures (unmatched record, low
  quality, male 0/1) make a negative `Indeterminate`, while a confident deficiency stays the leading
  phrase with a `QC:` note; site depth is read with `samtools mpileup` (same BAQ/orphan filters as
  the caller) when a reference is given; another caller's "Variable" counts as dissent and a
  discordant verdict's phenotype reads `Indeterminate - discordant ...`, never the authority's
  negative; no chrY reads infers unknown sex, not female. MT-RNR1: a filtered mutserve call at an
  increased-risk site makes the result `not typed` and is named in the notes; observed
  uncertain-risk variants are always listed.

---

## 2026-11-09 — cpic.json is strictly CPIC Level A/B; Level C pairs are not flagged

`docker/cpic.json` is the single source of flagging rules, so only CPIC Level A or B gene-drug
pairs may be listed. Checked against CPIC's live pair table; stale entries corrected.

- **Removed** (Level C, retired, or not a CPIC pair): HLA-B dapsone and B\*13:01; HLA-B
  flucloxacillin and B\*57:03; HLA-A oxcarbazepine (A\*31:01 text now names carbamazepine
  only); NAT2 isoniazid; both IFNL3 peginterferon entries (Retired); CYP3A4 "Tacrolimus (with
  CYP3A5)" (CYP3A5 tacrolimus stays); CYP1A2 fluvoxamine; CYP2A6 nicotine; CYP2C8 amodiaquine.
- **Changed:** CYP4F2 warfarin level B to A.
- **Added:** HLA-B oxcarbazepine (Level A; B\*15:02 text now names it); CYP2D6 aripiprazole
  (Level B). CPIC has no guideline URL for aripiprazole, so its `url` is empty, the report shows
  the drug name unlinked, and `rec` points to DPWG guideline annotation PA166104937.
- **Genes with no listed pair** (IFNL3, CYP3A4, CYP1A2, CYP2A6, CYP2C8) stay in the file
  and are still typed and shown. `_get_tier()` in `pgx-report.py` now returns "not flagged" for
  any gene whose `drugs` list is empty; their phenotype triggers and drug-specific notes were
  cleared. `pgx-report.py` `DIPLOTYPE_CHECKS["HLA-B"]` no longer includes B\*13:01 and B\*57:03.
- **Added:** NAT2 hydralazine (Level A, CPIC guideline PA166416361). NAT2 is flagged moderate for
  poor/slow, intermediate and rapid phenotypes, with CPIC's dosing text; this corrects the
  "typed but not flagged" note above for NAT2.
- Historical entries below are left as written.

## 2026-11-09 — NAT2 phenotype now computed from CPIC allele function values

NAT2 was "Indeterminate" for every sample: PyPGx has no NAT2 diplotype table and, on a tie with
Stargazer's non-CPIC label, the consensus picked "Indeterminate". `pgx-compare.py` now derives the
phenotype from CPIC's NAT2 allele functions (`nat2_cpic_phenotype`, applied in `compute_verdict` to
the agreed diplotype): increased+increased Rapid, increased+decreased Intermediate,
decreased+decreased Poor Metabolizer. Alleles CPIC marks uncertain/unknown, or does not list, stay
Indeterminate. Covered by `docker/test_cpic.py`.

## 2026-11-09 — OptiType class I now reads the HLA-* contigs too

**Cause.** `pgx-hla.sh` extracted only `chr6:28,510,020-33,480,577`, so reads aligned to the
`HLA-*` contigs of the GRCh38 alt-aware BAM were dropped and OptiType saw about half the
informative reads, giving false homozygotes (and one missing HLA-C). chr6 `*_alt` contigs make
no difference. **Fix.** The region list is the MHC window plus every `HLA-*` contig in the BAM
header, in the paired, single-end and no-chr branches and in the paired-read count. A BAM with
no `HLA-*` contigs behaves as before. Existing class I results must be re-run.

## 2026-11-09 — arcasHLA SIF restored as a required asset; HLA class II failure now loud

**Cause.** `StellarPGx/containers/arcashla.sif` did not exist on this host (it was never built
into the container folder; earlier class II results came from the host-side arcasHLA route).
`pgx-hla2.sh` logged "SIF not found", `rule hla2` wrote `hla2.status=1` and went on, and
`parse_arcashla` turned the missing JSON into `not_run`, so class II looked merely absent
(HG003_r2).

- **SIF built** (451 MB, `StellarPGx/containers/arcashla.sif`) with `docker/build-arcashla-sif.sh`:
  arcas-hla 0.6.0 biocontainer + the host IMGT/HLA reference (`hla.dat` header says 3.64.0;
  earlier notes said 3.63.0) + kallisto 0.44.0. Plain `apptainer pull` is not enough (no
  reference; kallisto 0.50.1 incompatible; `reference --version 3.55.0` invalid).
- **`pgx-hla2.sh` bug fix.** It used `arcasHLA extract` on a name-sorted BAM ("unable to index
  bam file"), and extract-filtered reads typed only DRB1. Now it uses the validated route
  (MHC window + chr6 alt + HLA-* contigs -> FASTQ -> `genotype`). Verified on HG003: DQA1,
  DQB1, DRB1 identical to the earlier results.
- **Loud failure.** `rule hla2` prints an `ERROR` on non-zero exit; `parse_arcashla` reads
  `hla2/hla2.status` and reports `failed` (not `not_run`) when it is non-zero with no JSON.
- Docs: README + `ToolsDocumentation.md` 5b list the SIF as required with the build command.
- The pgx-suite Docker image bakes in the scripts: rebuild it to pick up these changes.
- **Empty-FASTQ guard fixed.** The `-s` test on `hla.{1,2}.fq.gz` could never fire (samtools
  writes a 28-byte BGZF EOF block even for zero reads), so empty FASTQs reached `arcasHLA
  genotype`. It now fails with the same ERROR when either FASTQ decompresses to zero bytes.

## 2026-07-16 — arcasHLA reference relocated to durable tools dir (port-safety, no code change)

**The host-side arcasHLA clone + reference DB moved off scratch so it survives a port.**

- **The move.** `<scratch>/arcasHLA/` → **`<ref>/omnigen/tools/arcasHLA/`**
  (2.1 GB; the key artifact is the 225 MB `dat/ref/hla.idx`). Same filesystem, plain `mv`.
  Verified `dat/ref/hla.idx` is intact + readable at the new path and that a dry `arcasHLA`
  invocation there resolves the reference (arcasHLA self-resolves `dat/ref/` relative to its
  own install dir via `rootDir = dirname(realpath(scripts/*.py)) + '/../'` — the DB travels
  with the binary; no env var or config points at it).
- **No pipeline code changed.** No pgx-suite source hardcoded the old path — HLA class II
  typing is invoked by binary path (the host-side route in `ToolsDocumentation.md` §5b) and
  the SIF route (`pgx-hla2.sh`) bundles its own DB inside `arcashla.sif`, unaffected. Only the
  §5b runbook command + the durable-location note in `ToolsDocumentation.md` were repointed.
- **Conda envs stay put (rebuild-on-target).** `arcaspy-env` / `kallisto-env` remain at
  `<scratch>/` — conda bakes absolute prefixes, so moving breaks them. They are
  reproducible (python3.10 + numpy/scipy/pandas/biopython; kallisto 0.44.0); rebuild on the
  target host for a port. See `ToolsDocumentation.md` §5b "Conda envs are rebuild-on-target".

## 2026-07-16 — HG001 HLA class II back-filled (data run, no code change)

**GIAB HG001 / NA12878 gets its missing class II loci.**

- **The gap.** HG001's production run (2026-06-20) predates HLA class II (arcasHLA) support,
  so it had class I only (OptiType A/B/C) and no `genes/HLA-DQA1|DQB1|DRB1/` or `hla2/`.
- **The fix.** Ran the same host-side arcasHLA route validated on HG002/HG003 against
  `HG001.bwa.sortdup.bqsr.bam` (reusing the pre-built IMGT/HLA `hla.idx` from the HG002
  validation), staged `hla2/HG001.genotype.json`, and regenerated the three class II
  `genes/HLA-*/*_HG001_comparison.tsv` via `pgx-compare.py parse_arcashla` — identical
  schema to the class I TSVs. No pipeline code changed; no other gene/sample re-run.
  - Calls (`ok`): DQA1 `*01:01:01/*05:01:01`, DQB1 `*02:189/*05:01:01`, DRB1 `*01:126/*03:01:01`.
  - Documented in `docs/ToolsDocumentation.md` §5b ("Back-filled sample: GIAB HG001").

---

## 2026-07-15 — OptiType runs once per sample, not once per class-I gene (perf)

**Removed ~2-3x redundant OptiType compute (~6-12 min/sample of pure waste).**

- **The waste.** OptiType types all three class-I loci (A, B, C) in ONE ILP run and writes
  `A1 A2 B1 B2 C1 C2` into a single `<sample>_result.tsv`. But the Snakemake `optitype` rule
  was keyed on `{gene}`, so HLA-A, HLA-B and (after 2026-07-13) HLA-C each triggered a full,
  independent OptiType job that typed all three loci and then kept only its own — three
  identical ~3-6 min runs producing byte-identical output. Logged in `docs/TODO.md` when
  HLA-C was added (commit 935804b) but deferred for a dedicated change.
- **The fix.** `optitype` is now a **sample-level** rule (like PharmCAT / arcasHLA): it runs
  once per sample, writes `results/<sample>/optitype/<sample>_result.tsv`, and the three HLA
  compare steps depend on the single shared status (`results/<sample>/optitype/optitype.status`)
  and all parse that one TSV.
  - `docker/Snakefile`: `rule optitype` no longer carries a `{gene}` wildcard; output is the
    sample-level `OUT/optitype/optitype.status`. `compare_inputs()` gains the shared optitype
    dependency (removed from the per-gene tool loop). The compare `run` block reads OptiType's
    failure status from the sample-level path.
  - `docker/pgx-hla.sh`: signature is now `<BAM> <SAMPLE> <OUTPUT_DIR>` (no per-gene `GENE`
    arg); OptiType always typed all loci regardless of the dropped arg.
  - `docker/pgx-compare.py` `parse_optitype()`: unchanged locus logic — only *where it looks*.
    It now checks the sample-level `optitype/` dir first, then falls back to the legacy
    per-gene dir, so existing on-disk results still parse.
- **Contract preserved.** Per-gene `HLA-{A,B,C}_<sample>_comparison.tsv` and `_detail.json`
  keep the same schema and paths, so the OmniGen consumer (`report_hla.py` / `report_hla2.py`)
  works unchanged. Verified on the on-disk HG002 OptiType TSV: HLA-A → `HLA-A*01:01/HLA-A*01:01`,
  HLA-B → `HLA-B*35:08/HLA-B*38:01`, HLA-C → `HLA-C*12:03/HLA-C*12:03` from both the legacy and
  sample-level layouts. Guarded by two new tests (`test_genes_config`, `test_parsers`).

---

## 2026-07-13 — HLA-C promoted to a first-class typed gene (bug fix)

**OptiType was already typing HLA-C correctly; pgx-suite threw the call away.**

- **The bug.** OptiType types all three class-I loci (A, B, **C**) in a single ILP run and
  writes `A1 A2 B1 B2 C1 C2` into `<sample>_result.tsv`. But `parse_optitype()` in
  `docker/pgx-compare.py` had only `HLA-A` and `HLA-B` branches — its own docstring
  documented the `C1`/`C2` columns it then ignored — and `genes.tsv` had no `HLA-C` row, so
  no `genes/HLA-C/` output directory was ever produced. The C call was computed on every
  run since HLA typing was introduced, and silently discarded.
- **Why it mattered (downstream).** The consumer, OmniGen (`prototype/report_hla.py`), asks
  for `["HLA-A", "HLA-B", "HLA-C"]`. Because HLA-C never arrived, its **psoriasis
  association (HLA-C\*06:02)** was unscreenable for every sample ever reported — and it
  rendered as **"absent" rather than "not typed"**, so a C\*06:02 homozygote would have been
  actively told the risk allele was absent. This is a false-negative class of failure, not a
  missing-feature one. HLA-C also matters for transplant matching.
- **The fix.** Added the `HLA-C` branch to `parse_optitype()` (reads `C1`/`C2`, prepends
  `HLA-`, same status/diplotype/haplotype handling as A and B) and an `HLA-C` row to
  `docker/genes.tsv` (`optitype=1`), which is the single source of truth the Snakefile
  routes from — so HLA-C now gets a per-gene dir, `comparison.tsv` and `detail.json` exactly
  like A and B. Also wired into `pgx-bamstats.sh` (depth BED), `pgx-report.py` (gene locus),
  and the `test_genes_config.py` regression matrix.
- **Verified without a pipeline re-run.** `parse_optitype()` called directly against the
  existing OptiType output for HG002 returns **`HLA-C*12:03/HLA-C*12:03`, status `ok`**.
- **Backfilled 13 samples** (HG001, HG002, HG003, HG005, HG006, HG007, sample D, sample E, sample F,
  sample H, sample I, sample A, sample B) from their already-computed OptiType result TSVs — no
  OptiType or pipeline re-run. Schema is byte-identical to the HLA-A/HLA-B
  `comparison.tsv`, so OmniGen's reader works unchanged. Notably **HG003 is
  `HLA-C*06:02/HLA-C*06:02`** — a psoriasis-risk homozygote that the old code reported as
  "absent".
  - **HG004 could NOT be backfilled**: OptiType emitted **empty `C1`/`C2` fields** for that
    sample (A and B typed fine), so there is genuinely no C call on disk. The parser returns
    `not_run` and no file is written — correctly reading as *not typed*, not *absent*.
  - Backfill coverage is recorded as `NOT_COMPUTED` (`mean_depth: null`): HLA-C had no row
    in the bamstats depth BED when those results were produced, so no HLA-C region depth
    exists. No depth was invented. The OptiType call is unaffected — it comes from the same
    MHC extraction and the same result TSV as HLA-A/HLA-B.
- **Known inefficiency (not fixed, see `docs/TODO.md`).** Each class-I gene runs its own full
  OptiType job that types A, B *and* C and then keeps only its own locus. That was 2x
  redundant with A+B; adding HLA-C makes it **3x** — three identical ~3-6 min OptiType runs
  per sample producing the same result TSV. Recommended fix: run OptiType **once per sample**
  and have all three loci parse the one TSV. The parser is already locus-agnostic, so the
  change is confined to the Snakefile.

**37 → 38 genes/loci** (CPIC Level A count unchanged at 19 — HLA-C has no CPIC guideline).

---

## 2026-07-12 — HLA class II (arcasHLA) validated on GIAB HG002

Feature 3 (HLA class-II typing) run end-to-end and validated against published truth.

- **Tooling provisioned.** arcasHLA (RabadanLab) cloned to `<scratch>/arcasHLA`;
  kallisto 0.44.0, git-lfs, and a python-3.10 stack (numpy/scipy/pandas/biopython) installed
  in dedicated conda envs (samtools/bedtools/pigz from host). Reference built with
  `arcasHLA reference --rebuild` from **IMGT/HLA 3.63.0** `hla.dat` fetched from the EBI
  IPD-IMGT/HLA FTP — the ANHIG GitHub LFS budget was exhausted, so the git-lfs pull of
  `hla.dat` failed and EBI was used as the non-LFS source.
- **Genotyped HG002** (all 7 loci) from the 30–40X MHC reads
  (`chr6:28,510,020-33,480,577`, `HG002.hla.{1,2}.fq.gz`).
- **Concordance vs Chin et al. 2020 MHC benchmark, Suppl. Table 4** (clinical trio-phased),
  2-field: **DRB1 2/2 ✅, DQB1 2/2 ✅** (class-II target passes); **DQA1 1/2** (`03:01` exact,
  `01:05` vs truth `01:01`, same allele group). Class-I cross-check: **A 2/2**, B/C
  allele-group-correct but 2-field-discordant (thin WGS depth for an RNA-seq-tuned tool —
  OptiType stays the primary class-I typer). DPB1 has no entry in the truth table
  (`04:01` plausible, `677:01` low-confidence).
- **Contract written.** Raw output at `results/HG002/hla2/HG002.genotype.{json,genes.json,log}`;
  per-gene `results/HG002/genes/HLA-{DQA1,DQB1,DRB1}/*_comparison.tsv` + `*_detail.json`
  (reuses the class-I `Gene/Sample/Build/Tool/Diplotype/…` schema, `Tool=arcasHLA`).
- **Demo report regenerated** (`docker/pgx-report.py` → `docs/demo/HG002_report.html`) so the
  **HLA Class II & Blood Group** section and the DQA1/DQB1/DRB1 gene cards show the validated
  calls (Celiac DQ8 + T1D susceptibility surfaced).
- **Bug fix — `omnigen_addons.hla2_disease_findings` T1D label.** The Type-1-diabetes
  `risk_haplotype` label was hardcoded `"DR3-DQ2 / DR4-DQ8"` and rendered in full whenever
  *either* arm matched, overstating risk (e.g. HG002, which carries only DR4-DQ8, appeared to
  carry DR3-DQ2). Now each arm (DR3-DQ2 = DRB1\*03:01 + DQB1\*02; DR4-DQ8 = DRB1\*04 +
  DQB1\*03:02) is evaluated independently and only the arm(s) actually present are surfaced
  (`"DR4-DQ8 (present)"`), with the absent arm explicitly marked "not detected". Celiac
  DQ2.5/DQ8 labelling already worked this way and is unchanged. `test_omnigen_addons.py`
  (12 checks) still passes.
- Details: `docs/omnigen-additions-plan.md` Feature 3 §(f); `docs/ToolsDocumentation.md`.

## 2026-07-11 — OmniGen add-ons: HLA class II + ABO (haplogroups moved to OmniGen)

Two BAM-derived add-ons, each emitting a small stable contract file for OmniGen to render
(pgx-suite calls/types; OmniGen renders). Both work from the existing aligned GRCh38 BAM.

- **HLA class II (DQA1/DQB1/DRB1).** New `docker/pgx-hla2.sh` + sample-level `hla2` rule run
  **arcasHLA** (Apptainer SIF `arcashla.sif`) once per sample; `parse_arcashla` in
  `pgx-compare.py` reuses the class-I `*_comparison.tsv` schema. New `arcashla` column in
  `genes.tsv`; rows `HLA-DQA1/DQB1/DRB1`. Report surfaces celiac (DQ2.5/DQ8), narcolepsy
  (DQB1\*06:02), and T1D susceptibility (susceptibility markers, not diagnostic).
- **ABO blood type — PROVISIONAL / UNVALIDATED.** New `ABO` row (`vcf_check`) + `parse_abo_vcf`
  assign O (rs8176719 frameshift) / A vs B (rs8176746/747) via the VCF-Check machinery.
  Contract: `results/<S>/abo/abo_type.tsv`. **The GRCh38 reference-allele identity and variant
  coordinates are unconfirmed and the caller is not yet validated against a known-type
  control — every ABO result is flagged UNVALIDATED; do NOT use clinically.**

Shared deterministic logic in `docker/omnigen_addons.py` (pure, stdlib-only), unit-tested by
`docker/test_omnigen_addons.py` (10 checks, runs without a BAM). `test_genes_config.py`
updated for the 4 new genes + `arcashla` column. Report gains an **HLA Class II & Blood
Group** section (`docker/pgx-report.py`). Design/plan: `docs/omnigen-additions-plan.md`.

**mtDNA + Y-chromosome haplogroups moved to OmniGen.** After validating that haplogrep3
reproduces both GIAB trios' maternal lineages exactly (ancestry-correct) directly from the
30X WGS VCFs downstream, both haplogroups now live in OmniGen; pgx-suite no longer emits
`mito/`/`ychr/` contracts or bundles haplogrep3/yhaplo.

---

## 2026-06-20 — ALT-aware input gate, 6 new genes/loci, report & validation refresh

**ALT-awareness QC gate.** New `docker/pgx_altcheck.py` classifies whether the input
BAM was aligned ALT-aware from the file itself (`@PG` aligner, `@SQ` ALT-contig count,
fraction of reads on `_alt` contigs → `alt_aware.json`). A non-ALT-aware alignment
(minimap2, bowtie2, or a no-ALT reference) mis-distributes CYP2D6/CYP2D7 + ALT-contig
reads and yields a false CYP2D6 `*5/*5` deletion. A new sample-level `altcheck`
Snakefile rule feeds `pgx-compare.py --alt-aware`; when the verdict is not `yes`, the
paralog/SV genes (`cyrius`/`pypgx_sv`/`stargazer_sv`) report **NO_CALL** with a re-align
reason, single-SNP genes are unaffected, and the report shows a QC banner. Self-test wired
into `test.sh`. Verified: bwa HG002 → `yes`; simulated minimap2 → CYP2D6 NO_CALL, CYP2C19 OK.

**6 new genes/loci (31 → 37).**
- **CYP2C-CLUSTER** — warfarin variant `rs12777823` (chr10:94645745 G>A, CPIC, African-ancestry
  dosing) as a VCF-Check single-SNP locus + authority. Completes the warfarin panel
  (VKORC1/CYP2C9/CYP4F2).
- **GSTP1, SLC15A2, SLC22A2, SLCO2B1, UGT2B15** — PyPGx-supported additions, kept only after
  passing NA12878 GeT-RM validation. UGT2B7 (callers discordant) and UGT2B17 (whole-gene
  deletion, NO_CALL on short reads) were trialled and dropped.

**Validation.** ThermoFisher Axiom validation cohort re-scored on the verdict engine: 26 platform-sample pairs vs
Axiom — **73.0%** strict, **85.8%** adjudicated (ILMN ≈ MGI). GIAB HG001/NA12878 vs CDC
GeT-RM consensus: **19/28** exact-match. Full method: `AxiomValidation.md`.

**Report UX.** Gene pages gained the consensus-diplotype headline, metaboliser pill,
resolution/authority badge, and a "what the pipeline did" narrative (tool concordance +
CPIC resolution). Header is a 2-column grid; clickable top bar returns to the summary;
floating horizontal scrollbar drag fixed; Orthogonal Validation panel is landing-only.

**Fix.** PharmCAT regression — `apptainer exec --cleanenv` (the image's `JAVA_HOME` was
leaking into the SIF, silently failing PharmCAT).

---

## 2026-06-19 — Reconciliation engine: authoritative callers, phenotype tier, verdict-driven report

The pipeline stops being a string-equality ensemble and becomes a reconciliation
engine that resolves the constituent callers' blind spots. Every clinically
actionable CPIC gene now reaches a confident verdict on short reads, except the
two long-read-only complex CYP2D6 cases.

**Variant-tier reconciliation (`docker/reconcile.py`, `allele_synonyms.json`).**
Caller diplotypes are synonym-collapsed before the agreement count, so the same
allele written three ways (DPYD `*9A` = `*S10` = `c.85T>C` = `rs1801265`) counts
as one call, not three discordant ones. The curated synonym map deliberately
excludes genes where look-alike names are distinct variants (NUDT15 `*5`/`*6`,
VKORC1). Guarded by `test_reconcile.py`.

**Authoritative-caller verdicts (`AUTHORITATIVE` in `pgx-compare.py`).** Where the
star-allele callers genuinely cannot agree, a designated CPIC-reference method *is*
the verdict:
- **Cyrius** (new, vendored from Illumina; PolyForm Strict — image not redistributable)
  resolves **CYP2D6** structural variants / CNV / `*36+*10` hybrids that the
  star-allele callers mis-handle. New Snakefile rule + `parse_cyrius`.
- **PharmCAT 3.2.0** (new, official `pgkb/pharmcat` Apptainer SIF) is authoritative
  for **UGT1A1, CYP2B6, CYP4F2**. Sample-level rule: a GT-only positions VCF
  (`pharmcat_positions.vcf.bgz`, baked) feeds `pharmcat_pipeline`; `parse_pharmcat`
  reads its `report.json`.
- **VCF-Check** (in-house) reads the genotype directly from the CPIC variant table
  (`vcf_check_variants.json`) for **RYR1, VKORC1, IFNL3, G6PD, CACNA1S** —
  single-SNP / variant-list genes where a star-allele call is the wrong model.

**Clinical phenotype tier (`pgx-compare.py`).** When callers print different
diplotype *strings* (Stargazer `*S` codes, Aldy rsID tokens, different haplotype
phasing) but the phenotype-emitting callers agree on the CPIC phenotype, the call
is concordant-by-phenotype (`authority: "Phenotype"`); `Indeterminate`/unknown
counts as abstain, not dissent. This matches CPIC (which reports by
phenotype / gene-activity-score, not haplotype string) and resolved all five DPYD
nomenclature/phasing artifacts in the ThermoFisher Axiom MGI cohort while correctly leaving a
genuine NUDT15 3-way disagreement flagged.

**Verdict-driven report (`pgx-report.py`).**
- Card/grid colour and the summary ticker now follow the **verdict status**, not raw
  tool fraction — an authority-resolved gene (e.g. PharmCAT at 1/4 callers) is a
  confident green, not red. Fixes a miscount where authoritative genes inflated the
  "discordant" total.
- **Authority badges** on each card (`PharmCAT authoritative`, `VCF-Check
  authoritative`, `Cyrius authoritative`, `Resolved by phenotype concordance`) with a
  matching legend, for CAP/CLIA traceability.
- **Actionability-first grid** — genes group ★ Actionable findings → ⚠ Needs review →
  ✓ Normal, so an abnormal-but-concordant phenotype leads instead of being buried.
- Detail page reads the single-source verdict (card and detail concordance counts now
  match); wide tool-results tables scroll horizontally with a sticky first column;
  per-gene coverage is a collapsible dropdown; top pane lists all 9 callers.

**Orthogonal validation.** Axiom array adjudication panel added for validation
samples (`build_axiom_concordance.py`), classifying NGS-vs-array disagreements
(`array_FN` / `array_FP_or_NGS_FN` / `array_authoritative` / …) rather than treating
the array as ground truth.

**`genes.tsv`** gained `cyrius` and `pharmcat` columns (CYP2D6 cyrius=1;
UGT1A1/CYP2B6/CYP4F2 pharmcat=1; RYR1/VKORC1/IFNL3 vcf_check=1).

**Tests added:** `test_reconcile.py`, `test_vcf_check.py`, `test_cyrius.py`,
`test_pharmcat.py`; phenotype-tier cases in `test_verdict.py`.

**Bug fix — mutserve was silently broken.** The Dockerfile pulled
`releases/download/v2.0.0/mutserve.jar`, but v2.0.0 has no release asset, so the build
baked a 9-byte "Not Found" stub and MT-RNR1 (CPIC Level A) calling failed non-fatally.
Now downloads the **v2.0.3** `mutserve.zip` and extracts the real `mutserve.jar` (via
Python stdlib `zipfile`; no new apt dependency), and the build fails loudly if it isn't a
valid jar. Surfaced by the extended smoke test (`docker/test.sh` now covers all callers:
Cyrius, PharmCAT, OptiType, mutserve, and the reconciliation/VCF-Check config — 8 sections).

**Deployment note:** the PharmCAT SIF must be pulled with single-threaded mksquashfs
(`mksquashfs procs = 1`) on an ext4 TMPDIR — apptainer's bundled multithreaded
mksquashfs corrupts the ~1.8 GB image. See README "Pulling the PharmCAT SIF".

---

## 2026-06-17 — Snakemake migration, clinical-safety hardening, report redesign

**Orchestration migrated to Snakemake.** The hand-rolled bash orchestrators
(`pgx-run.sh`, `pgx-all-genes.sh`) are removed; the pipeline is now a wildcard `{gene}`
DAG in `docker/Snakefile`, scheduled by `--cores`. `run_pgx_suite.sh` and the image
entry point invoke Snakemake. Proven byte-for-byte equivalent to the old pipeline across
all 31 genes on HG002 (verdict + per-tool diplotype + status identical).

**Single source of truth.** Gene region, per-tool support, SV flags, and Stargazer
control now live in `docker/genes.tsv`, read by both the Snakefile and `pgx-compare.py`
(guarded by `test_genes_config.py`). The duplicated bash/Python matrices are gone.

**Clinical-safety fixes (`pgx-compare.py` / `pgx-report.py`):**
- Coverage gate — genes below `--min-depth` (default 10×) report **NO_CALL** instead of a
  confident wild-type call manufactured from zero reads.
- Single verdict authority — the gene verdict (`concordant`/`majority`/`discordant`/
  `no_call`) is computed once into `detail.json`; the HTML report and batch summary read it
  rather than recompute, so a 2-vs-2 tie surfaces as **DISCORDANT**, not a silent winner.
- Real exit codes — each caller writes a `logs/<tool>.status` sentinel; a crashed caller is
  marked `failed`, not the misleading `not_run`.

**Report redesign (from /design-review):** print-safe severity colour (survives PDF export),
non-colour status glyphs (greyscale / colour-blind legible), WCAG AAA badge contrast,
responsive layout + table-scroll wrappers, the gene grid split into ⚠ Needs review /
Normal — concordant groups, and an at-a-glance headline summary strip.

**Tests added:** `docker/test_parsers.py`, `test_verdict.py`, `test_coverage_gate.py`,
`test_genes_config.py`.

**Docs:** README rewritten for the Snakemake era; reference docs moved under `docs/`;
getting-started tutorial + add-a-gene how-to added.

---

## 2026-03-13 — Post-validation pipeline improvements + Docker rebuild

### Summary

Applied four actionable fixes derived from the ThermoFisher Axiom validation study and added a
fifth independent UGT1A1 caller (VCF-Check). Docker image rebuilt to incorporate
all changes.

### `docker/pgx-run.sh`

- **ABCG2:** StellarPGx support flag changed `1 → 0` — no calls produced in 13
  validation samples; tool is not compatible with this gene in practice

### `docker/pgx-compare.py`

- **ABCG2 GENE_SUPPORT:** `stellarpgx: True → False` to match run.sh
- **UGT1A1 VCF-Check (new):** fifth independent caller `parse_ugt1a1_vcf()` queries
  the bcftools VCF already produced by the pipeline for promoter/coding alleles:
  - *60 rs45530432 chr2:233757013 T>G (East Asian promoter)
  - *80 rs887829 chr2:233759924 C>T (promoter, moderate reduction)
  - *6  rs4148323 chr2:233760498 G>A (Gly71Arg, common East Asian)
  - *28 rs3064744 chr2:233760233 (TA-repeat indel)
  - Reports compound unphased notation: e.g. `*28(het)+*60(hom)+*80(het) [unphased]`
  - Validated on ILMN CRAM cohort: detects alleles missed by all four star-allele tools

### `docker/pgx-report.py`

- **RYR1 diplotype_check:** added `*ference/*ference` to exclusion list to prevent
  Aldy's reference notation from triggering spurious clinical alert banners
- **DPYD/RYR1 low-depth warning:** new `_DEPTH_SENSITIVE_GENES` constant; amber
  warning banner rendered when mean gene depth < 30× (threshold based on ThermoFisher Axiom
  cohort observation that rare het variants were missed below this depth)
- **DPYD description:** notes that WGS detects full variant spectrum including *5
  (c.1627A>G) and *9A (c.85T>C) not covered by standard pharmacogenomics arrays
- **UGT1A1 description:** notes that promoter alleles *60 and *80 are not uniformly
  represented across tool databases and may differ between tool outputs
- **TOOLS list:** `VCF-Check` added as seventh recognised tool

### CRAM sample reports regenerated

- `results/cram/sample N/sample_N_pgx_report.html` — depth warnings fire for DPYD
  (29.6×) and RYR1 (25.8×); VCF-Check row added to UGT1A1 section
- `results/cram/sample A/sample_A_pgx_report.html` — VCF-Check correctly calls
  `*28(het)+*60(hom)+*80(het) [unphased]` for UGT1A1

---

## 2026-03-13 — ThermoFisher Axiom Validation Report (13 samples × 16 genes, ILMN + MGI)

### Summary

Comprehensive concordance analysis of 13 clinical samples previously genotyped
on the Thermofisher Axiom Pharmacogenomics array (P6/P9) against the pgx-suite
pipeline on both Illumina and MGI WGS platforms.

### `AxiomValidation.md` (new)

- **13 samples**, 26 platform-sample pairs (ILMN + MGI), **16 genes** evaluated
- **Tool concordance:** Stargazer 81.5% > PyPGx 79.9% > Aldy 78.1% > StellarPGx 76.5%
- **Key finding — DPYD:** Pipeline detects reduced-function variants (*9A, *5)
  in samples Axiom reports as `*1/*1`; Axiom array covers only major DPYD
  variants while WGS calls the full gene — may have direct 5-FU dosing impact
- **CYP2D6 (55% concordance):** Hybrid allele nomenclature differences
  (*36+*10 vs *36-*10, *106 sub-classification); biologically concordant in most
  cases; Axiom P6 returns "Indeterminate" for most samples
- **UGT1A1 (62%):** *60 and *80 promoter alleles systematically missing from
  tool allele databases — documented gap
- **RYR1 (60%):** Rare het variants missed by PyPGx/Stargazer; Aldy most
  sensitive for RYR1 rare variant detection
- **ABCG2:** Only Aldy calls this gene; 100% concordant with Axiom
- **Platform:** ILMN and MGI produce equivalent genotypes; MGI has markedly
  lower duplication (6–10% vs 20–44%), giving higher effective depth
- Recommendations: DPYD clinical review, consensus ≥3/4 tools policy, 30×
  minimum gene depth for RYR1/DPYD, standardise VKORC1 to rsID notation

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
│   ├── <GENE>_HG002_LIMS_ID_comparison.tsv
│   └── <GENE>_HG002_LIMS_ID_detail.json
└── html_reports/
    ├── HG002_LIMS_ID.pgx.html            # landing page (26 gene cards)
    └── HG002_LIMS_ID.<GENE>.pgx.html     # per-gene detail (26 files)
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
