# Contributing

Thank you for helping. PGx Suite is a research tool; changes must not add clinical claims.

## Before you open a pull request

1. Open an issue first for anything larger than a typo, so the approach can be agreed.
2. Gene support is defined in `docker/genes.tsv`; see `docs/howto-add-a-gene.md`.
3. Run the unit tests (no Docker or data needed):

   ```bash
   for t in docker/test_*.py; do python3 "$t" || echo "FAIL $t"; done
   python3 docker/pgx_altcheck.py --selftest
   ```

4. If you changed the pipeline, build the image and run the smoke test:

   ```bash
   docker build -t pgx-suite:latest .
   tests/smoke/run.sh
   ```

   If a change alters an expected value in `tests/smoke/expected.tsv`, explain why in the pull request.
5. Add a dated entry to `docs/CHANGES.md`.

## Data rules

- Only public GIAB samples (HG001 to HG007) may appear in code, tests, docs, screenshots or example output.
- Never commit sample data other than the smoke-test slice, patient or family identifiers, host names,
  IP addresses, user names, internal file paths or credentials.
- Do not push the Docker image to a public registry: it contains Stargazer, Aldy and Cyrius, whose
  licences forbid redistribution or commercial use.
