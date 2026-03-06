# PGx Suite Docker Container

Single container for four pharmacogenomics star allele callers.
All tools configured for **GRCh38**.

## License Notice

Stargazer (UW) and Aldy (IURTC) are **non-commercial academic use only**.
This image must NOT be pushed to any public registry.

## Prerequisites

- Docker with `--privileged` support (local workstation)
- `StellarPGx/containers/stellarpgx-dev.sif` present (already in repo)
- ~3 GB disk for the built image
- GRCh38 reference FASTA (user-supplied, for StellarPGx full runs)

## Build

```bash
cd /data/alvin/PGxCallers
docker build -t pgx-suite:latest .
```

First build takes ~15-20 minutes. Subsequent builds use layer cache.

## Smoke Tests

```bash
# PyPGx, Aldy, Stargazer (fast, no volumes needed):
docker run --rm pgx-suite:latest bash /opt/pgx/test.sh

# Full suite including StellarPGx (requires volumes):
./docker/docker-run.sh -- bash /opt/pgx/test.sh
```

## Tools Installed

| Tool | Version | Invocation |
|------|---------|-----------|
| PyPGx | 0.26.0 | `pypgx` |
| Stargazer | 2.0.3 | `stargazer` |
| Aldy | 4.8.3 | `aldy` |
| StellarPGx | 1.2.7 | `nextflow run /pgx/stellarpgx/main.nf` |

## Volume Mounts (at runtime)

| Host Path | Container Path | Purpose |
|-----------|---------------|---------|
| `./StellarPGx` | `/pgx/stellarpgx` | StellarPGx pipeline scripts, DB, resources |
| `./StellarPGx/containers` | `/pgx/containers` | `stellarpgx-dev.sif` Singularity image |
| `./pypgx/pypgx-bundle` | `/pgx/bundle` | 1KGP VCFs + CNV data for PyPGx phasing |
| `/path/to/ref` | `/pgx/ref` | GRCh38 reference FASTA + .fai index |
| `/path/to/data` | `/pgx/data` | Input BAM/CRAM files |
| `/path/to/results` | `/pgx/results` | Output directory |

## Why --privileged?

StellarPGx uses Nextflow + Apptainer (Singularity) to run its pipeline steps
inside `stellarpgx-dev.sif`. Apptainer requires kernel namespace support
(`SYS_ADMIN` capability) to unpack SIF overlay filesystems.

## Phase 2 (coming — awaiting test BAM)

`pgx-run.sh <GENE> <BAM>` will run all supported tools for a given gene
and produce a side-by-side star allele comparison table.

## Troubleshooting

| Symptom | Cause | Fix |
|---------|-------|-----|
| `FATAL: could not open image` | SIF not mounted | Add `-v $(pwd)/StellarPGx/containers:/pgx/containers` |
| `FATAL: kernel too old` | Need privilege | Run with `--privileged` |
| Nextflow hangs at startup | JAR download | Ensure internet access during first run |
| `beagle.jar` error | Java not in PATH | Verify `java -version` in container |
