#!/usr/bin/env bash
# docker/docker-run.sh — Convenience wrapper for pgx-suite container
# Usage: ./docker/docker-run.sh [docker-run-args...] -- [container-command...]
# Example: ./docker/docker-run.sh -- bash /opt/pgx/test.sh
# Example: ./docker/docker-run.sh -v /data/ref:/pgx/ref -- bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

exec docker run --privileged --rm -it \
  -v "${PROJECT_DIR}/StellarPGx:/pgx/stellarpgx" \
  -v "${PROJECT_DIR}/StellarPGx/containers:/pgx/containers" \
  -v "${PROJECT_DIR}/pypgx/pypgx-bundle:/pgx/bundle" \
  "$@" \
  pgx-suite:latest
