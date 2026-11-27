#!/usr/bin/env bash
# docker/docker-run.sh — Convenience wrapper for pgx-suite container
# Usage: ./docker/docker-run.sh [docker-run-args...] -- [container-command...]
# Example: ./docker/docker-run.sh -- bash /opt/pgx/test.sh
# Example: ./docker/docker-run.sh -v /data/ref:/pgx/ref -- bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Runs as the invoking user (docker/run_as_user.sh), like run_pgx_suite.sh; OG_DOCKER_USER='' gives root.
source "${SCRIPT_DIR}/run_as_user.sh"
pgx_user_opts || exit 2

# shellcheck disable=SC2086  # PGX_USER_OPTS is deliberately word-split
exec docker run --privileged $PGX_USER_OPTS --tmpfs /pgx:mode=0777 --rm -it \
  -v "${PROJECT_DIR}/StellarPGx:/pgx/stellarpgx" \
  -v "${PROJECT_DIR}/StellarPGx/containers:/pgx/containers" \
  -v "${PROJECT_DIR}/pypgx/pypgx-bundle:/pgx/bundle" \
  "$@" \
  pgx-suite:latest
