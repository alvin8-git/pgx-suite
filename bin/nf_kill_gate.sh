#!/usr/bin/env bash
# usage: nf_kill_gate.sh RC   (called right after a caller step wrote its .status sentinel)
# A caller that merely failed (rc 1, 2, ...) keeps the sentinel convention: the task exits 0, the gene shows
# that caller as "failed". A caller KILLED (kill_codes.sh: 134 abort, 137 SIGKILL/OOM, 140 time limit, 143 SIGTERM)
# did not run to a verdict, so the task itself fails and is never silent. nextflow.config retries only 134 and 137
# (more memory can help); 143 (qdel, soft limit) and 140 (time) fail the run at once, without a retry.
for f in "$(dirname "$0")/../docker/kill_codes.sh" /opt/pgx/kill_codes.sh; do [ -r "$f" ] && { . "$f"; break; }; done
rc="${1:-0}"
if is_kill_code "$rc"; then
    case "$rc" in
        134|137) echo "ERROR: caller killed (exit $rc, out of memory); task fails so it is retried with more memory or the run stops" >&2 ;;
        *)       echo "ERROR: caller killed (exit $rc, signal or time limit); not retried, the run stops" >&2 ;;
    esac
    exit "$rc"
fi
