#!/usr/bin/env bash
# usage: nf_kill_gate.sh RC   (called right after a caller step wrote its .status sentinel)
# A caller that merely failed (rc 1, 2, ...) keeps the sentinel convention: the task exits 0, the gene shows
# that caller as "failed". A caller KILLED (134 abort, 137 SIGKILL/OOM, 143 SIGTERM, 140 time limit) did not
# run to a verdict, so the task itself fails: nextflow.config retries 134/137/143 with more memory and, when
# retries are used up, stops the run (non-zero pgx exit). 140 is not retried (more memory does not help).
# A killed caller never becomes a silent 3/3 instead of 4/4.
rc="${1:-0}"
case "$rc" in
    134|137|143)
        echo "ERROR: caller killed (exit $rc, out of memory or signal); task fails so it is retried with more memory or the run stops" >&2
        exit "$rc" ;;
    140)
        echo "ERROR: caller hit the time limit (exit 140); not retried, the run stops (raise --max_time)" >&2
        exit "$rc" ;;
esac
