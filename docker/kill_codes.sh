# Sourced by the shell wrappers and by bin/nf_kill_gate.sh: the ONE list of exit codes that mean "killed", not
# "failed" (134 abort, 137 SIGKILL/OOM, 140 time limit, 143 SIGTERM). The Python side is KILL_CODES in pgx-compare.py.
is_kill_code() { case "${1:-}" in 134|137|140|143) return 0 ;; esac; return 1; }
