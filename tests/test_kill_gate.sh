#!/usr/bin/env bash
# bin/nf_kill_gate.sh: 134/137/143 (kill) and 140 (time limit) fail the task with their own code;
# success and ordinary tool failures (1, 2) do not.
set -u
g="$(dirname "$0")/../bin/nf_kill_gate.sh"; fail=0
for rc in 0 1 2; do "$g" $rc 2>/dev/null || { echo "FAIL rc=$rc must pass"; fail=1; }; done
for rc in 134 137 140 143; do "$g" $rc 2>/dev/null; [ $? -eq $rc ] || { echo "FAIL rc=$rc must exit $rc"; fail=1; }; done
[ $fail -eq 0 ] && echo "nf_kill_gate: PASS"; exit $fail
