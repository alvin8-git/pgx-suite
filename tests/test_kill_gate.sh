#!/usr/bin/env bash
# bin/nf_kill_gate.sh: 134/137/143 (kill) and 140 (time limit) fail the task with their own code;
# success and ordinary tool failures (1, 2) do not.
set -u
g="$(dirname "$0")/../bin/nf_kill_gate.sh"; fail=0
for rc in 0 1 2; do "$g" $rc 2>/dev/null || { echo "FAIL rc=$rc must pass"; fail=1; }; done
for rc in 134 137 140 143; do "$g" $rc 2>/dev/null; [ $? -eq $rc ] || { echo "FAIL rc=$rc must exit $rc"; fail=1; }; done
# retry policy: only 134/137 retry (memory); 143 and 140 fail fast (still gated above)
cfg="$(dirname "$0")/../nextflow.config"
grep -q 'task.exitStatus in \[134, 137\] && task.attempt' "$cfg" || { echo "FAIL errorStrategy must retry only 134 and 137"; fail=1; }
grep -q '143' <(grep 'errorStrategy *=' "$cfg") && { echo "FAIL errorStrategy must not retry 143"; fail=1; }
# nested StellarPGx: the kill code of the task that failed, never a stale one
sk="$(dirname "$0")/../docker/stellar_kill_code.sh"; T=$(mktemp -d -t sk.XXXXXX)
H=$'task_id\thash\tnative_id\tname\tstatus\texit'
tr() { printf '%s\n' "$H" "$@" > "$T/t"; "$sk" "$T/t"; }
[ "$(tr $'1\ta\t1\tp (1)\tFAILED\t137' $'2\tb\t2\tp (1)\tCOMPLETED\t0')" = "" ] || { echo "FAIL retried-then-completed 137 must be ignored"; fail=1; }
[ "$(tr $'1\ta\t1\tp (1)\tFAILED\t137' $'2\tb\t2\tp (1)\tFAILED\t137')" = 137 ] || { echo "FAIL last-attempt 137 must be reported"; fail=1; }
[ "$(tr $'1\ta\t1\tq (1)\tFAILED\t137' $'2\tb\t2\tp (1)\tCOMPLETED\t0' $'3\tc\t3\tq (1)\tCOMPLETED\t0' $'4\td\t4\tr (1)\tFAILED\t1')" = "" ] || { echo "FAIL ordinary failure beside a stale kill must report none"; fail=1; }
[ "$(tr $'1\ta\t1\tp (1)\tFAILED\t1' $'2\tb\t2\tq (1)\tABORTED\t143')" = "" ] || { echo "FAIL aborted sibling must be ignored"; fail=1; }
[ -z "$("$sk" "$T/none")" ] || { echo "FAIL missing trace must print nothing"; fail=1; }
rm -rf "$T"
[ $fail -eq 0 ] && echo "nf_kill_gate: PASS"; exit $fail
