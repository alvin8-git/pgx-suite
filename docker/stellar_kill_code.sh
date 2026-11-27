#!/usr/bin/env bash
# usage: stellar_kill_code.sh TRACE   (TRACE = the nested StellarPGx run's `-with-trace` file)
# Prints the kill code (see kill_codes.sh) of the nested task that made the nested run fail: a task whose
# LAST attempt FAILED (a task retried and then completed, and siblings aborted by Nextflow when it stopped, are
# ignored). Prints nothing when that task failed for an ordinary reason or there is no trace.
. "$(dirname "$(readlink -f "$0")")/kill_codes.sh"
[ -s "${1:-}" ] || exit 0
code=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
     {n=$c["name"]; if(!(n in st)) order[++k]=n; st[n]=$c["status"]; ex[n]=$c["exit"]}
     END{for(i=1;i<=k;i++) if(st[order[i]]=="FAILED"){print ex[order[i]]; exit}}' "$1")
is_kill_code "$code" && echo "$code"
exit 0
