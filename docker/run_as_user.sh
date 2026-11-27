#!/usr/bin/env bash
# docker/run_as_user.sh — the one place that builds the "run the container as the invoking user" options.
# Source of truth for run_pgx_suite.sh, tests/smoke/run.sh and docker/docker-run.sh; nextflow.config
# calls it directly (bash docker/run_as_user.sh [DIR]) when nextflow is launched without run_pgx_suite.sh.
#
#   source docker/run_as_user.sh; pgx_user_opts [DIR]   # sets PGX_USER_OPTS (string, word-split on use)
#   bash docker/run_as_user.sh [DIR]                    # prints the same string
#
# DIR: a folder the run owns, for the generated passwd/group files (default: mktemp -d under TMPDIR, /var/tmp if
# unset); the folder used is left in PGX_IDS_DIR so the caller can remove it.
# OG_DOCKER_USER: unset = this user's uid:gid; empty = root (explicit opt-out, no options);
#   "--user ..." (starts with "--") = used as-is; "uid" or "uid:gid" = "--user" is prepended; else an error.
# Apptainer inside the image needs the uid in /etc/passwd ("unknown userid" otherwise). The files mounted are
# generated from getent (so LDAP/SSSD users work and no other host account is exposed): this user, this
# group and a root entry. Only if getent fails are the host's /etc/passwd and /etc/group mounted instead.
# HOME=/tmp is the container's own throw-away /tmp.
pgx_user_opts() {
    local dir="${1:-}" user pw gr
    PGX_USER_OPTS=""
    if [[ -v OG_DOCKER_USER ]]; then
        user="$OG_DOCKER_USER"
        if [[ -z "$user" ]]; then return 0; fi          # explicit opt-out: root
        if [[ "$user" == --* ]]; then :
        elif [[ "$user" =~ ^[0-9]+(:[0-9]+)?$ ]]; then user="--user $user"
        else
            echo "OG_DOCKER_USER='$user' is not valid: use '' (root), '--user ...', 'uid' or 'uid:gid'" >&2
            return 2
        fi
    else
        local u g; u="$(id -u)"; g="$(id -g)"
        if [[ -z "$u" || -z "$g" ]]; then echo "cannot determine uid/gid (id -u / id -g returned nothing)" >&2; return 2; fi
        user="--user $u:$g"
    fi
    if [[ -z "$dir" ]]; then dir="$(mktemp -d "${TMPDIR:-/var/tmp}/pgx_ids.XXXXXX")" || return 2; fi
    mkdir -p "$dir" || return 2
    PGX_IDS_DIR="$dir"
    pw="$(getent passwd "$(id -u)" 2>/dev/null)" || pw=""
    gr="$(getent group "$(id -g)" 2>/dev/null)" || gr=""
    if [[ -n "$pw" && -n "$gr" ]]; then
        { echo "root:x:0:0:root:/root:/bin/sh"; echo "$pw"; } > "$dir/passwd"
        { echo "root:x:0:"; echo "$gr"; } > "$dir/group"
        PGX_USER_OPTS="$user -v $dir/passwd:/etc/passwd:ro -v $dir/group:/etc/group:ro -e HOME=/tmp"
    else
        echo "getent failed for uid $(id -u) / gid $(id -g): mounting the host /etc/passwd and /etc/group" >&2
        PGX_USER_OPTS="$user -v /etc/passwd:/etc/passwd:ro -v /etc/group:/etc/group:ro -e HOME=/tmp"
    fi
}
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    pgx_user_opts "$@" || exit 2
    printf '%s\n' "$PGX_USER_OPTS"
fi
