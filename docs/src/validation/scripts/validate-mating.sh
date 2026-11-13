#!/usr/bin/env bash
# Mating schemes of FisherWright against SLiM (validation 4).
#
# Usage: validate-mating.sh SLIM OUTDIR R SEED [JOBS]
#   SLIM    path to the slim executable (SLiM 5)
#   OUTDIR  directory for the SLiM output and the report
#   R       replicates per scheme and simulator (a multiple of 100)
#   SEED    base seed; SLiM process j of a scheme uses SEED * 1000 + j
#   JOBS    SLiM processes run in parallel (default 8)
# The environment variable SCHEMES restricts the run to some schemes, for
# example SCHEMES="rm avoid_cousin5".
# Run from the root of the FisherWright repository.
set -euo pipefail
slim=$1 out=$2 reps=$3 seed=$4 jobs=${5:-8}
here=$(dirname "$0")
mkdir -p "$out"
schemes=${SCHEMES:-"rm rm5 fsgamma hier mono self avoid_half avoid_cousin5 avoid_mono"}
nproc=$((reps / 100))
"$slim" -v > "$out/slim-version.txt"
for s in $schemes; do
    rm -f "$out/$s.tsv"
    seq 1 "$nproc" | xargs -P "$jobs" -I{} sh -c \
        "'$slim' -s \$(( $seed * 1000 + {} )) -d \"SCHEME='$s'\" -d N=40 -d NT=30 -d K=100 \
         '$here/validate-mating.slim' | grep -P '^$s\t' > '$out/$s.{}.part'"
    cat "$out/$s".*.part > "$out/$s.tsv"
    rm -f "$out/$s".*.part
done
{
    echo "commit $(git rev-parse --short HEAD), $(head -1 "$out/slim-version.txt"), $(uname -n), $(date -u +%Y-%m-%dT%H:%MZ)"
    julia --startup-file=no --project -t 8 "$here/validate-mating.jl" "$out" "$reps" "$seed" "$schemes"
} | tee "$out/report.txt"
