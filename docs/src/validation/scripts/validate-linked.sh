#!/usr/bin/env bash
# Linked selection in FisherWright against SLiM (validation 7).
#
# Usage: validate-linked.sh SLIM OUTDIR R SEED [JOBS]
#   SLIM    path to the slim executable (SLiM 5)
#   OUTDIR  directory for the site tables, the SLiM output and the report
#   R       replicates per scenario and simulator (at most 9999)
#   SEED    base seed; SLiM replicate r of scenario i uses
#           SEED * 100000 + i * 10000 + r
#   JOBS    SLiM processes run in parallel (default 8)
# The environment variable SCENARIOS restricts the run to some scenarios,
# for example SCENARIOS="bgs"; ben and benfree run only when listed.
# Run from the root of the FisherWright repository.
set -euo pipefail
slim=$1 out=$2 reps=$3 seed=$4 jobs=${5:-8}
here=$(dirname "$0")
mkdir -p "$out"
out=$(cd "$out" && pwd)
scenarios=${SCENARIOS:-"neutral bgs sweep"}
jl="julia --startup-file=no --project"
$jl "$here/validate-linked.jl" setup "$out"
"$slim" -v > "$out/slim-version.txt"
i=0
for s in neutral bgs sweep ben benfree; do
    i=$((i + 1))
    case " $scenarios " in *" $s "*) ;; *) continue ;; esac
    seq 1 "$reps" | xargs -P "$jobs" -I{} sh -c \
        "'$slim' -s \$(( $seed * 100000 + $i * 10000 + {} )) -d \"SCENARIO='$s'\" \
         -d N=200 -d NT=1600 -d \"SITES='$out/sites-$s.tsv'\" -d \"BINS='$out/bins.tsv'\" \
         -d MUN=1e-6 -d MUS=2.5e-5 -d RBP=1e-7 -d CL=1000000 -d NCHR=4 \
         '$here/validate-linked.slim' | grep '^LINKED' > '$out/$s.{}.part'"
    cat "$out/$s".*.part > "$out/$s.tsv"
    rm -f "$out/$s".*.part
done
{
    echo "commit $(git rev-parse --short HEAD), $(head -1 "$out/slim-version.txt"), $(uname -n), $(date -u +%Y-%m-%dT%H:%MZ)"
    $jl -t 16 "$here/validate-linked.jl" compare "$out" "$reps" "$seed" "$scenarios"
} | tee "$out/report.txt"
