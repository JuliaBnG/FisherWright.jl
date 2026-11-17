#!/usr/bin/env bash
# Stabilizing selection in FisherWright against SLiM.
#
# Usage: validate-stabilizing.sh SLIM OUTDIR RUNS SEED [JOBS]
#   SLIM    path to the slim executable (SLiM 5)
#   OUTDIR  directory for the output of both simulators and the report
#   RUNS    runs per simulator (at most 99999)
#   SEED    base seed; SLiM run r uses SEED * 100000 + r, FisherWright
#           run r hash((SEED, r))
#   JOBS    SLiM processes, and FisherWright threads, in parallel (default 8)
# The environment variable STEPS restricts the run to some of "slim engine
# compare", for example to run the simulators on different machines and
# compare in one OUTDIR afterwards.
# Run from the root of the FisherWright repository.
set -euo pipefail
slim=$1 out=$2 runs=$3 seed=$4 jobs=${5:-8}
here=$(dirname "$0")
mkdir -p "$out"
out=$(cd "$out" && pwd)
steps=${STEPS:-"slim engine compare"}
jl="julia --startup-file=no --project"
case " $steps " in *" slim "*)
    alpha=$($jl "$here/validate-stabilizing.jl" alpha)
    "$slim" -v > "$out/slim-version.txt"
    seq 1 "$runs" | xargs -P "$jobs" -I{} sh -c \
        "'$slim' -s \$(( $seed * 100000 + {} )) -d RUN={} -d N=500 -d KC=2000 -d NCHR=20 \
         -d CL=100000000 -d MU=5e-7 -d OMEGA2=49.0 -d SIGMAE=1.0 -d BURN=10000 -d GAP=1000 -d NS=21 \
         -d \"ALPHA='$alpha'\" '$here/validate-stabilizing.slim' | grep '^STAB' > '$out/slim.{}.part'"
    cat "$out"/slim.*.part > "$out/slim.tsv"
    rm -f "$out"/slim.*.part ;;
esac
case " $steps " in *" engine "*)
    $jl -t "$jobs" "$here/validate-stabilizing.jl" engine "$out" "$runs" "$seed" ;;
esac
case " $steps " in *" compare "*)
    {
        echo "commit $(git rev-parse --short HEAD), $(head -1 "$out/slim-version.txt" 2>/dev/null || echo 'SLiM version unknown'), $(uname -n), $(date -u +%Y-%m-%dT%H:%MZ)"
        $jl "$here/validate-stabilizing.jl" compare "$out"
    } | tee "$out/report.txt" ;;
esac
