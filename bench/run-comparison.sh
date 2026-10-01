#!/usr/bin/env bash
# Comparison runner between FisherWright.jl, msprime, and SLiM across 10 autosomes
set -e

NE=${1:-2000}
NT=${2:-2000}
NCHR=${3:-10}
CHRLEN=${4:-100000000} # 100 Mb per autosome = 1 Gb total
MR=${5:-1.0}           # 1.0 -> 1e-8 / bp
SEED=${6:-42}
THREADS=${7:-8}

TOTAL_MB=$(( NCHR * CHRLEN / 1000000 ))

echo "=========================================================================="
echo " Benchmarking Forward & Backward Population Genetics Simulators"
echo " Parameters: Ne=$NE (nhaps=$((2*NE))), Nt=$NT, Autosomes=$NCHR × $((CHRLEN/1000000)) Mb = ${TOTAL_MB} Mb"
echo " Threads: $THREADS, Seed: $SEED"
echo "=========================================================================="
echo ""

echo ">>> [1/3] Running FisherWright.jl (Pure Julia Forward Simulator)..."
julia -t "$THREADS" --project=. bench/simulate-fisher-wright.jl "$NE" "$NT" "$NCHR" "$CHRLEN" "$MR" "$SEED"
echo ""

echo ">>> [2/3] Running msprime (Backwards DTWF / Coalescent in Python)..."
uv run --with msprime --with numpy bench/simulate-msprime.py \
    --ne "$NE" \
    --nt "$NT" \
    --nchr "$NCHR" \
    --chr-len "$CHRLEN" \
    --mut-rate 1e-8 \
    --rec-rate 1e-8 \
    --seed "$SEED" \
    --dtwf
echo ""

echo ">>> [3/3] Checking SLiM (Gold Standard Forward Simulator)..."
uv run --with tskit --with msprime --with numpy bench/simulate-slim.py
