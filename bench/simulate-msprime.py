#!/usr/bin/env python3
"""
Simulate a multi-autosome population using msprime and convert SNP genotypes to a bit-packed array.
"""

import sys
import time
import argparse
import resource
import numpy as np

try:
    import msprime
except ImportError:
    print("Error: msprime not found. Run with `uv run --with msprime --with numpy bench/simulate-msprime.py`")
    sys.exit(1)


def peak_rss_mb() -> float:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return rss / (1024 * 1024)
    return rss / 1024


def main():
    parser = argparse.ArgumentParser(description="Simulate population with msprime across autosomes.")
    parser.add_argument("--ne", type=int, default=2000, help="Diploid effective population size (default: 2000)")
    parser.add_argument("--nt", type=int, default=2000, help="Generations (default: 2000)")
    parser.add_argument("--nchr", type=int, default=10, help="Number of autosomes (default: 10)")
    parser.add_argument("--chr-len", type=int, default=100_000_000, help="Length per autosome in bp (default: 100,000,000)")
    parser.add_argument("--mut-rate", type=float, default=1e-8, help="Mutation rate per bp per generation (default: 1e-8)")
    parser.add_argument("--rec-rate", type=float, default=1e-8, help="Recombination rate per bp per generation (default: 1e-8)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    parser.add_argument("--dtwf", action="store_true", default=True, help="Use Discrete-Time Wright-Fisher model")
    args = parser.parse_args()

    nhaps = 2 * args.ne
    total_len = args.nchr * args.chr_len

    print("=================================================================")
    print(" msprime Simulation")
    print("=================================================================")
    print(f"Diploid Individuals (ne): {args.ne} (haplotypes: {nhaps})")
    print(f"Generations (nt)        : {args.nt}")
    print(f"Autosomes (nchr)        : {args.nchr} autosomes × {args.chr_len / 1e6:.1f} Mb = {total_len / 1e6:.1f} Mb total")
    print(f"Mutation rate           : {args.mut_rate} / bp / generation")
    print(f"Recombination rate      : {args.rec_rate} / bp / generation (1 cM/Mb)")
    print(f"Model                   : {'DTWF (Discrete-Time Wright-Fisher)' if args.dtwf else 'Hudson Standard Coalescent'}")
    print(f"Random seed             : {args.seed}")
    print("-----------------------------------------------------------------")

    t_start = time.perf_counter()
    total_ancestry_time = 0.0
    total_mutation_time = 0.0
    total_bit_time = 0.0
    total_sites = 0
    bit_chunks = []

    model = msprime.DiscreteTimeWrightFisher(duration=args.nt) if args.dtwf else "hudson"

    for c in range(args.nchr):
        chr_seed = args.seed + c * 10007

        # 1. Simulate Ancestry for this autosome
        t0 = time.perf_counter()
        ts = msprime.sim_ancestry(
            samples=args.ne,
            population_size=args.ne,
            sequence_length=args.chr_len,
            recombination_rate=args.rec_rate,
            model=model,
            random_seed=chr_seed,
        )
        total_ancestry_time += time.perf_counter() - t0

        # 2. Simulate Mutations
        t0 = time.perf_counter()
        mts = msprime.sim_mutations(ts, rate=args.mut_rate, random_seed=chr_seed, model=msprime.BinaryMutationModel())
        total_mutation_time += time.perf_counter() - t0
        total_sites += mts.num_sites

        # 3. Genotype extraction and bit-packing
        t0 = time.perf_counter()
        gt = mts.genotype_matrix()
        bit_chunk = np.packbits(gt, axis=1)
        bit_chunks.append(bit_chunk)
        total_bit_time += time.perf_counter() - t0

    # Combine bit arrays across autosomes
    t0 = time.perf_counter()
    combined_bit_array = np.vstack(bit_chunks)
    total_bit_time += time.perf_counter() - t0

    total_time = time.perf_counter() - t_start
    packed_mb = combined_bit_array.nbytes / (1024 * 1024)

    print(f"Ancestry runtime        : {total_ancestry_time:.3f} s")
    print(f"Mutation runtime        : {total_mutation_time:.3f} s")
    print(f"BitArray packing time   : {total_bit_time:.3f} s")
    print(f"Total runtime           : {total_time:.3f} s")
    print(f"Polymorphic SNPs found  : {total_sites}")
    print(f"BitArray shape (packed) : {combined_bit_array.shape} (dtype: uint8)")
    print(f"BitArray memory size    : {packed_mb:.2f} MiB (8 bits per byte, 1 bit/genotype)")
    print(f"Process peak RSS        : {peak_rss_mb():.2f} MiB")
    print("=================================================================")

    return combined_bit_array


if __name__ == "__main__":
    main()
