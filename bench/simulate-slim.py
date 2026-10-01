#!/usr/bin/env python3
"""
Run SLiM simulation of 10 autosomes, load the tree sequence, overlay mutations, and pack into a bit array.
"""

import os
import sys
import time
import shutil
import subprocess
import resource
import numpy as np

try:
    import tskit
    import msprime
except ImportError:
    print("Error: tskit or msprime not found. Run with `uv run --with tskit --with msprime --with numpy bench/simulate-slim.py`")
    sys.exit(1)


def peak_rss_mb() -> float:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return rss / (1024 * 1024)
    return rss / 1024


def main():
    slim_bin = os.environ.get("SLIM_BIN") or shutil.which("slim") or (
        os.path.abspath("SLiM/build/slim") if os.path.isfile("SLiM/build/slim") else None
    )
    slim_script = "bench/simulate-slim.slim"
    tree_file = "bench/slim-output.trees"

    print("=================================================================")
    print(" SLiM Forward Simulation (10 Autosomes × 1e8 bp = 1 Gb)")
    print("=================================================================")

    if not slim_bin:
        print("[Notice] `slim` executable was not found on PATH.")
        print("To run SLiM directly:")
        print("  1. Install SLiM: conda install -c bioconda slim (or build from source)")
        print(f"  2. Run recipe : slim {slim_script}")
        print("  3. Then re-run this script to overlay mutations and export the BitArray.")
        print("=================================================================")
        return None

    print(f"Using SLiM binary: {slim_bin}")
    print(f"Running recipe: {slim_script} ...")
    
    t0 = time.perf_counter()
    subprocess.run([slim_bin, slim_script], check=True)
    t_slim = time.perf_counter() - t0
    print(f"SLiM forward simulation finished in: {t_slim:.2f} s")

    # Load tree sequence
    t0 = time.perf_counter()
    ts = tskit.load(tree_file)
    t_load = time.perf_counter() - t0

    # Overlay neutral mutations
    t0 = time.perf_counter()
    mts = msprime.sim_mutations(ts, rate=1e-8, random_seed=42, model=msprime.BinaryMutationModel())
    t_mut = time.perf_counter() - t0

    # Extract genotype matrix and pack into bit array
    t0 = time.perf_counter()
    gt = mts.genotype_matrix()
    bit_array = np.packbits(gt, axis=1)
    t_bit = time.perf_counter() - t0

    packed_mb = bit_array.nbytes / (1024 * 1024)
    print(f"Tree load time          : {t_load:.3f} s")
    print(f"Mutation overlay time   : {t_mut:.3f} s")
    print(f"BitArray packing time   : {t_bit:.3f} s")
    print(f"Total post-process time : {t_load + t_mut + t_bit:.3f} s")
    print(f"Polymorphic SNPs found  : {mts.num_sites}")
    print(f"BitArray shape (packed) : {bit_array.shape} (dtype: uint8)")
    print(f"BitArray memory size    : {packed_mb:.2f} MiB")
    print(f"Peak RSS                : {peak_rss_mb():.2f} MiB")
    print("=================================================================")

    return bit_array


if __name__ == "__main__":
    main()
