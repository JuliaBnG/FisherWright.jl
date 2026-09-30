"""msprime DTWF reference statistics for bench/validate-popgen-theory.jl.

Simulates the same neutral diploid Wright-Fisher population as the Julia script
(whole population sampled, chromosomes simulated independently, which is free
recombination between them) and prints one line per replicate:

    S H w1..w5 s1..s5 num1..num5 den1..den5

S, H: segregating sites and Σ2pq in the whole population; w: whole-population
SFS groups; s: SFS groups in a random sample of `nsub` haplotypes; num/den:
ΣD² and Σp₁q₁p₂q₂ per LD distance bin (MAF filter as in the Julia script).

With --summary it prints instead the pooled σ²_d ± SE per LD bin, which is how
the reference values in test/runtests.jl were made:

    uv run --with msprime --with numpy python bench/validate-popgen-msprime.py --reps 100 --seed 1 --summary
"""

import argparse

import msprime
import numpy as np

WHOLE = [(1, 1), (2, 2), (3, 5), (6, 20), (21, None)]
SUB = [(1, 1), (2, 2), (3, 4), (5, 9), (10, 19)]
BINS = [(10_000, 30_000), (30_000, 100_000), (100_000, 300_000), (300_000, 1_000_000), (1_000_000, 3_000_000)]


def groups(counts, spec, n):
    out = []
    for lo, hi in spec:
        hi = n - 1 if hi is None else hi
        out.append(((counts >= lo) & (counts <= hi)).sum())
    return out


def ld_sums(G, pos, n, maf):
    p = G.sum(1) / n
    keep = (p >= maf) & (p <= 1 - maf)
    G, p, x = G[keep].astype(float), p[keep], pos[keep]
    D = G @ G.T / n - np.outer(p, p)
    V = np.outer(p * (1 - p), p * (1 - p))
    d = x[None, :] - x[:, None]
    num, den = [], []
    for lo, hi in BINS:
        m = (d >= lo) & (d < hi)
        num.append((D[m] ** 2).sum())
        den.append(V[m].sum())
    return np.array(num), np.array(den)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ne", type=int, default=100)
    ap.add_argument("--nchr", type=int, default=10)
    ap.add_argument("--chr-len", type=int, default=10_000_000)
    ap.add_argument("--mu", type=float, default=1e-8)
    ap.add_argument("--r", type=float, default=1e-8)
    ap.add_argument("--reps", type=int, default=50)
    ap.add_argument("--nsub", type=int, default=20)
    ap.add_argument("--maf", type=float, default=0.1)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--summary", action="store_true", help="print pooled LD per bin instead of replicate rows")
    a = ap.parse_args()
    n = 2 * a.ne
    rng = np.random.default_rng(a.seed)
    ratios, NUM, DEN = [], np.zeros(len(BINS)), np.zeros(len(BINS))
    for rep in range(a.reps):
        S = H = 0.0
        w = np.zeros(len(WHOLE))
        s = np.zeros(len(SUB))
        num = np.zeros(len(BINS))
        den = np.zeros(len(BINS))
        sub = rng.choice(n, a.nsub, replace=False)
        for c in range(a.nchr):
            seed = a.seed * 1_000_003 + rep * a.nchr + c + 1
            ts = msprime.sim_ancestry(
                samples=a.ne,
                population_size=a.ne,
                sequence_length=a.chr_len,
                recombination_rate=a.r,
                model=msprime.DiscreteTimeWrightFisher(),
                random_seed=seed,
            )
            ts = msprime.sim_mutations(ts, rate=a.mu, random_seed=seed)
            G = ts.genotype_matrix()
            bi = G.max(1) == 1  # biallelic sites only
            G, pos = G[bi], ts.tables.sites.position[bi]
            k = G.sum(1)
            f = k / n
            S += len(k)
            H += (2 * f * (1 - f)).sum()
            w += groups(k, WHOLE, n)
            ks = G[:, sub].sum(1)
            s += groups(ks[(ks > 0) & (ks < a.nsub)], SUB, a.nsub)
            x, y = ld_sums(G, pos, n, a.maf)
            num += x
            den += y
        if a.summary:
            ratios.append(num / den)
            NUM += num
            DEN += den
        else:
            print(" ".join(str(v) for v in [S, H, *w, *s, *num, *den]))
    if a.summary:
        se = np.std(ratios, axis=0, ddof=1) / np.sqrt(a.reps)
        print(f"msprime DTWF, ne = {a.ne}, {a.reps} reps x {a.nchr} chromosomes of {a.chr_len} bp, mu = {a.mu}, r = {a.r}, MAF >= {a.maf}")
        for (lo, hi), v, e in zip(BINS, NUM / DEN, se):
            print(f"  {lo / 1e6:g}-{hi / 1e6:g} Mb  sigma2_d = {v:.4f} +/- {e:.4f}")


if __name__ == "__main__":
    main()
