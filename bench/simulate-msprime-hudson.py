# msprime Hudson coalescent to full coalescence, same layout as bench/simulate-msprime.py
# (per-chromosome ancestry, BinaryMutationModel, seed + c*10007). Needed because the
# --dtwf flag in bench/simulate-msprime.py (v0.3.8) cannot be switched off.
import msprime, time, sys
ne, nchr, L, rate, seed = 2000, 10, 100_000_000, 1e-8, 42
print(f"msprime {msprime.__version__}; Hudson; ne={ne}, {nchr} x {L/1e6:.0f} Mb, mu=r={rate}, seed={seed}")
ta = tm = 0.0; S = 0
for c in range(nchr):
    s = seed + c * 10007
    t0 = time.perf_counter()
    ts = msprime.sim_ancestry(samples=ne, population_size=ne, sequence_length=L,
                              recombination_rate=rate, model="hudson", random_seed=s)
    ta += time.perf_counter() - t0; t0 = time.perf_counter()
    mts = msprime.sim_mutations(ts, rate=rate, random_seed=s, model=msprime.BinaryMutationModel())
    tm += time.perf_counter() - t0
    S += mts.num_sites
    assert all(t.num_roots == 1 for t in ts.trees()), "not fully coalesced"
print(f"Ancestry runtime        : {ta:.3f} s")
print(f"Mutation runtime        : {tm:.3f} s")
print(f"Total runtime           : {ta+tm:.3f} s")
print(f"Sites (binary model)    : {S}")
