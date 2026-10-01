# Population genetics validation

This page records how `FisherWright.jl` has been checked against population
genetics theory and against an independent simulator, and what the checks
found. It includes the test code, the driver scripts and the saved output, so
the validation doesn't need to be repeated unless the simulation code changes
(see [When to rerun](@ref)).

The code blocks below are read from the repository files when the documentation
is built, so they always show the current code. The results are saved output
files; each starts with the date, commit, Julia and msprime versions it was
produced with.

## Download and run full validation

[Download the source archive](https://github.com/JuliaBnG/FisherWright.jl/archive/refs/heads/main.zip)
and extract it to obtain the validation scripts. Install FisherWright once,
then run the full validation from the directory containing the downloaded
scripts:

```bash
julia --startup-file=no -e 'using Pkg; Pkg.add("FisherWright")'
cd FisherWright.jl/docs/src/validation/scripts
julia -t 8 validate-popgen-theory.jl [reps=50] [ne=100] [seed=2026]
```

The full validation requires `uv` on the `PATH`; it downloads msprime and numpy
into a temporary environment. The [Power of the CI tests](@ref) script also
needs the repository's test suite and Git history, so run that script from a
source checkout as shown in [When to rerun](@ref).

## Summary

| Check | Reference | Result (v0.3.6) |
|---|---|---|
| Crossovers per meiosis | Poisson(L/M) per chromosome; 0.5 switch between chromosomes | Mean, variance/mean and switch probability match (\|z\| ≤ 1.8) |
| Segregating sites S | Watterson's θ·aₙ; msprime | S/θaₙ = 1.006 (msprime 1.010) |
| Heterozygosity Σ2pq | θ(1 − 1/n); msprime | 390.0 ± 3.4 vs. 398.0 (msprime 393.4 ± 3.6) |
| SFS, whole population | msprime DTWF | All 5 frequency classes match |
| SFS, sample of 20 haplotypes | coalescent θ/k | All 5 frequency classes match |
| LD decay σ²_d, 10 kb–3 Mb | msprime DTWF | All 5 distance bins match |

All 24 checks of the full validation pass for two independent seeds. The CI
tests pass for five seeds and fail for each of three deliberately introduced
errors (see [Power of the CI tests](@ref)).

These checks found one real bug, fixed in v0.3.5: the generation loop mutated
parents before mating, so the returned population lacked its newest
singletons and had about 8% fewer segregating sites than θ·aₙ.

## Model and statistics

All checks use the same neutral model:

- ``N = 100`` diploids (``n = 2N = 200`` haplotypes), constant size, random
  mating with two sexes;
- ``20N = 2000`` generations from a population with no variation, which is
  ``10 \times 2N`` generations, enough to reach mutation–drift equilibrium;
- 10 chromosomes of 10 Mb, ``\mu = r = 10^{-8}`` per bp per generation, free
  recombination between chromosomes;
- the whole population is the sample.

With ``\theta = 4N\mu L`` for total genome length ``L``:

| Statistic | Definition | Expectation used |
|---|---|---|
| ``S`` | number of segregating sites | ``\theta a_n``, ``a_n = \sum_{i=1}^{n-1} 1/i`` |
| ``\Sigma 2pq`` | ``\sum_\text{sites} 2p(1-p)`` with population frequencies | ``\theta (1 - 1/n)`` in the full validation; ``\theta`` in CI, where the 0.5% difference is far inside the bounds |
| ``\xi_k`` | number of sites with ``k`` derived copies | ``\theta / k`` for a sample much smaller than the population |
| ``\sigma^2_d`` | ``\sum D^2 / \sum p_1 q_1 p_2 q_2`` over SNP pairs in a distance bin, MAF ≥ 0.1 | msprime DTWF under the same model and filter |
| crossovers | count per chromosome per meiosis | Poisson with mean ``L_\text{chr}/M`` |

### Choice of references

- **Whole-population SFS.** The coalescent ``\theta/k`` assumes a sample much
  smaller than the population. With all ``2N`` haplotypes sampled, a discrete
  Wright–Fisher population has about 12% more singletons than ``\theta``, in
  FisherWright and msprime alike. So the whole-population SFS is compared with
  msprime's discrete-time Wright–Fisher model (DTWF), and ``\theta/k`` is
  tested on a random sample of 20 haplotypes.
- **LD.** Sved's ``E[r^2] = 1/(1 + 4Nc)`` overestimates ``r^2`` at short
  distances (0.94 vs. about 0.59 observed at 10–30 kb) and ignores the MAF
  filter. Ohta and Kimura's ``\sigma^2_d = (10 + \rho)/(22 + 13\rho + \rho^2)``,
  ``\rho = 4Nc``, is about 25% too low here for the same reason. Only msprime
  gives the same statistic under the same filter. Both formulas are printed
  for reference.
- **msprime settings.** `DiscreteTimeWrightFisher()` with no duration, so
  genealogies run to full coalescence; each chromosome is simulated separately
  (free recombination); only biallelic sites are used.

## CI tests

`test/runtests.jl` contains two population genetics testsets, run by
`Pkg.test()` in about 10 s:

- **Neutral Wright-Fisher theory**: 20 equilibrium replicates with a fixed seed
  feed three sub-testsets:
  - 0.97 < S/θaₙ < 1.06 and 0.90 < Σ2pq/θ < 1.10 (SE about 0.7% and 1.3%);
  - each SFS class of a 20-haplotype sample within 0.90–1.10 of θ/k (SE 2–3%);
  - σ²_d at 0.3–1 Mb and 1–3 Mb within 0.85–1.15 of the msprime reference
    values (SE 2–3%), which come from [LD reference values](@ref).
- **Crossovers per meiosis matches Poisson and assortment**: 10,000 draws from
  `cobp!` on chromosomes of 50 and 100 Mb; bounds are about 3.5 SE.

Each bound was set from replicate runs so that it passes for the correct model
with several seeds, and fails for the errors listed in the next section.

```@eval
using Markdown, FisherWright
src = read(joinpath(pkgdir(FisherWright), "test", "runtests.jl"), String)
a = first(findfirst("# BEGIN population genetics validation", src))
b = last(findfirst("# END population genetics validation", src))
Markdown.MD(Markdown.Code("julia", src[a:b]))
```

Latest `Pkg.test()` run:

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "ci-tests.txt"), String)))
```

## Power of the CI tests

A test that can't fail proves nothing. `docs/src/validation/scripts/check-test-power.jl` runs the
validation block of `test/runtests.jl` as it is, for five seeds, and then with
one known error at a time:

| Case | Expected outcome | Result |
|---|---|---|
| Correct model, seeds 2026, 1, 7, 42, 99 | all pass | all pass |
| Simulated Ne halved (theory assumes 100) | diversity, SFS and LD fail | all 9 fail: ratios ≈ 0.5; LD 1.35, 1.64 |
| Recombination rate doubled | LD and crossover counts fail | LD 0.65, 0.54; both crossover means fail |
| Pre-v0.3.5 generation loop | S fails | S/θaₙ = 0.919 |

Each error is caught by the test aimed at it, and only by tests that should be
sensitive to it. For example, the pre-fix loop only removed rare variants, so
Σ2pq, the sample SFS and LD are unaffected.

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "test-power.txt"), String)))
```

Script (a copy used for this documentation):

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("julia", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "scripts", "check-test-power.jl"), String)))
```

## Full validation against msprime

`docs/src/validation/scripts/validate-popgen-theory.jl` runs 50 FisherWright replicates and 50
msprime DTWF replicates of the model above (msprime through
`docs/src/validation/scripts/validate-popgen-msprime.py`, run with `uv`). Every comparison is a
z-score; a check fails when |z| > 3.5, and the script then exits with status 1.
With 24 checks, the chance that at least one fails by chance is about 1%.

Run it from the directory containing the downloaded scripts; it takes about
25 s with 8 threads:

```bash
julia -t 8 validate-popgen-theory.jl [reps=50] [ne=100] [seed=2026]
```

### Results, seed 2026

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-popgen-seed2026.txt"), String)))
```

### Results, seed 7

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-popgen-seed7.txt"), String)))
```

### Reading the results

- **Diversity.** S is within 1% of θ·aₙ in both runs, as for msprime.
  Σ2pq was slightly below msprime in both runs (z = −0.67 and −1.97). This is
  not significant; more replicates would settle it.
- **SFS.** The whole-population classes match msprime, including the excess
  of singletons over θ, which both simulators share. The sample SFS matches
  θ/k.
- **LD.** FisherWright is slightly below msprime in all bins with seed 2026
  and above it at 0.1–3 Mb with seed 7. The bins share simulations, so their
  z-scores are correlated; the change of sign shows no consistent bias.

### Driver scripts

`docs/src/validation/scripts/validate-popgen-theory.jl`:

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("julia", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "scripts", "validate-popgen-theory.jl"), String)))
```

`docs/src/validation/scripts/validate-popgen-msprime.py`:

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("python", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "scripts", "validate-popgen-msprime.py"), String)))
```

## LD reference values

The CI LD test compares with fixed msprime values, since CI can't run msprime.
They come from 100 replicates of 10 chromosomes:

```bash
cd docs/src/validation/scripts
uv run --with msprime --with numpy python validate-popgen-msprime.py --reps 100 --seed 1 --summary
```

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "ld-reference-msprime.txt"), String)))
```

The test uses the 0.3–1 Mb and 1–3 Mb values, 0.2869 and 0.1398. If the model
settings in the CI testset change, regenerate these values with matching
arguments.

## History

- **Diversity deficit (fixed in v0.3.5).** A comparison with θ·aₙ and msprime
  showed about 8% fewer segregating sites, while Σ2pq was unaffected. The
  generation loop mutated parents and then reproduced, so the returned
  population had gone through one round of reproduction since its last
  mutations. Expected loss: ``\theta \sum_k e^{-k}/k \approx 0.46\theta``, which
  matched the observed deficit. The loop now mates first and mutates the
  offspring.
- **Earlier SFS and LD tests replaced.** An earlier version compared shapes
  with a Pearson correlation, which doesn't change when Ne or the
  recombination rate is scaled. It used one replicate, compared LD with Sved's
  curve despite a 40–45% gap, and its bench script printed "PASSED"
  unconditionally. Over 10 seeds, its SFS test passed 6/10 for the correct
  model and 6/10 with Ne halved; its LD test passed 2/10 for the correct model
  and 5/10 with recombination doubled.

## Not covered

- The crossover test checks `cobp!` only. Transmission of crossovers through
  `recombine` is covered indirectly: doubling the recombination rate fails
  the LD test. `recombine` itself is checked against a reference
  implementation in the testset "recombine correctness".
- Only a constant-size, neutral, panmictic population is checked; the package
  doesn't model demography, selection or subdivision.
- Benchmark-scale runs (N = 2,000, 1 Gb) are checked only through their SNP
  count against msprime.

## When to rerun

The CI tests run on every `Pkg.test()`. Rerun the full validation and the
power check, and replace the result files in `docs/src/validation/`, when any
of these change:

- the generation loop in `src/fwp.jl`, `recombine`, `cobp!` or the
  recombination map, `random_mate!`, or `muts2bitarray`;
- the validation block in `test/runtests.jl` (bounds, settings or seeds);
- the model settings or statistics in
  `docs/src/validation/scripts/validate-popgen-theory.jl`.

```bash
D=docs/src/validation
julia -t 4 --project=. docs/src/validation/scripts/check-test-power.jl > $D/test-power.txt
julia -t 8 --project=. docs/src/validation/scripts/validate-popgen-theory.jl 50 100 2026 > $D/validate-popgen-seed2026.txt
julia -t 8 --project=. docs/src/validation/scripts/validate-popgen-theory.jl 50 100 7 > $D/validate-popgen-seed7.txt
```
