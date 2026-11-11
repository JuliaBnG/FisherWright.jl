# Validation under a changing population size

`docs/src/validation/scripts/validate-demography.jl` checks `popsize`
against an exact expectation and against msprime DTWF with the same
demography. The model is that of the
[Population genetics validation](@ref), ``N = 100`` for 2,000
generations from a population without variation; each replicate is then
continued through `init` in two scenarios, and the whole final
population is the sample:

| Scenario | Trajectory after the burn-in | Sample |
|---|---|---|
| Bottleneck | 20 individuals for 20 generations, then 100 for 20 generations | 200 haplotypes |
| Expansion | 400 individuals for 100 generations | 800 haplotypes |

The two scenarios share their burn-in, so they are correlated with each
other; every check is within one scenario.

**Exact expectation of Σ2pq.** With ``U = \mu L`` new mutations per
haplotype per generation (infinite sites) the recursion of
[Effective size](@ref) becomes

```math
H_B' = (1 - a)H_B + \tfrac{a}{2}H_W + 2U, \qquad H_W' = H_B + 2U,
```

with ``a = E[1/(4N_m) + 1/(4N_f)]`` over the sex counts of the parents
of that generation, and ``\Sigma 2pq = H_W/(2N) + (1 - 1/N)H_B`` for
the ``N`` individuals of a cohort. It is run from ``H = 0`` through the
same generations as the simulation, so it assumes neither equilibrium
nor an effective size. At the end of the burn-in it gives 395.9, the
``\theta(1 - 1/n)`` of the constant-size validation.

**msprime with the same demography.** Two differences between the
models have to be mapped:

- *Effective size.* msprime's DTWF picks both parents of an individual
  from one pool of ``N'`` individuals, so a cohort of ``N`` two-sex
  parents is represented by ``N' = N_e(N)``: 19.38, 99.48 and 399.50
  for ``N`` = 20, 100 and 400. ``N'`` must be an integer, so each
  msprime chromosome draws ``\lfloor N_e \rfloor`` or
  ``\lceil N_e \rceil`` for every epoch, with probabilities that make
  the mean ``N_e``. The statistics are smooth in ``N'``, so their mean
  over replicates is the value at ``N_e`` up to terms of order
  ``(1/N')^2``, below 0.1% here. No statistic is rescaled.
- *Timing.* Generation ``g`` of FisherWright draws the parents of its
  ``N_g`` offspring from the ``N_{g-1}`` individuals of generation
  ``g - 1``. msprime draws the parents of the individuals at time
  ``t - 1`` from the size set for time ``t``; this was checked with a
  population reduced to one individual at time ``T``, where all
  lineages collapse into two genomes at time ``T``, not ``T + 1``
  (msprime 1.4.5). Cohort ``g`` is at time ``nt - g``, so an epoch
  that starts with cohort ``g_0`` ends, backwards, at time
  ``nt - g_0``, and the earlier size applies from time
  ``nt - g_0 + 1``.

Each scenario has 18 checks at |z| ≤ 3.5: Σ2pq against the recursion,
and Σ2pq, ``S``, five classes of the whole-population SFS, five classes
of the SFS of 20 random haplotypes and five bins of ``\sigma^2_d``
against msprime. There is no coalescent ``\theta/k`` reference out of
equilibrium, so the sample SFS is compared with msprime as well. With 36
checks, the chance that at least one fails by chance is about 2%.

It takes about 2 minutes with 8 threads for 200 replicates:

```bash
julia -t 8 validate-demography.jl [reps=100] [seed=2026]
```

## Results, 200 replicates, seed 2026

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-demography-seed2026.txt"), String)))
```

## Results, 200 replicates, seed 7

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-demography-seed7.txt"), String)))
```

## Reading the results of both runs

- All 36 checks pass for both seeds; the largest |z| is 2.45 (seed
  2026) and 2.34 (seed 7), in different statistics.
- **Σ2pq.** The bottleneck leaves 277.9 of the 395.9 at equilibrium,
  and the expansion raises it to 537.2; FisherWright is within 1.7 SE
  of both in both runs.
- **LD.** The bottleneck raises ``\sigma^2_d`` in every bin (0.24
  against 0.14 at 1–3 Mb at constant size) and the expansion lowers it
  (0.043); both simulators agree in all bins.
- **Resolution.** The standard error of Σ2pq is about 0.7% for the
  bottleneck and 0.3% for the expansion. The last lines of each run
  give what that resolves: ten more generations of expansion are 7.7
  SE away, but one more generation of bottleneck, or 21 instead of 20
  individuals, only 2 to 3 SE. An error of one generation in the
  bottleneck would therefore not be caught reliably here; the CI test
  "Heterozygosity through a bottleneck and an expansion" is the check
  aimed at that, and separates parents from offspring by more than 12
  SE.
- Not covered: growth that changes every generation (`growth`), fixed
  sex ratios combined with a changing size (the CI test of
  heterozygosity decay covers Σ2pq only), and sizes below 20.

The driver is
[`validate-demography.jl`](https://github.com/JuliaBnG/FisherWright.jl/blob/main/docs/src/validation/scripts/validate-demography.jl);
it reuses the statistics of `validate-popgen-theory.jl` and calls
`validate-popgen-msprime.py` with `--epochs`, both shown under
[Driver scripts](@ref).
