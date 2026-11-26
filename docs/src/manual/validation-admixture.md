# Validation of an admixture pulse

FisherWright has no model of subdivision. A single admixture pulse is
built from separate calls through `init`, as described in
[Merging populations (admixture)](@ref). This page checks that recipe
against two expectations that hold for any neutral pulse: the ancestry
proportion and the decay of linkage disequilibrium (LD).

## Model

- N = 500 diploids in every population, with two sexes (sex is Binomial(N, ½)
  in each generation, conditioned on both sexes) and random mating.
- 2 chromosomes of 50 Mb at 1e8 bp per Morgan, which is 0.5 M each.
  Crossovers are Poisson, so Haldane's map ``r = (1 - e^{-2d})/2`` holds
  exactly. Chromosomes recombine freely, and μ = 1e-8 per bp.
- **Ancestral population:** 10N = 5000 generations from no variation.
- **Split:** populations A and B each start from the ancestral haplotypes and
  evolve for 100 generations with mutation, in separate calls with
  `result = false`.
- **Merge:** 150 individuals drawn at random from A and 350 from B, as whole
  individuals, so α = 0.3. A second variant shuffles the 1000 merged
  haplotypes among individuals.
- **Continuation:** 20 generations at N = 500 with `mutate = false`, one call
  per generation, so the same sites are followed throughout.
- **Sites:** polymorphic in the merged population with a minor-allele
  frequency of at least 0.05, thinned evenly to 1000 per chromosome.
  All statistics are computed on the whole population.

Every replicate has its own ancestral history. The two variants share the
history and the merged individuals within a replicate, so each check is
within one variant.

## References

**Ancestry proportion.** Let ``p_A`` and ``p_B`` be the allele frequencies
of the merged individuals from A and B, and ``\delta = p_A - p_B``. The
least-squares estimate is

```math
\hat\alpha_t = \frac{\sum (p_t - p_B)\,\delta}{\sum \delta^2}.
```

It equals α at the merge exactly. For later t, ``E[\hat\alpha_t] = \alpha``,
because drift without mutation or selection leaves ``E[p_t] = p_0``.

**LD decay.** Let ``D_t`` be the LD of a pair of sites at recombination
fraction r, and ``C_0`` the covariance of the two sites between the two
haplotypes of one individual at the merge:

```math
C_0 = \frac{1}{2N}\sum_i (x_{i1} y_{i2} + x_{i2} y_{i1}) - p_x p_y .
```

A gamete is parental with probability ``1 - r``. Otherwise it joins the
alleles of the parent's two haplotypes. This gives
``E[D_1] = \lambda[(1 - r) D_0 + r C_0]``. From generation 1 on, individuals
are random unions of gametes, so ``E[D_{t+1}] = \lambda(1 - r)E[D_t]``.
Together:

```math
E[D_t] = \lambda^t (1 - r)^{t-1}\,[(1 - r) D_0 + r C_0], \qquad t \ge 1,
```

where ``\lambda = 1 - 1/(2N_e)`` and ``N_e`` comes from the two-sex
recursion of [Population genetics validation](@ref).

- **Whole individuals:** ``C_0`` is about ``\alpha(1-\alpha)\delta_x\delta_y``,
  the between-source part of ``D_0``. That part therefore does not decay in
  the first generation, whose offspring are F1.
- **Shuffled haplotypes:** ``C_0 \approx 0``, and the usual
  ``\lambda^t(1 - r)^t D_0`` holds.

For every distance bin, the check compares ``S_t = \sum D_t D_0`` with
``P_t = \sum E[D_t] D_0`` over the pairs in the bin. Since
``E[S_t \mid \text{merge}] = P_t``, the difference ``S_t - P_t`` has mean zero
over replicates. The prediction neglects terms of order ``r/N`` in the drift
factor.

Each comparison is a z-score over replicates, and a check fails when
|z| > 3.5. A run has 84 checks: 2 variants × 6 times × (1 ancestry + 6
distance bins). If they were independent, the chance that at least one fails
by chance would be about 4%. They are not independent: the times of one
replicate share its history. The resolution section gives the z-score of the
prediction that ignores ``C_0``, for the merge of whole individuals.

200 replicates take about 30 minutes on 16 threads:

```bash
julia -t 16 validate-admixture.jl [reps=40] [seed=2026]
```

## Results, 200 replicates, seed 2026

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-admixture-seed2026.txt"), String)))
```

## Results, 200 replicates, seed 7

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-admixture-seed7.txt"), String)))
```

## Reading the results

READING

## Not covered

- **Selection:** plain haplotypes carry no selected alleles, so a merge is
  neutral. Histories run with `fitness` are not merged.
- **Continued gene flow:** exchanging migrants between one-generation calls
  is possible but not checked here.
- **More than two sources, or several pulses:** these follow from the same
  recipe, but are not checked.
- **Estimators of real data:** the ancestry estimate uses the true source
  frequencies of the merged individuals. Methods such as ADMIXTURE, which
  estimate them, are a question for the data, not for the simulator.

## Script

[`validate-admixture.jl`](https://github.com/JuliaBnG/FisherWright.jl/blob/main/docs/src/validation/scripts/validate-admixture.jl)
uses the helpers of `validate-popgen-theory.jl` in the same directory.
