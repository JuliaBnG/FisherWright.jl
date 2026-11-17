# Validation of selection

Selection in `fisher_wright` (`fitness = SiteFitness(loci)`) is checked
by the fixation probability of one new copy of a selected allele, which
depends on every part of the mechanism: the fitness of each genotype,
the draw of parents in proportion to fitness, and the transmission of
the selected allele through meiosis.

## Model

- ``N`` diploid individuals, discrete non-overlapping generations, two
  sexes. Sex is Binomial(``N``, ½) per generation, conditioned on both
  sexes being present, or exactly ``N/2`` males where stated.
- One selected site; a genotype has fitness ``1``, ``1 + hs`` or
  ``1 + s`` with zero, one or two copies of the derived allele.
- Soft selection within each sex: every offspring draws its sire among
  the males and its dam among the females with probability proportional
  to fitness, with replacement.
- A replicate starts from one copy on a random haplotype of a population
  without other variation, has no further mutation, and runs to loss or
  fixation of the allele.

## Reference

Kimura's diffusion formula and the exact Markov chain of an ideal
Wright–Fisher population describe a single pool of ``2N`` genes. The
two-sex model differs from that by terms of order ``1/N``:

- the binomial sex ratio lowers the variance effective size to about
  ``4N_mN_f/N``, and the probability of fixation of a favoured allele is
  proportional to ``N_e/N``;
- a rare allele competes only within the sex of its carrier, where its
  frequency is twice its frequency in the population, so its advantage
  over the mean of its competitors is reduced by about ``s/(2N)`` more
  than in a single pool.

The reference for the implementation is therefore an independent
simulation of the same model, `fixation-reference.jl`. It follows one
locus as genotype counts and shares no code with the package: sex by
coin flips, parents by rejection sampling on fitness, one Bernoulli
draw per transmitted gene. The ideal Wright–Fisher value (`wf_exact`,
the exact chain for ``2N`` genes) is printed beside it.

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "fixation-reference.txt"), String)))
```

- At ``N = 30`` and ``s = 0.2`` the two-sex model is 8% below the ideal
  Wright–Fisher value for ``h = ½``, 13% for a recessive and 6% for a
  dominant allele. A neutral allele fixes with probability
  ``1/(2N) = 0.01667``, as it must (0.01669 ± 0.00009).
- A deleterious allele (``s = -0.1``) fixes more often than in the
  ideal population, 0.00031 against 0.00021: selection is weaker in
  both directions.
- With ``2Ns = 10`` the deficit is 5.0% at ``N = 50`` and 2.4% with
  exactly ``N/2`` males, so about half of it is the binomial sex ratio.
  At ``N = 500`` both are within 2.4% of the ideal value and within
  1.1 standard errors of it: the difference vanishes as ``N`` grows.

## Engine against the reference

`validate-selection.jl` runs the same cases with the package's engine
and compares each with the naive reference by a z-score of the
difference of two proportions; a check fails when |z| > 3.5.

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-selection.txt"), String)))
```

All eight cases pass, with |z| ≤ 1.64, at 200,000 engine replicates
and 2,000,000 reference replicates per case. The standard error of the
engine is 0.5% of the value for ``s = 0.2`` and 0.7% for ``2Ns = 10``,
so the 5 to 13% between the two-sex model and the ideal population are
resolved many times over: the engine agrees with the two-sex reference
and not with the ideal Wright–Fisher column. The reference values of
this run differ from those of the first table by sampling error only;
the two runs use different seeds.

## Recurrent mutation at many sites

The fixation probability checks one new allele at one site.
`validate-selection-balance.jl` checks the full multi-site path of
`fisher_wright`: mutation at the fitness loci, selection on the product
of fitness over sites, and the extraction of selected substitutions.

**Design.** ``N = 30``, ``K`` selected sites, each on its own
chromosome of 100 bp, so that the sites recombine freely and the
neutral process is switched off. The population starts without
variation and runs for ``T = 300`` generations (``10N``); at the end,
each site is classified by its derived-allele count, with substituted
sites in the class "fixed". The reference is a naive simulation of one
site of the same model (the mutation probability per gene and
generation, ``1 - e^{-\mu}``, follows from the engine's Poisson hits).
Sites in one population share its pedigree, so the unit of replication
of the engine is the run; standard errors are taken over runs.

| Case | ``K`` | ``s`` | ``h`` | ``\mu`` | What it covers |
|---|---:|---:|---:|---:|---|
| 1 | 10 | −0.1 | 0.5 | 0.002 | mutation–selection–drift balance, additive |
| 2 | 10 | −0.1 | 0.1 | 0.001 | nearly recessive deleterious alleles |
| 3 | 50 | 0 | | 0.002 | neutral sites: drift and neutral substitutions |
| 4 | 20 | 0.05 | 0.5 | 0.0005 | beneficial alleles: sweeps and substitutions |
| 5 | 50 | −0.1 | 0.5 | 0.002 | as 1 with more sites: interference (not checked) |

The reference omits the other ``K - 1`` sites. Selection at unlinked
sites adds variance in offspring number and lowers the effective size
at each site (Robertson 1961), by a factor of about ``1/(1 + 4C^2)``,
``C^2`` the variance of relative fitness (Santiago and Caballero 1995).
The checked cases keep ``C^2`` at 0.002 or less, where the effect is
below 1% and below the resolution. Case 5 is case 1 with 50 sites,
``C^2 = 0.010``, and is reported but not checked; it is compared with
the reference at ``N`` and at ``N - 1``.

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-selection-balance.txt"), String)))
```

- All 28 checked comparisons pass, with |z| ≤ 1.57, at about 200,000
  engine sites and 400,000 reference replicates per case. The
  standard error of the mean derived frequency is 0.5% of its value
  for the deleterious cases.
- The neutral case also checks neutral substitutions: 32.6% of sites
  are fixed after 300 generations in both. In the beneficial case 24.0%
  are fixed, and the sites still segregating match class by class.
- Case 5 shows the interference between sites. With 50 sites, the
  mean derived frequency is 3.9% above the single-site reference
  (z = 5.9) and deleterious fixations are 12.5% more frequent. The
  predicted loss of effective size is ``1/(1 + 4 \cdot 0.0103) = 0.960``,
  and the engine agrees with the reference at ``N = 29`` (mean
  frequency 0.0677 against 0.0678). This is a property of the model,
  not of the implementation: in runs with many selected sites, the
  effective size at each site is below that of a neutral population.
- The deterministic values are not reached at this population size
  (``4Nh|s| = 6`` and ``1.2``): drift raises the frequency of deleterious
  alleles above ``\mu/(h|s|)`` (0.065 against 0.040 for case 1) and
  lowers the mean fitness below ``e^{-2K\mu}`` (0.936 against 0.961).

## Weak selection at large N: Poisson random field

The checks above compare the engine with exact simulations of the same
model at ``N = 30``. `validate-prf.jl` compares it with diffusion
theory at large ``N`` and weak selection: the stationary site frequency
spectrum of the Poisson random field (Sawyer and Hartl 1992).

**Design.** ``N = 1{,}000`` (binomial sex,
``N_e = 1/a + ½ = 999.5``), run for ``T = 20N`` generations from a
population without variation, so that ``e^{-10}`` of the slowest
neutral mode is left. There are ``K = 500{,}000`` selected sites,
evenly spaced on 20 chromosomes of 1 Morgan, in five interleaved
classes of 100,000 sites with ``\alpha = 2N_e s \in \{-10, -2, 0, 2,
10\}`` and ``h = ½``. Mutation is ``4N\mu = 10^{-4}`` per site, and
neutral mutation is off.

**Infinite sites.** The fitness loci are finite: a fixed site stays
fixed, and a hit on a carried site has no effect. The low per-site
rate keeps recurrent hits at about ``10^{-4}`` of the spectrum. At most
0.44% of the sites fix (``\alpha = 10``), so the influx is that of the
sites not fixed at the end, ``K_a``, to well under 0.1% over the
lifetime of the variants seen. About 400 sites segregate at a time,
so interference between them is negligible.

**Statistic and reference.** For each run and class, the expected
spectrum ``\xi_i`` of a sample of ``n = 20`` genes, drawn without
replacement from the ``2N``, is computed exactly (hypergeometric) from
the derived-allele counts. It is divided by ``K_a`` and by
``\theta = 4N_e\mu``, and summed over bins of ``i``. The reference is

```math
\frac{\xi_i}{\theta} = \int_0^1 \binom{n}{i} x^i (1 - x)^{n-i}
\frac{1 - e^{-\alpha(1-x)}}{(1 - e^{-\alpha})\,x(1 - x)}\,dx,
```

with ``1/i`` for ``\alpha = 0``. This is the diffusion limit for genic
selection, ``\Delta x = (s/2)x(1 - x)``, with drift variance
``x(1 - x)/(2N_e)``. Its error at ``N = 1{,}000`` is of order ``1/N``.
The unit of replication is the run: 1,000 runs per seed, standard
errors over runs, and failure at ``|z| > 3.5``.

**Results, seed 2026** (64 threads, 2 hours):

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-prf-seed2026.txt"), String)))
```

**Results, seed 7** (32 threads, 3 hours):

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-prf-seed7.txt"), String)))
```

- All 60 bins pass, with ``|z| \le 1.72``. The total
  ``\sum_i \xi_i/\theta`` is within 0.7% of the PRF in every class and seed, at a standard
  error of 0.4% to 0.7%. Pooled over the two seeds and five classes,
  the totals are 0.14% below the PRF, at ``z \approx -1.0``. This is
  within the ``1/N`` error expected of the diffusion.
- The classes are far apart: the total is 1.19, 2.69, 3.55, 4.40 and
  5.91 for ``\alpha`` from −10 to 10. Near ``\alpha = 2``, the total
  changes by about 0.43 per unit of ``\alpha``. So an error of 10% in
  the effective strength of selection would shift it by 2%, about four
  standard errors in one seed.
- The shape is resolved as well. Beneficial alleles shift the spectrum
  to high frequencies (bin 10–19: 2.46 against 0.72 neutral), and
  deleterious alleles to singletons (0.70 of a total of 1.19 at
  ``\alpha = -10``). Every bin matches within its standard error; the
  largest relative departure, 6% (``z = 0.65``), is in the sparse bin
  10–19 at ``\alpha = -10``.
- Selection within each sex at ``N_e = 999.5``, with
  ``s = \alpha/(2N_e)``, therefore behaves as genic selection in the diffusion limit. The
  deviations of order ``1/N`` from the ideal Wright–Fisher model that
  the fixation check shows at small ``N`` are below the resolution
  here.

## Selection on a trait

`validate-trait.jl` checks `TraitFitness` against quantitative-genetic
theory. The mechanism itself, that one generation of selection moves
the mean genotypic value by exactly the fitness-weighted mean of the
parents within each sex, is tested exactly in CI ("Selection on a
trait").

**Neutral genetic variance.** ``K = 20{,}000`` additive loci
(``h = ½``, so the allele substitution effect is ``\alpha = a/2``) with
effects ``a \sim N(0, 0.2^2)`` drawn anew for every replicate, mutation
``\mu = 10^{-5}`` per site (``U = K\mu = 0.2`` per haplotype and
generation), constant fitness, a population without variation at the
start. For allele indicators, the sample variance of the dosage among
the ``N`` individuals has expectation ``2H_B - H_W``; with independent
effects of zero mean the loci add, so

```math
E[\operatorname{var} G] = E[\alpha^2]\,(2H_B - H_W),
```

with ``H_B`` and ``H_W`` from the two-sex recursion with ``2U`` new
alleles per pair of genes (infinite sites). The expectation is exact at
every generation, not only at equilibrium, where it is
``2N_eV_m`` with ``N_e = 1/a`` and ``V_m = 2UE[\alpha^2]`` (Lynch and
Hill 1986). With ``\sigma_e > 0``, also
``E[\operatorname{var} P] = E[\operatorname{var} G] + \sigma_e^2``.

**Directional response, for comparison.** With fitness ``e^{\beta P}``
and Gaussian ``G`` and ``P``, the fitness-weighted mean of ``G`` is
shifted by ``\beta V_G`` and its variance is unchanged (no Bulmer
effect), so the response per generation is about ``\beta V_G``. In a
finite population the weighted mean over ``n`` parents of a sex is a
ratio estimator, with a first-order bias of
``-\beta V_G e^{\beta^2 V_P}/n``.

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-trait.txt"), String)))
```

- All 22 checks of the neutral variance pass, |z| ≤ 1.93, at 1,000
  replicates per case, from the second generation (one round of
  mutation) to near equilibrium (``0.385 \pm 0.005`` against
  ``0.389`` after 500 generations at ``N = 50``, where ``2N_eV_m`` is
  ``0.392``). ``\operatorname{var} P - \operatorname{var} G`` is
  ``0.998 \pm 0.008`` for ``\sigma_e^2 = 1``.
- The cumulative response over 30 generations of directional selection
  is ``0.923 \pm 0.008`` of ``\beta\sum V_G``. The finite-``N`` bias
  of the weighted mean explains most of the difference (prediction
  0.943); the rest, ``2.0 \pm 0.8``%, is within what finite loci and
  drift within a generation would give, and is not tested.

## Stabilizing selection: matched SLiM model

`validate-stabilizing.jl` checks the genetic variance maintained by
mutation, Gaussian stabilizing selection and drift (`TraitFitness` with
`Stabilizing`). The theory for this balance is approximate, so the test
reference is a matched model in SLiM (`validate-stabilizing.slim`), an
independent implementation with SLiM's own meiosis and mutation; the
diffusion theory is reported alongside for comparison.

**Model.** ``N = 500`` diploids. In every generation the roles of sire
and dam are assigned anew (binomial sex, conditioned on both;
``N_e = 499.5``), and every offspring draws its sire among the males
and its dam among the females in proportion to fitness, with
replacement. ``P = G + e``, ``e \sim N(0, 1)`` drawn anew, fitness
``\exp(-P^2/2\omega^2)`` with ``\omega^2 = 49`` and the optimum at 0.
There are ``K = 8{,}000`` additive trait loci (``h = ½``), evenly
spaced on 20 chromosomes of 1 Morgan, in four interleaved classes of
2,000 with substitution effect ``|\alpha|`` fixed per class and a sign
drawn anew for each run. The sizes give
``B = 2N_e\alpha^2/V_s \approx 1, 4, 16, 64`` at
``V_s = \omega^2 + V_P \approx 51``. Mutation to the derived allele
is at ``4N\mu = 10^{-3}`` per site; the sites are biallelic, with no
back mutation, and a fixed site stays fixed (SLiM: stacking policy
"f", no conversion to substitutions). Neutral mutation is off. SLiM
runs the model as nonWF with explicit reproduction, so that the sexes,
the phenotypes and the weighted draw are coded exactly as above.

**Protocol and statistics.** From a population without variation,
``20N`` generations of burn-in, then 21 samples ``2N`` apart (30,000
generations per run). Per sample: the genic variance
``V_a = \sum 2pq\alpha^2`` of each class, divided by the number of
sites not fixed, ``K_a``, and by ``4\mu V_s``; the fixed sites per
class; the mean genotypic value; and the means of ``V_P``, ``V_G`` and
the variance of log fitness over the preceding window. The samples of
a run are averaged and standard errors taken over 2,000 runs per
simulator and seed. A comparison of the two simulators fails at
``|z| > 3.5``; the mean ``G`` is also compared with 0, its exact
expectation under random signs.

**Theory, for comparison.** Forward mutation at biallelic loci makes
every allele unique: this is the house-of-cards family (Turelli 1984;
Bürger 2000), not the Gaussian regime (Kimura 1965; Lande 1975), which
needs a continuum of alleles per locus. With the mean at the optimum a
site is underdominant, ``\Delta x = -(\alpha^2/2V_s)x(1-x)(1-2x)``
(Robertson 1956). The single-locus diffusion gives the expected genic
variance per site, in units of ``4\mu V_s``,

```math
B \int_0^1 \frac{S(1) - S(x)}{S(1)\,\psi(x)}\,dx,
\qquad \psi(x) = e^{Bx(1-x)},\quad S(x) = \int_0^x \psi,
```

which is ``B/2`` (neutral) for ``B \to 0`` and 1 (Turelli's house of
cards) for ``B \to \infty``, and above 1 in between. Bürger's
stochastic house-of-cards formula, ``(B/2)/(1 + B/2)`` in these units,
is printed too. The diffusion leaves out the linkage disequilibrium
that selection builds: carriers of an allele of effect ``\alpha`` sit
on backgrounds shifted by about ``-\alpha k/(1 + k)``,
``k = V_G/V_s``, which weakens selection on the allele to about
``B(1 - k)^2``. The table prints the diffusion with this correction
("LD-corr."). Neither includes the drift of the mean about the optimum.

**Results, seed 2026** (SLiM 64 processes, 40 min; FisherWright 32
threads, 1.6 hours):

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-stabilizing-seed2026.txt"), String)))
```

**Results, seed 7:**

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-stabilizing-seed7.txt"), String)))
```

- FisherWright and SLiM agree on all 26 checks: |z| ≤ 1.74. The
  genic variance per class matches within 0.8% at a standard error of
  about 0.45% per simulator. ``V_P``, ``V_G`` and the variance of log
  fitness match within 0.13%. Both means of ``G`` are consistent with
  0. (No site of the ``B = 64`` class fixes, so its fixed-site row is
  empty.)
- The classes are far apart: 0.43, 1.10, 1.27 and 1.09 for
  ``B = 1, 4, 16, 64``. An error in the strength of selection would
  show at ``B = 1`` and 4, where the variance changes steeply with
  ``B``.
- ``V_G = 0.776`` is the summed genic variance, 0.790, less the
  Bulmer disequilibrium, about ``V_a^2/V_s = 0.012``.
- Against theory, both simulators exceed the plain diffusion by up to
  6% (``B = 64``), as expected of a model with ``k = 0.015``. The
  LD-corrected diffusion leaves 0.3–2.7% at ``B \ge 4``, and nothing
  at ``B = 1``. Bürger's formula is 11–39% low; it is not meant for
  biallelic, underdominant loci. These departures belong to the model,
  since SLiM shows them too; they are not tested.

## CI test

The testset "Selection under recurrent mutation" runs case 1 of the
balance check with 1,000 runs against the reference values (about 10 s
on one thread); a dominance of ``h = 0.1`` instead of ½ is 5.3 standard
errors away.

The testset "Fixation probability under selection" runs four of the
fixation (``N = 30``; ``s = 0.2`` with ``h`` = 0, ½ and 1, and ``s = 0``)
with 4,000 replicates each and compares them with the reference values
of the first table. The three values for ``s = 0.2`` are far apart, so
an error in the fitness of the heterozygote would be caught. Other CI
tests check the parts separately: the fitness of an individual against
the product over a dense genotype, the joint meiosis of neutral and
selected alleles against the recombinant of their union, the mutation
rate at the fitness loci, and the weighted draw of parents.

## Not covered

- Stabilizing selection in the Gaussian regime (a continuum of alleles
  per locus): the loci are biallelic, so only the house-of-cards
  regime is compared, against SLiM; the diffusion theory is not
  tested.
- A distribution of fitness effects rather than a few fixed values of
  ``s``, and dominance other than ½ at large ``N``.
- Strong interference between selected sites (Hill–Robertson effects).
  The effect of selected sites on linked neutral diversity is
  compared with SLiM on [Validation of linked selection](@ref).

## Scripts

[`fixation-reference.jl`](https://github.com/JuliaBnG/FisherWright.jl/blob/main/docs/src/validation/scripts/fixation-reference.jl),
[`validate-selection.jl`](https://github.com/JuliaBnG/FisherWright.jl/blob/main/docs/src/validation/scripts/validate-selection.jl)
and
[`validate-selection-balance.jl`](https://github.com/JuliaBnG/FisherWright.jl/blob/main/docs/src/validation/scripts/validate-selection-balance.jl),
[`validate-prf.jl`](https://github.com/JuliaBnG/FisherWright.jl/blob/main/docs/src/validation/scripts/validate-prf.jl),
[`validate-trait.jl`](https://github.com/JuliaBnG/FisherWright.jl/blob/main/docs/src/validation/scripts/validate-trait.jl)
and the stabilizing-selection driver
[`validate-stabilizing.sh`](https://github.com/JuliaBnG/FisherWright.jl/blob/main/docs/src/validation/scripts/validate-stabilizing.sh)
with
[`validate-stabilizing.jl`](https://github.com/JuliaBnG/FisherWright.jl/blob/main/docs/src/validation/scripts/validate-stabilizing.jl)
and
[`validate-stabilizing.slim`](https://github.com/JuliaBnG/FisherWright.jl/blob/main/docs/src/validation/scripts/validate-stabilizing.slim).
Rerun them, and replace the result files and the reference values in
`test/runtests.jl`, when the draw of parents, `SiteFitness`,
`TraitFitness`, the
mutation or the meiosis of the selected stream changes:

```bash
D=docs/src/validation
julia -t 16 $D/scripts/fixation-reference.jl 2000000 2026 > $D/fixation-reference.txt
julia -t 8 --project=. $D/scripts/validate-selection.jl 200000 2000000 2026 > $D/validate-selection.txt
julia -t 16 --project=. $D/scripts/validate-selection-balance.jl 4000 400000 2026 > $D/validate-selection-balance.txt
julia -t 8 --project=. $D/scripts/validate-trait.jl 1000 2026 > $D/validate-trait.txt
# about 2 hours on 64 threads
julia -t 64 --project=. $D/scripts/validate-prf.jl 1000 2026 > $D/validate-prf-seed2026.txt
# SLiM 5 and FisherWright, 2,000 runs each: about 40 min on 64 and
# 1.6 hours on 32 cores; STEPS="slim", "engine" or "compare" splits
# the work between machines
$D/scripts/validate-stabilizing.sh path/to/slim stab-2026 2000 2026 64
cp stab-2026/report.txt $D/validate-stabilizing-seed2026.txt
```
