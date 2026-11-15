# Validation of linked selection

`docs/src/validation/scripts/validate-linked.sh` compares FisherWright
with SLiM 5 (`validate-linked.slim`) where selection at many linked
sites matters: background selection, selective sweeps and the
interference between them, which the single-site checks of
[Validation of selection](@ref) cannot reach and for which no exact
theory exists. The reference is an independent implementation of the
same model in SLiM.

**Model.** ``N = 200`` diploid individuals, two sexes with exactly 100
males in every generation (SLiM's WF sexual model;
`RandomMating(nmale = 100)`), soft selection within each sex, 1,600
generations (``8N``) from a monomorphic start. Four chromosomes of
1 Mb and 0.1 Morgan (``M = 10^7`` bp; crossovers Poisson along a
chromosome, independent assortment between chromosomes). Neutral
mutations at ``10^{-6}`` per bp per haplotype per generation. 2,000
selected sites in 20 genes of 10 kb (100 sites each, positions drawn
once), five genes per chromosome with centres 200 kb apart; new
derived alleles at ``2.5 \times 10^{-5}`` per site per haplotype, so
``U = 0.1`` per diploid genome. Fitness is multiplicative over sites,
``1``, ``1 + hs``, ``1 + s``.

| Scenario | Selected sites |
|---|---|
| `neutral` | ``s = 0`` at every site: the machinery of selected sites without selection |
| `bgs` | ``s \sim -\text{Gamma}(0.5, 0.04)`` per site (mean −0.02, median −0.009; ``2Ns`` below −20 at a tenth of the sites), ``h = 0.25`` |
| `sweep` | as `bgs`, with 40 sites chosen at random beneficial, ``s = 0.05``, ``h = 0.5`` |
| `ben` | the 40 beneficial sites of `sweep` only, in place |
| `benfree` | the same 40 sites, each on a chromosome of its own (1 bp, after the 4 Mb), unlinked from the rest |

`ben` and `benfree` split the sweep scenario into hitchhiking and
selection without linkage (see
[Reading the results](@ref linked-reading)); they run only when listed
in `SCENARIOS`. Neutral statistics are always taken over the 4 Mb,
with the distance bins of the 2,000 sites.

**Matching the site model.** FisherWright fixes ``s`` per site, and a
site is biallelic: a hit on a haplotype that carries the derived
allele, or on a fixed site, changes nothing. SLiM draws a new
mutation object at every hit. The SLiM model makes the two the same: a
`mutation()` callback gives a new allele the coefficient of its site,
reuses the allele object of the site when one is already in the
population, and drops a hit on a haplosome that carries it; selected
alleles are not converted to substitutions, so a fixed site keeps
blocking hits. The report gives how often each case occurs (about 9%
of hits fall on a segregating site and 1–2% on a carried or fixed
one); without the reuse, those hits would make second alleles of the
same site with a fitness of ``(1 + hs)^2`` instead of ``1 + s`` in a
carrier of both.

Neutral sites are biallelic in FisherWright too, and SLiM is matched
by its stacking policy and its statistics: a hit on a haplosome that
carries an allele at the position is dropped, neutral alleles are not
converted to substitutions, and the frequency of a neutral site is the
sum over the allele objects at its position. The difference matters
at this resolution. When a non-carrier is hit at a segregating site,
FisherWright makes the two genes identical, where SLiM without the
pooling counts a new, distinct allele, as under infinite sites; the
biallelic ``\pi_0`` below is 0.2% lower than the infinite-sites one,
and under sweeps, which carry many neutral alleles to high frequency,
the effect is larger still. An earlier run with SLiM's per-object
accounting put FisherWright 0.65% below SLiM in ``\pi`` under `sweep`
(``z \approx -3`` in each of two seeds); with the pooling, what is
left is 0.23% over four seeds (see
[Reading the results](@ref linked-reading)). The one remaining
difference between the models, coefficients kept in single precision
in SLiM, is negligible.

**Statistics.** At generation 1,600, over all ``2N`` haplotypes:
``\pi`` (``\Sigma 2pq`` per bp, with replacement) and the number of
segregating neutral sites ``S``, by distance to the nearest selected
site on the same chromosome (below 1 kb, which holds the genes, 1–5,
5–20, 20–50 and at least 50 kb); ``\pi`` over the genome and the
effective size ``N\pi/\pi_0`` it implies; neutral substitutions; the
selected sites that segregate and that are fixed (deleterious,
beneficial), the sum of derived frequencies, and the mean fitness
and mean log fitness of the population, fixed alleles included (the
load is ``1 - \bar w`` under `bgs`). Here ``\pi_0`` is the neutral
expectation per bp at generation 1,600 for a biallelic site (ancestral
to derived at rate ``\mu``, no back mutation), from the drift
recursion of [Validation of the mating schemes](@ref) carried over to
pairs of genes. For a pair of genes in different individuals (B) or
in one individual (W), let ``d`` be the probability that exactly one
is derived and ``b`` that both are, and ``q`` the probability that a
gene is derived. With ``a = 1/N``, drift gives

```math
d_B^* = (1 - a)d_B + \tfrac{a}{2}d_W, \qquad
b_B^* = (1 - a)b_B + \tfrac{a}{2}(q + b_W), \qquad
(d_W^*, b_W^*) = (d_B, b_B),
```

(a coalesced pair is both derived with probability ``q``), and then
mutation, for either pair,

```math
d' = d(1 - \mu) + 2\mu(1 - \mu)(1 - d - b), \qquad
b' = b + \mu d + \mu^2(1 - d - b), \qquad q' = q + \mu(1 - q).
```

From a monomorphic start, ``\pi_0 = d_W/(2N) + (1 - 1/N)d_B = 7.837
\times 10^{-4}`` at generation 1,600 (the infinite-sites recursion
gives ``7.853 \times 10^{-4}``, and ``4N\mu = 8 \times 10^{-4}`` at
equilibrium). The `neutral` scenario is also compared with it. Each
simulator runs 2,000 replicates per scenario; the table gives
mean ± SE and the z-score of the difference.

A seed takes about 3.5 hours, with 32 SLiM processes and FisherWright
on 16 threads:

```bash
docs/src/validation/scripts/validate-linked.sh path/to/slim outdir 2000 2026 32
```

To run some scenarios only, list them in `SCENARIOS`:

```bash
SCENARIOS="neutral sweep" docs/src/validation/scripts/validate-linked.sh path/to/slim outdir 2000 11 32
```

## Results, 2,000 replicates, seed 2026

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-linked-seed2026.txt"), String)))
```

## Results, 2,000 replicates, seed 7

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-linked-seed7.txt"), String)))
```

## Neutral and sweep, 2,000 replicates, seed 11

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-linked-seed11.txt"), String)))
```

## Sweep, 8,000 replicates, seed 17

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-linked-sweep-seed17.txt"), String)))
```

## Beneficial sites alone, 4,000 replicates, seed 13

```@eval
using Markdown, FisherWright
Markdown.MD(Markdown.Code("text", read(joinpath(pkgdir(FisherWright), "docs", "src", "validation", "validate-linked-ben-seed13.txt"), String)))
```

## [Reading the results](@id linked-reading)

**Neutral.** Over three seeds (6,000 replicates per simulator),
FisherWright's ``\pi`` agrees with SLiM's (pooled ``z = 0.06``) and
with ``\pi_0`` (``z`` between −0.45 and 0.79 per seed; ``N_e`` from
``\pi`` 199.9 to 200.2). The largest single differences, in seed 2026
(``\pi`` below 1 kb, ``z = -2.30``; neutral substitutions, ``z =
2.99``), do not recur in seeds 7 and 11.

**Background selection.** In two seeds, ``\pi`` over the genome
differs by ``z = -0.24`` and 0.04, and no statistic by more than
``|z| = 1.85``. Background selection lowers ``\pi`` to 0.95 of
``\pi_0`` (``N_e`` 190.5 to 191.0), to 0.93 near the genes; the load
``1 - \bar w`` is 0.116, and 15.5 to 15.7 deleterious sites fix per
run, in both simulators.

**Sweeps.** About 17 beneficial sites fix per run, and ``\pi`` falls to
0.875 of ``\pi_0`` (``N_e = 175``), from 0.82 near the genes to 0.90 at
50 kb and beyond. Over four seeds (14,000 replicates per simulator,
weighted by precision), FisherWright's ``\pi`` is 0.23% below SLiM's
(pooled ``z = -2.79``). All four seeds are negative (``z = -2.19``,
−2.66, −0.57, −0.96), and other statistics shift with it (neutral
substitutions +0.13%, ``z = 2.09``). The estimate fell as replicates
were added: two seeds of 2,000 gave −0.40%, the run of 8,000 −0.11%.

To locate a difference, the sweeps were run without the deleterious
sites. With the 40 beneficial sites alone, ``\pi`` differs by −0.01%
(``z = -0.08``) when they are linked (`ben`) and by −0.19%
(``z = -1.40``) when they are not (`benfree`); hitchhiking and
selection without linkage agree. If the residual under `sweep` is real,
it comes from the interplay of sweeps with linked deleterious sites;
a review of FisherWright's meiosis (crossovers, assortment, one strand
for neutral and selected sites), fitness, sampling of parents and
mutation found nothing to cause it. These tables hold about 60
comparisons, so a pooled ``z`` of −2.8 is not significant once that is
allowed for; its consistency over seeds is recorded here rather than
dismissed. In terms of effective size, it amounts to 0.4 of 175.

FisherWright's linked selection agrees with SLiM's: within noise under
neutrality, background selection and sweeps alone, and within 0.25% in
``\pi`` and ``N_e`` when sweeps and background selection act together.

## Scripts

```@eval
using Markdown, FisherWright
dir = joinpath(pkgdir(FisherWright), "docs", "src", "validation", "scripts")
Markdown.MD(Markdown.Code("text", read(joinpath(dir, "validate-linked.slim"), String)),
    Markdown.Code("julia", read(joinpath(dir, "validate-linked.jl"), String)),
    Markdown.Code("bash", read(joinpath(dir, "validate-linked.sh"), String)))
```
