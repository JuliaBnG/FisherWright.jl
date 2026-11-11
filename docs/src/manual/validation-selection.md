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

## CI test

The testset "Fixation probability under selection" runs four of these
cases (``N = 30``; ``s = 0.2`` with ``h`` = 0, ½ and 1, and ``s = 0``)
with 4,000 replicates each and compares them with the reference values
of the first table. The three values for ``s = 0.2`` are far apart, so
an error in the fitness of the heterozygote would be caught. Other CI
tests check the parts separately: the fitness of an individual against
the product over a dense genotype, the joint meiosis of neutral and
selected alleles against the recombinant of their union, the mutation
rate at the fitness loci, and the weighted draw of parents.

## Not covered

- Mutation–selection balance, the substitution rate under selection and
  the site frequency spectrum at selected sites.
- Linked selection: the effect of selected sites on neutral diversity
  (background selection, sweeps) against SLiM.
- Dominance in multi-locus runs, and deleterious alleles beyond the one
  case ``s = -0.1``.

## Scripts

[`fixation-reference.jl`](https://github.com/JuliaBnG/FisherWright.jl/blob/main/docs/src/validation/scripts/fixation-reference.jl)
and
[`validate-selection.jl`](https://github.com/JuliaBnG/FisherWright.jl/blob/main/docs/src/validation/scripts/validate-selection.jl).
Rerun both, and replace the two result files and the reference values in
`test/runtests.jl`, when the draw of parents, `SiteFitness`, or the
meiosis of the selected stream changes:

```bash
D=docs/src/validation
julia -t 16 $D/scripts/fixation-reference.jl 2000000 2026 > $D/fixation-reference.txt
julia -t 8 --project=. $D/scripts/validate-selection.jl 200000 2000000 2026 > $D/validate-selection.txt
```
