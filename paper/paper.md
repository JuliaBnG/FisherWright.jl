---
title: 'FisherWright.jl: forward-time Wright–Fisher simulation of whole-genome founder populations in Julia'
tags:
  - Julia
  - population genetics
  - forward simulation
  - Wright–Fisher model
  - breeding simulation
  - genomic prediction
authors:
  - name: Xijiang Yu
    orcid: 0000-0001-6508-9251
    corresponding: true
    affiliation: 1
affiliations:
  - name: Department of Animal and Aquacultural Sciences, Faculty of Biosciences, Norwegian University of Life Sciences (NMBU), Ås, Norway
    index: 1
    ror: 04a1mvv97
date: 7 October 2026
bibliography: paper.bib
---

# Summary

Stochastic simulations of breeding programs and complex traits require
founder populations whose allele frequencies and linkage disequilibrium
(LD) determine marker–QTL architecture, genomic prediction accuracy, and
long-term selection response. `FisherWright.jl` generates such founders
by forward-time simulation of a finite, diploid Wright–Fisher
population. The model assumes a constant census size $N$, two sexes,
random union of a sire and dam sampled with replacement,
Poisson-distributed crossovers with independent assortment among
chromosomes, and neutral mutation under a biallelic finite-sites model.
Starting from monomorphic genomes, the simulation runs for a specified
number of generations, optionally moves fixed mutations out of the
active haplotypes, and exports segregating variants. Output formats
include a bit-packed `BitMatrix` with a linkage map, subsampled SNP-chip
panels, and native `BnGStructs.jl` haplotype containers designed for
downstream breeding simulations. `FisherWright.jl` is implemented in
pure Julia [@bezanson2017julia], MIT-licensed, and available through the
Julia General Registry.

# Statement of need

Simulations of genomic selection require founder genomes with site
frequency spectra (SFS) and LD profiles that arise from an explicit
evolutionary process. Naive Bernoulli sampling produces genotypes in
linkage equilibrium, whereas heuristic approximations, such as targeting
Sved's expectation $E[r^2] = 1/(1 + 4N_e c)$, do not reproduce the LD
decay of a finite population after minor allele frequency (MAF)
filtering. Coalescent simulators generate neutral equilibrium variation
efficiently, but their integration into Julia workflows typically
requires inter-process communication, foreign-function wrappers, or
file exchange with external C and Python executables.

`FisherWright.jl` addresses this gap with an in-process, pure-Julia
founder engine. It eliminates cross-language dependencies, returns the
data types used by downstream breeding pipelines, and provides a
transparent forward-simulation loop that can be inspected, extended,
and taught without external tooling.

# State of the field

Existing genomic simulation software spans several paradigms. Coalescent
engines such as `msprime` [@kelleher2016efficient;
@baumdicker2022efficient] and `MaCS` [@chen2009fast] are the fastest
means of generating neutral equilibrium genomes; `MaCS`, for example,
is the default founder generator in `AlphaSimR`
[@gaynor2021alphasimr]. General forward-time simulators such as `SLiM`
[@haller2019tree; @haller2026slim] and `fwdpy11`
[@thornton2019polygenic] accommodate arbitrary selection, complex
demography, and spatial structure via the Eidos scripting language or
Python/C++ interfaces. Specialized breeding simulators, including
`QMSim` [@sargolzaei2009qmsim], `MoBPS` [@pook2020mobps], and `XSim`
[@cheng2015xsim; @chen2022xsim], focus on pedigree management, mating
designs, and multi-generation selection programs. Among these, `XSim`
version 2 is also implemented in Julia.

`FisherWright.jl` occupies a distinct niche through its deliberately
limited scope. It does not implement selection, population structure, or
non-stationary demography; workflows requiring those features remain
better served by `SLiM` or `fwdpy11`. For neutral equilibrium founders,
`msprime` is substantially faster when a Julia-only workflow is not
needed. We built a new package rather than contribute to these tools for
three reasons. First, the `JuliaBnG` packages needed a founder generator
without non-Julia dependencies that produces the `BnGStructs.jl`
haplotype type used throughout the pipeline. Second, it needed to be
small enough to validate thoroughly against population-genetic theory.
Third, it was refactored from the author's livestock breeding simulator
`xyBnG.jl` [@yu2026xybng], which is not publicly available, so that the
founder step could be released, tested, and cited independently.

# Software design

Haplotypes are represented as sorted `Vector{UInt32}` coordinates of
derived alleles on a concatenated genome of up to $2^{32} - 1$ bp. This
sparse representation avoids the prohibitive memory footprint of dense
genotype arrays at gigabase scale. Tree-sequence recording would improve
the efficiency of forward simulation, but would introduce an additional
dependency and a second data model.

Meiosis is implemented as a single two-pointer traversal of the two
parental haplotypes, switching strands at each crossover breakpoint
(`recombine`). De novo mutations are added by a sorted merge into a
preallocated buffer (`merge_sorted!`). These kernels and the crossover
sampler (`cobp!`) write to preallocated thread-local buffers; regression
tests verify that they perform no heap allocations after warm-up. Within
each generation, meioses and mutation draws run in parallel using
`Threads.@threads`; results are reproducible for a fixed random seed and
thread count.

Under neutrality, fixed derived alleles accumulate at the substitution
rate $\mu L$ per generation across $L$ sites and, if retained,
progressively reduce meiosis throughput. When tracking is enabled
(`result = true`), a sweep at a user-defined interval (every generation
by default) removes positions fixed across all $2N$ haplotypes and
records them in the substitution log of a `FisherWrightResult` struct.
The active haplotype vectors then contain only segregating sites, which
bounds memory consumption and traversal overhead during long burn-in
periods.

Data export is decoupled from forward simulation. `muts2bitarray`
converts polymorphic sites to a 1-bit-per-allele `BitMatrix` with an
accompanying `DataFrame` linkage map. For large genomes, where dense
matrices are unwieldy, `extract_chip_bitarray` and
`to_haplotype(result, positions)` evaluate a sorted subset of marker
coordinates directly against sparse haplotypes to generate target SNP
panels without instantiating intermediate whole-genome arrays.

# Validation

The validation suite compares `FisherWright.jl` (v0.3.8) with
analytical neutral expectations and matched simulations from `msprime`'s
discrete-time Wright–Fisher (DTWF) model. The validation setting models
$N = 100$ individuals across 10 chromosomes of 10 Mb each ($L = 100$
Mb), evolving for $20N = 2{,}000$ generations at $\mu = r = 10^{-8}$ per
base pair per generation, with 50 independent replicates per simulator.

The suite comprises 24 prespecified statistical checks against theory,
`msprime`, or both: crossover counts per meiosis (the mean and
variance-to-mean ratio for each chromosome) and independent assortment;
segregating sites ($S$) against Watterson's $\theta a_n$; gene diversity
($\sum 2pq$) against $\theta(1 - 1/n)$; whole-population and sampled
($n = 20$) site frequency spectra; and LD decay ($\sigma^2_d$, MAF
$\geq 0.1$) across five distance bins from 10 kb to 3 Mb. A check fails
when $|z| > 3.5$, where $z$ is calculated from means over the 50
replicates. All 24 checks passed in each of two independent runs (seeds
7 and 2026; maximum $|z|$ values of 2.24 and 2.49, respectively).

A lightweight statistical suite runs in continuous integration for every
commit. It passes for five seeds and fails when the simulated $N$ is
halved, the recombination rate is doubled, or an earlier bug in
final-generation mutation retention (fixed in v0.3.5) is reintroduced.
The online documentation includes the complete validation scripts, raw
output logs, and a comparison with release v0.1.6.

# Research impact statement

`FisherWright.jl` is the founder-population generator for the
`JuliaBnG` ecosystem and directly produces `BnGStructs.jl` haplotype
containers. An earlier release (v0.1.x) simulated the historical base
population for an evaluation of enteric methane breeding objectives in
the Norwegian White Sheep program [@yu2026effects]. That release
predates the corrections listed in the package changelog; their effects
on the founders are documented in the manual.

The package repository includes reproducible benchmark drivers. In one
run with $N = 2{,}000$ individuals evolved for 2,000 generations
($t = N$, before equilibrium) on a 1 Gb genome ($10 \times 100$ Mb
chromosomes), forward simulation required 59.8 s on 8 threads of an AMD
Ryzen 9 3900X. Converting the resulting 567,975 segregating sites across
4,000 haplotypes to a 271 MiB `BitMatrix` required a further 9.0 s. By
comparison, `SLiM` with tree-sequence recording required 51.6 s for
forward simulation and 11.2 s for variant export, whereas `msprime`
using its DTWF model completed the equivalent steps in 9.6 s and 7.6 s,
respectively.

# AI usage disclosure

Generative AI tools were used during development and writing:

- GitHub Copilot, for code suggestions and for commits after v0.1.6;
  some of these list Copilot as co-author in the Git history, which was
  corrected in later commits;
- GPT-5, for an early code review;
- Gemini, for drafts of earlier statistical tests, which were later
  replaced, and for manuscript editing;
- Google Antigravity, for manuscript wording and code review;
- Claude (Anthropic), for code and manuscript review, validation and
  comparison scripts, documentation, and a first draft of this paper.

The author reviewed, edited, and validated all AI-assisted output and
made the design decisions. Every reported number is drawn from logged
runs in the repository.

# Acknowledgements

The author acknowledges financial support from the Green ERA-Hub
programme within the EU’s Horizon Europe research and innovation
programme (grant no. 101056828) and the Research Council of Norway
(grant no. 351087).

# References
