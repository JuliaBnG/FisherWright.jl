# Roadmap: variable population size, non-random mating, and selection

Status: in progress (drafted 2026-10-07 against v0.3.8; selection added
2026-11-02 against v0.4.0). v0.4.0 shipped without the refactor below,
so the staging is renumbered from 0.5.0. Implemented so far, all from
the 0.5.0 row: the golden-output regression test and the `_recombine!`
kernel (v0.4.2); the `mutate` and `init` keywords, the internal
refactor, and exogenous `popsize` (2026-11-10); the `mating` keyword
with a fixed number or proportion of males, and the census record in
results with `effective_size` (2026-11-11). None of the 2026-11 work is
released yet.

The refactor (`src/simulate.jl`): a `SimState` struct holds the
haplotypes, buffers and per-slot streams; `fisher_wright` builds it and
calls `_simulate!(st, demog, mating, fitness, …)` once, specialised on
the concrete types. Buffers are allocated for the largest census size
of the run, and the current cohort occupies their leading entries.
Demography is `ConstantSize` or `SizeTrajectory` (from `popsize`, a
vector or a function evaluated before the run). Mating is
`RandomMating`, exported and passed as `mating`: binomial sex by
default, or exactly `nmale` males, or a proportion `sexratio` of males,
drawn as a uniformly random subset of each parental cohort. Fitness has
only its default (`Neutral`); the demography and fitness types are not
exported. `FisherWrightResult` carries the census record of the run
(`census`, N_0 to N_nt, and `males`, the male parents of each
generation), appended to that of an `init` result; `effective_size`
returns 4N_mN_f/(N_m+N_f) per generation
or its harmonic mean. Still to do in 0.5.0: the trajectory helpers
(`bottleneck`, `growth`) and validation 3 against msprime.

## Current model (v0.4.0)

- Diploid, dioecious, discrete non-overlapping generations, constant
  census N = `ne`.
- Sex is redrawn each generation as iid Bernoulli(½) (`rand!(rng, sex)`
  in `random_mate!`), so N_m ~ Binomial(N, ½), not exactly N/2. Then
  E[4N_mN_f/N] = N − 1: an O(1/N) bias, negligible at the validated
  `ne = 100` but relevant for small N or any sex-ratio feature.
- Each offspring draws its sire uniformly from males and its dam
  uniformly from females, with replacement. Family sizes are multinomial
  (≈ Poisson, V_k ≈ 2); selfing is impossible, and full-sib and half-sib
  matings occur at random rates.
- Strictly neutral. A mutation is identified only by its `UInt32`
  position: sites are biallelic, a repeat hit is a no-op, and there is
  no back-mutation. There is no per-mutation effect to select on.
- Order of events in generation *g*: mating on `prt`, meiosis, mutation.
  `prt` at the start of a generation is therefore the adult cohort,
  which is where fitness must be evaluated.
- Mating is drawn serially from `rng`; meiosis and mutation use one
  `Xoshiro` stream per slot (`srng`), seeded once from `rng`. A run is
  independent of the thread count.
- All working buffers are sized by `ne` once (`prt`, `off`, `mbuf`,
  `nbuf` of 2N; `cbuf`, `pm`, `sex` of N). Variable N needs this
  restructured first.
- Every run starts monomorphic; no parentage is kept between
  generations.

## Design principles

1. **Default path is bit-identical.** With the same seed, the default
   call (`RandomMating()`, constant N, `Neutral()`) reproduces v0.4.0
   output exactly. This is enforced by a regression test, so the default
   scheme must keep the current RNG draw order.
2. **Type stability through a function barrier.** `fisher_wright` stays
   the entry point and dispatches once to an internal
   `_simulate!(state, demog::D, mating::M, fitness::F)` specialised on
   concrete types. No dynamic dispatch inside the generation loop;
   `Neutral()` compiles to the current loop.
3. **Allocation-free inner loop.** Buffers are allocated at max N when
   the trajectory is known in advance; otherwise they grow geometrically
   with `resize!`. The zero-allocation tests are extended to the new
   mating schemes and fitness models.

## Proposed API

```julia
fisher_wright(ne, nt, chr, mr;
    popsize = nothing,          # nothing (constant ne) | AbstractVector{<:Integer} (length nt) | g -> N_g
    mating  = RandomMating(),   # <: AbstractMating
    fitness = Neutral(),        # <: AbstractFitness
    init    = nothing,          # Vector{Vector{UInt32}} or FisherWrightResult, to continue a run
    record  = :none,            # :none | :census | :pedigree
    …existing keywords…)
```

### Demography

Generation *g* uses N_{g−1} parents to produce N_g offspring.

- **Exogenous:** a vector or a function of *g* gives the census
  trajectory (bottleneck, exponential growth, step changes). Helper
  constructors such as `bottleneck(N0, Nb, t0, dur)` and `growth(N0, r)`
  build these trajectories.
- **Emergent:** when family sizes are drawn, N follows from them. This
  requires density regulation (a cap at K by random culling, or
  Beverton–Holt); otherwise a critical branching process goes extinct
  almost surely.
- **Extinction policy:** if N_g < 2 or one sex is missing, throw a typed
  error that reports *g*. Silent retries would bias the conditioning.

### Mating schemes

`AbstractMating`; each scheme implements
`mate!(pm, sex, w, state, rng)`, where `w` is the vector of relative
fitnesses (all ones, and ignored, under `Neutral()`).

| Scheme | Mechanism | Expected Ne (discrete generations) |
|---|---|---|
| `RandomMating()` | Current behaviour, the default | ≈ N − 1 under binomial sex |
| `RandomMating(; nmale, sexratio)` | Fixed N_m, N_f or a fixed proportion | 4N_mN_f/(N_m+N_f) (Wright) |
| `FamilySize(dist_f; dist_m)` | Offspring counts per dam (and sire usage) from Poisson, NegBinomial(r, p), or a constant (equal family size) | 8N/(V_m+V_f+4) (Hill 1972); equal families ≈ 2N |
| `Hierarchical(nsire, dams_per_sire, offspring_per_dam)` | Livestock-style nested design | Hill/Gowe formulas; compare against SLiM |
| `WeightedParents(w_m, w_f)` | Skewed parental usage (Gamma or Zipf weights, or user-supplied) | Through V_k |
| `Monoecious(; selfing = s)` | Hermaphrodite mode with selfing rate s | N/(1+F), F = s/(2−s) |
| `AvoidRelatives(inner; level = :fullsib/:halfsib)` | Rejection sampling over a 1–2 generation pedigree | Design-dependent; empirical |
| `CustomMating(f)` | User callback `(g, sex, pedigree, rng) -> pm` | User's responsibility |

`WeightedParents` and selection share one sampler: a preallocated Vose
alias table, rebuilt per sex per generation in O(N), with O(1) draws.
`Distributions.sampler(Categorical(p))` allocates on every call and is
not used in the loop.

### Selection

Soft selection in the Wright–Fisher sense: fitness is evaluated on the
adult cohort `prt`, and parents are drawn with probability proportional
to relative fitness w within each sex. N is unaffected. Hard selection
(absolute fitness reducing N, extinction) belongs with emergent
demography.

**Fitness loci.** Selected sites form an explicit, static table built
once per run:

```julia
struct FitnessLoci
    pos::Vector{UInt32}    # sorted, unique selected sites
    s::Vector{Float64}     # derived-allele selection coefficient (or additive effect a)
    h::Vector{Float64}     # dominance
    μ::Float64             # per-site mutation rate
end
```

Constructors take an explicit list (known QTL, candidate genes) or draw
sites and effects: `fitness_loci(chr, nsites; dfe, h, regions, rng)`,
with `h` either constant or a function h(s). A static table rather than
a per-mutation draw keeps effects consistent across the per-slot RNG
streams, makes runs reproducible, and lets the true causal effects be
exported for GWAS and genomic-prediction studies. Memory is about 16 B
per site. Should a genome-wide DFE ever be needed, the fallback is a
counter-based draw `s = quantile(dfe, u(hash(pos, seed)))`, computed
once per new mutation and cached.

**Two mutation streams.** Selected alleles live in a second
per-haplotype vector `sel` (with `soff`, `smbuf`, `snbuf`), holding
positions of carried derived alleles at fitness loci. Neutral mutations
stay in `prt` with the existing genome-wide model; the streams never
mix, so a neutral mutation at a selected coordinate is unambiguous.

- Meiosis applies the same crossover list and the same starting strand
  to both streams. `recombine` currently draws the strand internally
  (`o = rand(rng, Bool)`); it is split into a kernel
  `_recombine!(h1, h2, ho, co, o::Bool)`. The public `recombine` draws
  `o` and calls the kernel, so its behaviour and the default RNG order
  are unchanged.
- Mutation at fitness loci: per haplotype, `Poisson(K·μ)` hits on sites
  drawn uniformly from `1:K`, merged with `merge_sorted!` as now.
  Already carried sites are no-ops (biallelic, no back-mutation).
- Fitness evaluation scans only `sel`, which is far shorter than `prt`.
  Intersecting each full haplotype with `pos` would cost O(|hap|) per
  generation, comparable to recombination itself.

**Fitness evaluation.** Threaded over individuals, allocation-free. For
individual *i*, merge-walk `sel[2i-1]` and `sel[2i]`: a site on both
haplotypes adds `log1p(s)`, a site on one adds `log1p(h*s)`; effects are
looked up with `searchsortedfirst(pos, x)`. Then
`w[i] = exp(lw[i] - maximum(lw))`.

**Fitness models.**

```julia
abstract type AbstractFitness end
struct Neutral <: AbstractFitness end                            # default; compiles away
struct SiteFitness <: AbstractFitness; loci::FitnessLoci; end    # multiplicative DFE
struct TraitFitness{F} <: AbstractFitness                        # QTL → phenotype → w
    loci::FitnessLoci   # s read as additive effect a, h as dominance
    σe::Float64
    wfun::F             # Gaussian stabilizing (θ, ω²), directional exp(βP), …
end
```

`TraitFitness` is the quantitative-genetic layer: G = Σ dosage·a (plus
dominance), P = G + e with e ~ N(0, σe²), w = f(P). Environmental draws
come from a serial stream so the run stays thread-independent.
Pleiotropy (several traits on shared loci) and sex-specific fitness are
later extensions.

**Fixation bookkeeping** (`result = true`).

- `SiteFitness`: a fixed site multiplies every w by the same constant,
  which cancels on normalisation, so removing it is exact. Selected
  substitutions are recorded separately as `(pos, s, h, g_fixed)`.
- `TraitFitness`, or any non-multiplicative combination across loci: a
  fixation shifts G for everyone, which matters relative to an optimum
  θ. A `fixed_offset::Float64` accumulates the fixed effects and is
  added to G.
- Lost selected alleles need no maintenance: the site table is static.

### Out of scope

- Overlapping generations.
- Subdivision and migration (a separate feature).
- Minimum-coancestry or optimum-contribution mating, and truncation
  selection on EBVs: O(N²) or model-fitting per generation; belongs with
  a breeding-programme layer (Breeding.jl).
- Phenotypic assortative mating: becomes possible once `TraitFitness`
  provides phenotypes, but is not scheduled.

## Results and bookkeeping

- `FisherWrightResult` gains optional fields: census trajectory,
  realised N_m and N_f, per-generation V_k for each sex, and a pedigree
  of the last *k* generations when `record = :pedigree`. Changing a
  public struct before 1.0 calls for a minor-version bump. The census
  trajectory and realised N_m are in (`census`, `males`), always
  recorded; the three-argument constructor is kept. V_k and the
  pedigree wait for the `record` keyword and the schemes that need
  them.
- Under selection it also carries the `FitnessLoci` table, the active
  `sel` haplotypes, selected substitutions, `fixed_offset`, and a
  per-generation trajectory of mean and variance of log w (and of G, P
  under `TraitFitness`). `to_haplotype` and `muts2bitarray` gain a way
  to export the selected stream with its true effects.
- `init` allows running to equilibrium once and then applying a
  demographic event or switching selection on from that state. It is
  also needed for the heterozygosity-decay and fixation-probability
  validations below.

## Validation

Each check is a z-score at the existing |z| ≤ 3.5 threshold, with
replicates and seeds; logs go under `docs/src/validation/`.

1. **Mating mechanics (unit tests in CI):** exact sex counts; realised
   family-size distribution matches the request (moment z-tests and χ²
   goodness of fit); realised selfing rate; zero relative-avoidance
   violations; the bit-identical default-path regression test. Done for
   the fixed-sex schemes ("Mating with fixed sex numbers": exact counts,
   and χ² tests that the males are a uniform random subset and that
   sires and dams are used uniformly).
2. **Ne under constant N, against theory:**
   - Start from a polymorphic population via `init`, run without new
     mutation (`mutate = false`), and fit H_t/H_0 = (1 − 1/(2Ne))^t.
     Done for `RandomMating()` as a CI test against the exact two-sex
     recursion ("Heterozygosity decay without mutation", N = 20); the
     other schemes follow with their implementation.
   - Equilibrium π and S against θ = 4Neμ with Ne from the formulas
     above.
   - Cases: skewed sex ratio (10♂/90♀ → Ne = 36), NegBinomial family
     sizes, equal families, selfing s ∈ {0.2, 0.5}. The sex-ratio case
     is done as a CI test of heterozygosity decay against the exact
     recursion ("Heterozygosity decay under a fixed sex ratio": 3♂/17♀,
     10♂/10♀, and a proportion of 0.2 males through a change in size);
     its equilibrium π and S check is still open.
3. **Variable N.** Done as a CI test: heterozygosity through a
   bottleneck and an expansion against the exact two-sex recursion
   (30 → 6 for 8 generations → 60 for 12), which also shows that the
   drift of generation g follows its N_{g−1} parents. Still to do,
   **against msprime DTWF with the same demography:** SFS, π
   and σ²_d per distance bin after a bottleneck and after an expansion,
   reusing the existing msprime harness. msprime DTWF is monoecious and
   Poisson-like, so this matches only `RandomMating`; the O(1/N) dioecy
   bias must be quantified or Ne mapped explicitly.
4. **Non-random schemes, against SLiM:** SLiM WF models with
   `mateChoice()` and `modifyChild()` callbacks express hierarchical
   designs, weighted parents and relative avoidance. Build matched
   scripts in `SLiM/`.
5. **Selection mechanics (unit and Supposition property tests in CI):**
   fitness from the merge-walk equals a brute-force calculation on dense
   genotypes; realised parent usage matches w/Σw (χ²); `_recombine!` on
   both streams equals `recombine` on their union, split afterwards.
6. **Selection against single-locus theory:**
   - Fixation probability of one introduced copy (via `init`) against
     Kimura's diffusion result: closed form at h = ½, numerical integral
     for general h; 4Ns ∈ {−2, 0, 2, 10}, about 10⁵ replicates, binomial
     confidence intervals. Use the dioecious Ne (≈ N − 1 under binomial
     sex), not N.
   - Mutation–selection balance: w̄ ≈ e^{−U} (Haldane–Muller, h > 0,
     loosely linked sites) and per-site q ≈ μ/(hs).
   - Substitution-rate ratio K_sel/K_neu = 2N·u(1/2N).
   - SFS at selected sites against the Poisson random field
     (Sawyer–Hartl) expectation per 4Ns class.
7. **Linked selection, against SLiM:** extend `bench/simulate-slim.slim`
   with an `m2` mutation type and genomic elements at the same sites,
   the same DFE, and `mutationStackPolicy = "f"`. Compare neutral π and
   S as a function of distance to selected sites (background selection
   and sweeps, beyond the range of Hudson–Kaplan and Nordborg
   approximations under strong Hill–Robertson interference), load, and
   substitution counts. SLiM draws s per mutation while this design
   fixes s per site: choose K so that recurrent hits at segregating
   sites are negligible, and report that rate. Under selection,
   offspring-number variance is no longer Poisson; report Ne estimated
   from neutral sites rather than reusing the neutral validation.
8. **Benchmarks:** re-run the 1 Gb transient benchmark. Release
   criterion: the default path is within run-to-run noise of v0.4.0; any
   regression is reported with numbers. Report selection-path overhead
   as a function of K and U.

## Release staging

| Release | Content |
|---|---|
| **0.5.0** | Internal refactor (state struct, function barrier over demography, mating and fitness, resizable buffers), `_recombine!` kernel split, `init` continuation, mutation switch for validation, exogenous `popsize`, fixed sex ratio, census and Ne metadata in results, bit-identical default path. Prerequisite for the rest. |
| **0.6.0** | `FitnessLoci`, `SiteFitness`, `sel` stream, alias-table weighted mating, selected substitutions; validations 5–6. |
| **0.7.0** | `TraitFitness` (stabilizing and directional), true-effect export; `FamilySize`, `Hierarchical`, `WeightedParents`, density regulation and extinction policy; Hill-formula and msprime validation. |
| **0.8.0** | Pedigree recording, `Monoecious(selfing)`, `AvoidRelatives`, `CustomMating`; SLiM matched validation, including linked selection (validation 7). |

## Open decisions

1. Mutation-free validation: resolved 2026-11-02 in favour of a
   `mutate = false` keyword; the `mr` bounds are unchanged and still
   checked when mutation is off.
2. Default sex model: keep binomial sex for compatibility, or switch to
   exactly N/2 males in 0.5.0 (documented breaking change; makes Ne ≈ N
   exact)?
3. JOSS paper: mention this roadmap as future work, or keep it scoped to
   v0.4.x?
4. Selected-site semantics: fixed s per site (this design) or s per
   mutation as in SLiM? Per-site is reproducible and exportable;
   per-mutation needs an effect table that grows and is pruned at loss.
5. Should selection precede the remaining mating schemes (as staged), or
   follow them? They share the alias sampler and the `mate!` signature,
   so either order works once 0.5.0 lands.
