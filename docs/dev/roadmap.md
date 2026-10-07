# Roadmap: variable population size and non-random mating

Status: planned, not implemented (drafted 2026-10-07, against v0.3.8).

## Current model (v0.3.8)

- Diploid, dioecious, discrete non-overlapping generations, constant census
  N = `ne`.
- Sex is redrawn each generation as iid Bernoulli(½) (`rand!(rng, sex)` in
  `random_mate!`), so N_m ~ Binomial(N, ½), not exactly N/2. Then
  E[4N_mN_f/N] = N − 1: an O(1/N) bias, negligible at the validated `ne = 100`
  but relevant for small N or any sex-ratio feature.
- Each offspring draws its sire uniformly from males and its dam uniformly from
  females, with replacement. Family sizes are multinomial (≈ Poisson,
  V_k ≈ 2); selfing is impossible, and full-sib and half-sib matings occur at
  random rates.
- All working buffers are sized by `ne` once (`prt`, `off`, `mbuf`, `nbuf` of
  2N; `cbuf`, `pm`, `sex` of N). Variable N needs this restructured first.
- Every run starts monomorphic; no parentage is kept between generations.

## Design principles

1. **Default path is bit-identical.** With the same seed and thread count, the
   default call reproduces v0.3.8 output exactly. This is enforced by a
   regression test, so the default scheme must keep the current RNG draw order.
2. **Type stability through a function barrier.** `fisher_wright` stays the
   entry point and dispatches once to an internal
   `_simulate!(state, demog::D, mating::M)` specialised on concrete types. No
   dynamic dispatch inside the generation loop.
3. **Allocation-free inner loop.** Buffers are allocated at max N when the
   trajectory is known in advance; otherwise they grow geometrically with
   `resize!`. The zero-allocation tests are extended to the new mating schemes.

## Proposed API

```julia
fisher_wright(ne, nt, chr, mr;
    popsize = nothing,          # nothing (constant ne) | AbstractVector{<:Integer} (length nt) | g -> N_g
    mating  = RandomMating(),   # <: AbstractMating
    init    = nothing,          # Vector{Vector{UInt32}} or FisherWrightResult, to continue a run
    record  = :none,            # :none | :census | :pedigree
    …existing keywords…)
```

### Demography

Generation *g* uses N_{g−1} parents to produce N_g offspring.

- **Exogenous:** a vector or a function of *g* gives the census trajectory
  (bottleneck, exponential growth, step changes). Helper constructors such as
  `bottleneck(N0, Nb, t0, dur)` and `growth(N0, r)` build these trajectories.
- **Emergent:** when family sizes are drawn, N follows from them. This requires
  density regulation (a cap at K by random culling, or Beverton–Holt);
  otherwise a critical branching process goes extinct almost surely.
- **Extinction policy:** if N_g < 2 or one sex is missing, throw a typed error
  that reports *g*. Silent retries would bias the conditioning.

### Mating schemes

`AbstractMating`; each scheme implements `mate!(pm, sex, state, rng)`.

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

### Out of scope

- Phenotypic assortative mating: the simulator is neutral and has no trait
  model; this would need a separate trait layer.
- Overlapping generations.
- Subdivision and migration (a separate feature).
- Minimum-coancestry or optimum-contribution mating: O(N²) per generation;
  belongs with a breeding-programme layer.

## Results and bookkeeping

- `FisherWrightResult` gains optional fields: census trajectory, realised
  N_m and N_f, per-generation V_k for each sex, and a pedigree of the last *k*
  generations when `record = :pedigree`. Changing a public struct before 1.0
  calls for a minor-version bump.
- `init` allows running to equilibrium once and then applying a demographic
  event from that state. It is also needed for the heterozygosity-decay
  validation below.

## Validation

Each check is a z-score at the existing |z| ≤ 3.5 threshold, with replicates
and seeds; logs go under `docs/src/validation/`.

1. **Mating mechanics (unit tests in CI):** exact sex counts; realised
   family-size distribution matches the request (moment z-tests and χ²
   goodness of fit); realised selfing rate; zero relative-avoidance violations;
   the bit-identical default-path regression test.
2. **Ne under constant N, against theory:**
   - Start from a polymorphic population via `init`, run without new mutation,
     and fit H_t/H_0 = (1 − 1/(2Ne))^t. Requires a way to switch mutation off
     (currently `mr > 0.01` is enforced).
   - Equilibrium π and S against θ = 4Neμ with Ne from the formulas above.
   - Cases: skewed sex ratio (10♂/90♀ → Ne = 36), NegBinomial family sizes,
     equal families, selfing s ∈ {0.2, 0.5}.
3. **Variable N, against msprime DTWF with the same demography:** SFS, π and
   σ²_d per distance bin after a bottleneck and after an expansion, reusing
   the existing msprime harness. msprime DTWF is monoecious and Poisson-like,
   so this matches only `RandomMating`; the O(1/N) dioecy bias must be
   quantified or Ne mapped explicitly.
4. **Non-random schemes, against SLiM:** SLiM WF models with `mateChoice()` and
   `modifyChild()` callbacks express hierarchical designs, weighted parents and
   relative avoidance. Build matched scripts in `SLiM/`.
5. **Benchmarks:** re-run the 1 Gb transient benchmark. Release criterion: the
   default path is within run-to-run noise of v0.3.8; any regression is
   reported with numbers.

## Release staging

| Release | Content |
|---|---|
| **0.4.0** | Internal refactor (state struct, function barrier, resizable buffers), `init` continuation, exogenous `popsize`, fixed sex ratio, census and Ne metadata in results, bit-identical default path. Prerequisite for the rest. |
| **0.5.0** | `FamilySize`, `Hierarchical`, `WeightedParents`, density regulation and extinction policy; Hill-formula and msprime validation. |
| **0.6.0** | Pedigree recording, `Monoecious(selfing)`, `AvoidRelatives`, `CustomMating`; SLiM matched validation. |

## Open decisions

1. Mutation-free validation: relax the `mr` lower bound, or add a
   `mutate = false` validation keyword?
2. Default sex model: keep binomial sex for compatibility, or switch to exactly
   N/2 males in 0.4.0 (documented breaking change; makes Ne ≈ N exact)?
3. JOSS paper: mention this roadmap as future work, or keep it scoped to
   v0.3.x?
