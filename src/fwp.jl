_mean_length(v) = isempty(v) ? 0.0 : sum(length, v) / length(v)

function _fmt_seconds(t::Real)
    s = round(Int, t)
    h, s = divrem(s, 3600)
    m, s = divrem(s, 60)
    return h > 0 ? string(h, "h", lpad(m, 2, '0'), "m") :
           m > 0 ? string(m, "m", lpad(s, 2, '0'), "s") : string(s, "s")
end

# One progress line on stderr. On a terminal the line is redrawn in place;
# otherwise (log files, batch jobs) each report goes on its own line.
function _report_progress(io::IO, g, nt, t0, prt, nsub, result::Bool)
    elapsed = time() - t0
    eta = elapsed / g * (nt - g)
    msg = string(
        "Generation ", g, " / ", nt,
        " | mean muts/hap: ", round(_mean_length(prt); digits=2),
        result ? string(" | substitutions: ", nsub) : "",
        " | elapsed: ", _fmt_seconds(elapsed),
        " | ETA: ", _fmt_seconds(eta),
    )
    if io isa Base.TTY
        print(io, "\r\e[K", msg)
        g == nt && println(io)
    else
        println(io, msg)
    end
    flush(io)
    return nothing
end

"""
    fisher_wright(
        ne::Integer,
        nt::Integer,
        chr::Vector{<:Integer},
        mr::Float64;
        M = 1e8,
        mut_base = 1e8,
        result::Bool = false,
        fixation_interval::Int = 1,
        mutate::Bool = true,
        init = nothing,
        popsize = nothing,
        mating = RandomMating(),
        fitness = Neutral(),
        verbose::Bool = false,
        progress_interval::Integer = 100,
        rng::AbstractRNG = Random.default_rng(),
    ) -> Tuple{Vector{Vector{UInt32}}, Vector{UInt32}} | FisherWrightResult

Simulate a diploid Fisher-Wright population of size `ne` over `nt` generations under
random mating, crossover recombination, and recurrent mutation, optionally with a changing
census size, a fixed sex ratio, and selection at given sites.

Mutations are tracked as sorted, unique vectors of `UInt32` genomic coordinates, supporting
total genome lengths up to ``2^{32} \\approx 4.29 \\times 10^9`` base pairs (bp).

# Arguments
- `ne::Integer`: Diploid population size (number of individuals, `ne > 1`), constant unless `popsize` is given; `2 * ne` haplotypes are simulated. With `popsize`, `ne` is the size of the starting population.
- `nt::Integer`: Number of generations to simulate (`nt > 0`).
- `chr::Vector{<:Integer}`: Chromosome lengths in base pairs (all entries must be positive). Total genome length `sum(chr)` must be `< 2^32`.
- `mr::Float64`: Expected mutation rate per `mut_base` base pairs per meiosis (must satisfy `0.01 < mr < 20.0`).

# Keywords
- `M = 1e8`: Number of base pairs per Morgan (`M > 0`). A chromosome of length `L` experiences an expected `L / M` crossovers per meiosis.
- `mut_base = 1e8`: Base pairs per unit of mutation rate `mr` (`mut_base > 0`). Decoupled from `M` so adjusting recombination does not alter mutation intensity.
- `result::Bool = false`:
  - If `false` (default), returns `(prt, cbp)` where `prt` is the vector of `2 * ne` haplotypes and `cbp` is cumulative chromosome end coordinates. Fixed mutations remain in all haplotypes.
  - If `true`, returns a [`FisherWrightResult`](@ref), periodically extracting and tracking fixed substitutions to save memory and runtime.
    For the same `rng` state this is the same realization as `result = false`: `FisherWright._with_fixed(res)` equals `prt`.
- `fixation_interval::Int = 1`: Frequency (in generations) to scan for and extract fixed substitutions when `result = true`.
- `mutate::Bool = true`: If `false`, no new mutations are generated: the population evolves by drift and recombination only, and `mr` and `mut_base` are not used (`mr` must still lie in its valid range). Mating is drawn from `rng` exactly as with `mutate = true`, but meiosis shares the per-slot streams with the skipped mutation draws, so for the same seed the two settings give different realizations. Intended for validation runs that start from a polymorphic population given by `init`, for example heterozygosity decay ``H_t/H_0 = (1 - 1/(2N_e))^t``; without `init`, a run with `mutate = false` stays monomorphic.
- `init = nothing`: Starting population; `nothing` starts monomorphic. Either a `Vector{Vector{UInt32}}` of `2 * ne` haplotypes (each sorted, unique, with positions in `1:sum(chr)`), or a [`FisherWrightResult`](@ref) whose `chromosome_ends` equal `cumsum(chr)`. Haplotypes `2i - 1` and `2i` form individual `i`, the convention of the output, and these individuals are the parents of generation 1. `init` is copied and left unmodified. The substitutions of a `FisherWrightResult` are carried over when `result = true` and merged back into every haplotype when `result = false`; given as plain haplotypes, sites carried by all of them are extracted at the first fixation scan when `result = true`. See *Continuing a run* below.
- `popsize = nothing`: Census size per generation; `nothing` keeps it at `ne`. Either a vector of `nt` integers or a function `g -> N_g` of the generation `g = 1, …, nt`, which is called once per generation before the run; every size must exceed 1. Generation `g` draws the parents of its `N_g` offspring from the `N_{g-1}` individuals of the previous generation (`N_0 = ne`), so the run returns `2 * N_nt` haplotypes. [`bottleneck`](@ref) and [`growth`](@ref) build common trajectories. See *Changing population size* below.
- `mating = RandomMating()`: Mating scheme. The default redraws sex as Bernoulli(½) per individual in every generation; [`RandomMating`](@ref)`(nmale = m)` or `RandomMating(sexratio = p)` fixes the number or the proportion of males instead. See *Sex ratio* below. [`FamilySize`](@ref) controls the distribution of family sizes, with one sire per dam, and [`Hierarchical`](@ref) is the nested design of `s` sires each mated to `d` dams.
- `fitness = Neutral()`: Fitness model. [`SiteFitness`](@ref)`(loci)` selects on the sites of a [`FitnessLoci`](@ref) table, and [`TraitFitness`](@ref)`(loci; σe, fitness)` on a quantitative trait with the loci of a [`TraitLoci`](@ref) table; each draws parents with probability proportional to fitness within each sex (with [`FamilySize`](@ref), multiplies the split weights by fitness), and requires `result = true`. See *Selection* below.
- `verbose::Bool = false`: If `true`, writes a start banner and a progress line to `stderr` every `progress_interval` generations and at the last generation, reporting mean mutations per haplotype, accumulated substitutions (when `result = true`), elapsed wall time, and a linear ETA. On a terminal the line is updated in place; when `stderr` is redirected each report is a separate line.
- `progress_interval::Integer = 100`: Generations between progress reports when `verbose = true` (must be positive).
- `rng::AbstractRNG = Random.default_rng()`: Random number generator; see *Reproducibility* below.

# Returns
- When `result = false`:
  - `prt::Vector{Vector{UInt32}}`: Active haplotypes (`2 * ne`, or `2 * N_nt` with `popsize`), each containing sorted, unique 1-based genomic mutation positions.
  - `cbp::Vector{UInt32}`: Cumulative chromosome end positions in base pairs (`cumsum(chr)`).
- When `result = true`:
  - A [`FisherWrightResult`](@ref) containing `active_haplotypes`, `chromosome_ends`, accumulated `substitutions`, and the record of the run: the census size of every generation (`census`) and the realized number of males among its parents (`males`), from which [`effective_size`](@ref) gives the realized effective size. With selection it also holds the selected alleles (`selected_haplotypes`, `selected_substitutions`) and, per generation, the mean fitness of the parents and the mean and variance of its logarithm (`mean_fitness`, `mean_log_fitness`, `var_log_fitness`).

# Details
- **Recombination & Independent Assortment**: Recombination is modeled as a Poisson process along each chromosome with rate `chr[i] / M`. Independent assortment between chromosomes occurs with probability 0.5 at each chromosome boundary.
- **Mutation model**: Each haplotype receives `Poisson(sum(chr) / mut_base * mr)` new mutations per generation at uniformly drawn positions. A new mutation at a position the haplotype already carries, or at a substituted (fixed and removed) position when `result = true`, has no effect, so both return modes describe the same population.
- **Multithreading & Reproducibility**: Meiosis and mutation generation are parallelized using `Threads.@threads`. Each meiosis/haplotype slot owns a `Xoshiro` stream seeded once from `rng` at the start, and mating is sampled serially from `rng`, so a run is determined by the state of `rng` alone, independent of `Threads.nthreads()` and task scheduling.
  This holds for a fixed Julia version and package environment. The realization also depends on the variate algorithms of Julia's `Random` and of Distributions.jl, which can change between versions: `randn` differs between Julia 1.10 and 1.13, so Poisson draws with mean ≥ 6 (Distributions' `PoissonADSampler`; in practice the per-haplotype mutation count of most genome-scale runs) can differ, and the same seed then gives a different population. Record the Julia version and `Manifest.toml` with simulated data.
- **Continuing a run**: With `init`, the population evolves for `nt` further generations from the given state, for example to burn in once and then start several replicates, or to follow a set of standing variants with `mutate = false`. The per-slot streams are seeded anew from `rng` on every call, so two calls of `nt₁` and `nt₂` generations are a valid realization of `nt₁ + nt₂` generations but not the one a single call would give for the same seed. When `init` is a `FisherWrightResult` and `result = true`, the census record of the new generations is appended to that of `init`, so the result describes all `nt₁ + nt₂` generations.
- **Changing population size**: With `popsize`, the census size follows the given trajectory, for example a bottleneck or an expansion; mating, meiosis and mutation are otherwise as described above, and sex is redrawn in every generation, so the effective size of generation `g` follows the `N_{g-1}` parents. Working buffers and the per-slot streams are allocated once for the largest size of the run, so the realization for a given seed depends on that maximum: a trajectory that stays at `ne` reproduces the run without `popsize`, while two trajectories with different maxima share no random numbers.
- **Sex ratio**: Sex is not inherited: it is assigned anew to each cohort when it reproduces, and every offspring draws a sire and a dam with replacement. By default the number of males among `N` parents is Binomial(`N`, ½) conditioned on both sexes being present, which gives ``N_e \\approx N - 1``. With `mating = RandomMating(nmale = m)` exactly `m` individuals, a uniformly random subset, are males in every generation, and with `RandomMating(sexratio = p)` exactly `round(Int, p * N)` (at least one of each sex); the effective size is then ``4 N_m N_f / (N_m + N_f)``, so `sexratio = 0.5` gives ``N_e \\approx N`` for even `N`. The fixed-sex schemes draw sex differently from `rng`, so they share no realization with the default for the same seed.
- **Selection**: With `fitness = SiteFitness(loci)`, the sites of `loci` are under selection: a genotype has fitness `1`, `1 + h*s` or `1 + s` with zero, one or two copies of the derived allele, multiplied over sites. Selection is soft and acts within each sex: every offspring draws its sire among the males and its dam among the females with probability proportional to fitness, with replacement, so the census size is unchanged and only relative fitness within a sex matters. The derived alleles at the fitness loci are kept apart from the neutral mutations: they are returned in `selected_haplotypes` and, once fixed, `selected_substitutions`, they arise by mutation at rate `loci.μ` per site per haplotype per generation (never by the neutral mutation process, which may hit the same coordinates without effect on fitness), and they recombine with the neutral mutations as one haplotype. A selected allele of the starting population is given through `init`, as a `FisherWrightResult` with `selected_haplotypes`. If no member of a sex has positive fitness the run stops with an error. For the same seed, a run with selection shares no realization with a neutral run. With `fitness = TraitFitness(loci; σe, fitness)` the same holds for a quantitative trait: the genotypic value `G` sums `0`, `h*a` and `a` over the sites of `loci` (substituted sites included), the phenotype is `P = G + e` with `e ~ Normal(0, σe²)` drawn serially from `rng`, and fitness is `fitness(P)`, for example [`Stabilizing`](@ref)`(θ, ω²)` or [`Directional`](@ref)`(β)`.
- **Allocation Efficiency**: Working buffers (mutations, crossovers, mating assignments) are allocated once at start and reused, achieving zero steady-state per-generation allocations.

# Examples
```julia
using FisherWright

# Simulate 50 individuals for 10 generations with two 100 kb chromosomes
haps, cbp = fisher_wright(50, 10, [100_000, 100_000], 1.0)
length(haps) # 100 haplotypes

# Return a structured result tracking fixed mutations
res = fisher_wright(50, 10, [100_000, 100_000], 1.0; result=true)
res.chromosome_ends

# Continue that population for 5 generations of drift without new mutations
res2 = fisher_wright(50, 5, [100_000, 100_000], 1.0; init=res, mutate=false, result=true)

# A bottleneck of 10 individuals for 5 generations, then recovery to 50
sizes = [fill(10, 5); fill(50, 15)]
haps, _ = fisher_wright(50, 20, [100_000, 100_000], 1.0; popsize=sizes)
length(haps) # 100 haplotypes

# 5 sires and 45 dams in every generation: Ne ≈ 4 * 5 * 45 / 50 = 18
haps, _ = fisher_wright(50, 20, [100_000, 100_000], 1.0; mating=RandomMating(nmale=5))

# Exponential growth from 50, as a function of the generation
haps, _ = fisher_wright(50, 20, [100_000, 100_000], 1.0; popsize=g -> round(Int, 50 * 1.05^g))

# Selection at 100 sites with s = 0.05, h = 0.5, new alleles at rate 1e-4 per site
loci = fitness_loci([100_000, 100_000], 100; s=0.05, μ=1e-4)
res = fisher_wright(50, 100, [100_000, 100_000], 1.0; fitness=SiteFitness(loci), result=true)
res.selected_substitutions

# The same two trajectories from the helpers
haps, _ = fisher_wright(50, 20, [100_000, 100_000], 1.0; popsize=bottleneck(50, 10, 1, 5))
haps, _ = fisher_wright(50, 20, [100_000, 100_000], 1.0; popsize=growth(50, 0.05))
```
"""
function fisher_wright(
    ne::T1,
    nt::T2,
    chr::Vector{T3},
    mr::Float64;
    M=1e8,
    mut_base=1e8,
    result::Bool=false,
    fixation_interval::Int=1,
    mutate::Bool=true,
    init::Union{Nothing,Vector{Vector{UInt32}},FisherWrightResult}=nothing,
    popsize=nothing,
    mating::AbstractMating=RandomMating(),
    fitness::AbstractFitness=Neutral(),
    verbose::Bool=false,
    progress_interval::Integer=100,
    rng::AbstractRNG=Random.default_rng(),
) where {T1<:Integer,T2<:Integer,T3<:Integer}
    if !(ne > 1 && nt > 0 && all(chr .> 0) && 0.01 < mr < 20.0)
        throw(ArgumentError("Invalid parameter(s)"))
    end
    fixation_interval > 0 || throw(ArgumentError("fixation_interval must be positive"))
    mut_base > 0 || throw(ArgumentError("mut_base must be positive"))
    progress_interval > 0 || throw(ArgumentError("progress_interval must be positive"))

    tg = sum(chr)
    tg < 2^32 || error("Total genome length must be < 2^32 bp for UInt32 storage")
    cbp = UInt32.(cumsum(chr))
    tbp = cbp[end]

    p_mut = Poisson(tg / mut_base * mr)
    recomb_map = uniform_recombination_map(chr; M=M)
    span = UInt32(1):UInt32(tbp)

    demog = _demography(Int(ne), Int(nt), popsize)
    _check_mating(mating, demog, Int(ne), Int(nt))
    _check_fitness(fitness, tbp)
    _check_selection(mating, fitness)
    if _selected(fitness)
        result || throw(ArgumentError(
            "selection needs result = true: the selected alleles are returned in a FisherWrightResult"))
    end
    prt, substitutions = _initial_state(init, 2 * Int(ne), cbp, result)
    sel, ssubs = _initial_selected(init, fitness, 2 * Int(ne))
    st = SimState(prt, substitutions, _maxsize(demog, nt), rng, sel, ssubs)

    verbose &&
        @info "Fisher-Wright population simulation start" ne nt total_bp = Int(tg) threads =
            Threads.nthreads() progress_interval

    _simulate!(
        st, demog, mating, fitness, Int(nt), recomb_map, p_mut, span;
        mutate=mutate,
        track=result,
        fixation_interval=fixation_interval,
        verbose=verbose,
        progress_interval=Int(progress_interval),
        rng=rng,
    )
    # Drop the spare buffers beyond the final cohort.
    haps = resize!(st.prt, 2 * st.n)
    result || return (haps, cbp)
    census, males, fmean, flmean, flvar, gmean, gvar, pmean, pvar, ppairs, mpairs = _record(init, st)
    return FisherWrightResult(haps, cbp, st.substitutions, census, males,
        resize!(st.sel, 2 * st.n), st.ssubs, fmean, flmean, flvar, gmean, gvar, pmean, pvar, ppairs, mpairs)
end

"""
    fisher_wright(
        sp::Species,
        nt::Integer,
        mr::Float64 = 1.0;
        M = sp.M,
        mut_base = 1e8,
        result::Bool = false,
        fixation_interval::Int = 1,
        mutate::Bool = true,
        init = nothing,
        popsize = nothing,
        mating = RandomMating(),
        fitness = Neutral(),
        verbose::Bool = false,
        progress_interval::Integer = 100,
        rng::AbstractRNG = Random.default_rng(),
    ) -> Tuple{Vector{Vector{UInt32}}, Vector{UInt32}} | FisherWrightResult

Simulate a diploid Fisher-Wright population using species parameters from a
`BnGStructs.Species` instance (such as `Cattle`, `Pig`, `Chicken`, or `GenericSpecies`).

The population size `ne`, chromosome lengths `chr`, and default base pairs per Morgan `M`
are automatically extracted from `sp`:
- `ne = Int(sp.nid)`
- `chr = sp.chromosome`
- `M = sp.M`

# Arguments
- `sp::Species`: Species definition from `BnGStructs` specifying diploid population size (`nid`), chromosome lengths (`chromosome`), and recombination scale (`M`).
- `nt::Integer`: Number of generations to simulate (`nt > 0`).
- `mr::Float64 = 1.0`: Mutation rate per `mut_base` base pairs per meiosis (`0.01 < mr < 20.0`).

# Keywords
- `M = sp.M`: Number of base pairs per Morgan. Defaults to `sp.M`.
- `mut_base = 1e8`: Base pairs per unit of `mr`.
- `result::Bool = false`: If `true`, returns a [`FisherWrightResult`](@ref); otherwise returns `(prt, cbp)`.
- `fixation_interval::Int = 1`: Interval (in generations) to scan for and extract fixed substitutions when `result = true`.
- `mutate::Bool = true`: If `false`, no new mutations are generated; see the method above.
- `init = nothing`: Starting population of `2 * sp.nid` haplotypes, or a [`FisherWrightResult`](@ref); see the method above.
- `popsize = nothing`: Census size per generation, as a vector of `nt` integers or a function of the generation; see the method above.
- `mating = RandomMating()`: Mating scheme, for example `RandomMating(nmale = 10)`; see the method above.
- `fitness = Neutral()`: Fitness model, for example `SiteFitness(loci)` or `TraitFitness(loci; σe, fitness)`; requires `result = true`; see the method above.
- `verbose::Bool = false`: If `true`, reports progress to `stderr`; see the method above.
- `progress_interval::Integer = 100`: Generations between progress reports when `verbose = true`.
- `rng::AbstractRNG = Random.default_rng()`: Random number generator; see the method above.

# Examples
```julia
using BnGStructs, FisherWright

sp = GenericSpecies("Example", Int32(20), UInt32[100_000, 100_000], UInt32(50_000_000))
res = fisher_wright(sp, 10; result = true)
res.chromosome_ends
```
"""
function fisher_wright(
    sp::Species,
    nt::Integer,
    mr::Float64 = 1.0;
    M = sp.M,
    mut_base = 1e8,
    result::Bool = false,
    fixation_interval::Int = 1,
    mutate::Bool = true,
    init::Union{Nothing,Vector{Vector{UInt32}},FisherWrightResult} = nothing,
    popsize = nothing,
    mating::AbstractMating = RandomMating(),
    fitness::AbstractFitness = Neutral(),
    verbose::Bool = false,
    progress_interval::Integer = 100,
    rng::AbstractRNG = Random.default_rng(),
)
    return fisher_wright(
        Int(sp.nid),
        nt,
        sp.chromosome,
        mr;
        M = M,
        mut_base = mut_base,
        result = result,
        fixation_interval = fixation_interval,
        mutate = mutate,
        init = init,
        popsize = popsize,
        mating = mating,
        fitness = fitness,
        verbose = verbose,
        progress_interval = progress_interval,
        rng = rng,
    )
end
