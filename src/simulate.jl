# Internal simulation engine: the state of a run and the generation loop.
#
# `fisher_wright` validates its arguments, builds a `SimState` and calls
# `_simulate!` once. `_simulate!` is specialised on the concrete demography,
# mating and fitness types, so the generation loop has no dynamic dispatch
# and the default models compile to the plain neutral random-mating loop
# (docs/dev/roadmap.md, design principle 2). Nothing here is public API yet:
# each axis has only its default model.

"""
    AbstractDemography

Census size per generation. Generation `g` uses the `N_{g-1}` individuals
of the previous cohort as parents and produces `N_g` offspring; `N_0` is the
starting size. A subtype implements `_popsize(d, g)`, which returns `N_g`
for `g ≥ 1`, and `_maxsize(d, nt)`, the largest of `N_0, …, N_nt`, by which
the working buffers are allocated.
"""
abstract type AbstractDemography end

"""
    AbstractMating

Assignment of a sire and a dam to each offspring. A subtype implements
`_mate!(m, f, st, no, rng)`, which assigns sex to the `st.n` parents and
writes the sire and dam of each of `no` offspring to rows `1:no` of `st.pm`,
drawn according to the fitness model `f`.
"""
abstract type AbstractMating end

"""
    AbstractFitness

Fitness of the adult cohort. A subtype implements `_fitness!(f, st)`, called
on the parents before mating, which writes the relative fitness of parent
`i` to `st.w[i]` for `i = 1:st.n`: finite, non-negative, and positive for at
least one parent of each sex. Parents are then drawn with probability
proportional to fitness within each sex (soft selection: the census size is
not affected). [`Neutral`](@ref) writes nothing and keeps the uniform draw.
"""
abstract type AbstractFitness end

"""
    ConstantSize(n)

Constant census size of `n` diploid individuals.
"""
struct ConstantSize <: AbstractDemography
    n::Int
end
_popsize(d::ConstantSize, ::Integer) = d.n
_maxsize(d::ConstantSize, ::Integer) = d.n

"""
    SizeTrajectory(n0, sizes)

Census size given in advance: `n0` individuals at the start and `sizes[g]`
offspring in generation `g`.
"""
struct SizeTrajectory <: AbstractDemography
    n0::Int
    sizes::Vector{Int}
end
_popsize(d::SizeTrajectory, g::Integer) = @inbounds d.sizes[g]
_maxsize(d::SizeTrajectory, nt::Integer) = max(d.n0, maximum(view(d.sizes, 1:Int(nt))))

# Demography from the `popsize` keyword of `fisher_wright`: nothing, a vector
# of `nt` sizes, or a callable `g -> N_g`, which is evaluated once for
# `g = 1:nt` so that the trajectory is known before the run.
_demography(ne::Int, ::Int, ::Nothing) = ConstantSize(ne)

function _demography(ne::Int, nt::Int, popsize::AbstractVector{<:Integer})
    length(popsize) == nt ||
        throw(ArgumentError("popsize has $(length(popsize)) entries; nt = $nt are required"))
    sizes = Vector{Int}(undef, nt)
    for (g, n) in enumerate(popsize)
        n > 1 || throw(ArgumentError("popsize must exceed 1 in every generation; it is $n in generation $g"))
        sizes[g] = n
    end
    return SizeTrajectory(ne, sizes)
end

function _demography(ne::Int, nt::Int, popsize)
    applicable(popsize, 1) || throw(ArgumentError(
        "popsize must be nothing, a vector of nt integers, or a function of the generation"))
    sizes = Vector{Int}(undef, nt)
    for g = 1:nt
        n = popsize(g)
        n isa Integer ||
            throw(ArgumentError("popsize($g) returned $(repr(n)); an integer is required"))
        sizes[g] = n
    end
    return _demography(ne, nt, sizes)
end

"""
    bottleneck(n, nb, start, duration; recovery = n) -> Function

Census trajectory `g -> N_g` for the `popsize` keyword of [`fisher_wright`](@ref):
`nb` individuals are born in each of the `duration` generations
`start, …, start + duration - 1`, `n` in every generation before them and
`recovery` in every generation after them.

Generation `g` draws its parents from the `N_{g-1}` individuals of the previous
generation, so the `nb` individuals of a bottleneck generation are the parents
of the next one: drift is increased in generations `start + 1, …, start + duration`.

All sizes must exceed 1, and `start` and `duration` must be positive. The
function can also be broadcast over generations to build or combine
trajectories as vectors, for example `bottleneck(100, 10, 21, 5).(1:60)`.

# Examples
```julia
using FisherWright

# 100 individuals, reduced to 10 in generations 21 to 25, then 100 again
res = fisher_wright(100, 60, [1_000_000], 1.0; popsize = bottleneck(100, 10, 21, 5), result = true)
res.census[20:28]   # N_19, …, N_27: 100, 100, 10, 10, 10, 10, 10, 100, 100
```

# See also
[`growth`](@ref)
"""
function bottleneck(n::Integer, nb::Integer, start::Integer, duration::Integer; recovery::Integer=n)
    (n > 1 && nb > 1 && recovery > 1) ||
        throw(ArgumentError("n, nb and recovery must exceed 1"))
    (start > 0 && duration > 0) ||
        throw(ArgumentError("start and duration must be positive"))
    n, nb, recovery = Int(n), Int(nb), Int(recovery)
    first, last = Int(start), Int(start) + Int(duration) - 1
    return g -> g < first ? n : g <= last ? nb : recovery
end

"""
    growth(n0, r; start = 0, limit = nothing) -> Function

Census trajectory `g -> N_g` for the `popsize` keyword of [`fisher_wright`](@ref):
geometric growth at rate `r` per generation from `n0` individuals,
`N_g = round(Int, n0 * (1 + r)^(g - start))` for `g > start` and `n0` up to
generation `start`. `r > -1`; a negative `r` is a decline.

With `limit`, the size stops there: it is the ceiling of a growing population
(`r > 0`) and the floor of a declining one (`r < 0`). Without it, a decline
that reaches a size below 2 is rejected by `fisher_wright`, which names the
generation, and growth is bounded only by memory: the run allocates for the
largest size of the trajectory.

`n0` should be the starting size `ne` of the run, so that the first generation
continues from it. The function can be broadcast over generations to obtain
the trajectory as a vector, `growth(50, 0.05).(1:40)`.

# Examples
```julia
using FisherWright

# 5% growth per generation from 50 individuals, up to 500
res = fisher_wright(50, 60, [1_000_000], 1.0; popsize = growth(50, 0.05; limit = 500), result = true)
res.census[[1, 2, 11, 61]]   # 50, 52, 81, 500

# Constant for 20 generations, then a 10% decline per generation down to 10
decline = growth(50, -0.1; start = 20, limit = 10)
```

# See also
[`bottleneck`](@ref)
"""
function growth(n0::Integer, r::Real; start::Integer=0, limit::Union{Nothing,Integer}=nothing)
    n0 > 1 || throw(ArgumentError("n0 must exceed 1"))
    (isfinite(r) && r > -1) || throw(ArgumentError("the growth rate r must be finite and exceed -1"))
    start >= 0 || throw(ArgumentError("start must not be negative"))
    limit === nothing || limit > 1 || throw(ArgumentError("limit must exceed 1"))
    if limit !== nothing && r != 0 && (r > 0) != (limit > n0)
        throw(ArgumentError("limit = $limit is on the wrong side of n0 = $n0 for r = $r"))
    end
    n0, start, factor = Int(n0), Int(start), 1 + Float64(r)
    bound = limit === nothing ? nothing : Float64(limit)
    return function (g::Integer)
        g <= start && return n0
        x = n0 * factor^(g - start)
        if bound !== nothing
            x = factor > 1 ? min(x, bound) : max(x, bound)
        end
        x < 2.0^62 || throw(ArgumentError("growth gives a census size beyond Int in generation $g"))
        return round(Int, x)
    end
end

"""
    RandomMating(; nmale = nothing, sexratio = nothing)

Random mating between the sexes: each offspring draws its sire uniformly from
the males and its dam uniformly from the females of the parental cohort, with
replacement. Sex is assigned anew to every cohort when it reproduces, by one
of three models:

- `RandomMating()`: Bernoulli(½) per individual, conditioned on both sexes
  being present, so the number of males among `N` parents is Binomial(`N`, ½)
  truncated to `1:N-1`. This is the default of [`fisher_wright`](@ref).
- `RandomMating(nmale = m)`: exactly `m` males and `N - m` females, the males
  being a uniformly random subset of the cohort. Every parental cohort must
  have more than `m` individuals.
- `RandomMating(sexratio = p)`: exactly `round(Int, p * N)` males, limited to
  `1:N-1`, for a proportion of males `0 < p < 1`; the number of males then
  follows a changing census size.

With `N_m` males and `N_f` females the effective size of that generation is
approximately `4 N_m N_f / (N_m + N_f)` (Wright 1931).
"""
struct RandomMating <: AbstractMating
    nmale::Int          # fixed number of males; 0 if not fixed
    sexratio::Float64   # fixed proportion of males; 0 if not fixed
    function RandomMating(; nmale::Union{Nothing,Integer}=nothing, sexratio::Union{Nothing,Real}=nothing)
        nmale === nothing || sexratio === nothing ||
            throw(ArgumentError("give nmale or sexratio, not both"))
        nmale === nothing || nmale > 0 || throw(ArgumentError("nmale must be positive"))
        sexratio === nothing || 0 < sexratio < 1 ||
            throw(ArgumentError("sexratio is the proportion of males and must lie strictly between 0 and 1"))
        return new(something(nmale, 0), something(sexratio, 0.0))
    end
end

# Number of males among `n` parents; 0 stands for binomial sex.
function _nmale(m::RandomMating, n::Int)
    m.nmale > 0 && return m.nmale
    m.sexratio > 0 && return clamp(round(Int, m.sexratio * n), 1, n - 1)
    return 0
end

# Every parental cohort (generations 0 to nt - 1) must be able to hold the
# requested sexes; checked before the run so that it cannot fail midway.
function _check_mating(m::RandomMating, demog::AbstractDemography, n0::Int, nt::Int)
    m.nmale > 0 || return nothing
    n = n0
    for g = 1:nt
        m.nmale < n || throw(ArgumentError(
            "nmale = $(m.nmale) leaves no female among the $n parents of generation $g"))
        n = _popsize(demog, g)
    end
    return nothing
end

"""
    Neutral()

No selection: every parent has the same fitness, and nothing is evaluated.
"""
struct Neutral <: AbstractFitness end

"""
    SimState

Mutable state of one run. Individual `i` is haplotypes `2i - 1` and `2i`.
Every buffer is allocated for the largest census size of the run, and the
current cohort occupies the leading entries.

- `n`: number of individuals in the current cohort.
- `prt`: haplotypes; the first `2n` are the current cohort (`_cohort`), each
  sorted and unique, and the rest are spare buffers.
- `off`: target of the next meiosis; swapped with `prt` every generation.
- `substitutions`: sorted positions fixed and removed from `prt`.
- `mbuf`, `nbuf`: per-haplotype scratch for mutation insertion and for the
  new mutations of one generation.
- `cbuf`: per-offspring crossover buffer.
- `pm`, `sex`, `sires`, `dams`: mating table (sire and dam of each
  offspring) and its scratch.
- `w`: relative fitness of each parent, written by `_fitness!`; unused under
  `Neutral`.
- `tm`, `tf`: alias tables for the weighted draw of sires and dams.
- `srng`: one stream per meiosis/haplotype slot, seeded serially from the
  caller's rng, so a run does not depend on the thread count or scheduling.
- `sel`, `soff`, `ssubs`, `smbuf`, `snbuf`: the selected stream, laid out as
  `prt`, `off`, `substitutions`, `mbuf` and `nbuf`: the derived alleles each
  haplotype carries at the fitness loci. Empty unless the fitness model
  keeps selected alleles.
- `fmean`, `flmean`, `flvar`: per generation, the arithmetic mean of the
  fitness of the parents and the mean and variance of its logarithm.
- `fcnt`, `mcnt`, `fidx`, `bsires`, `bdams`: scratch of mating schemes with
  family sizes (`bsires`, `bdams`: the breeding sires and dams).
- `gval`, `pval`: genotypic value and phenotype of each parent under a
  trait model; `gmean`, `gvar`, `pmean`, `pvar`: their mean and variance per
  generation (NaN without a trait).
- `pcnt`: scratch for the number of offspring of each parent.
- `swt`: scratch for the split weights of parents under selection.
- `ppairs`, `mpairs`: per generation, the number of ordered pairs of
  offspring that share their sire and that share their dam.
- `census`, `males`: record of the run. `census` starts with the size of the
  starting cohort and gains the number of offspring of every generation;
  `males` gains the number of males among its parents.
"""
mutable struct SimState
    n::Int
    prt::Vector{Vector{UInt32}}
    off::Vector{Vector{UInt32}}
    substitutions::Vector{UInt32}
    const mbuf::Vector{Vector{UInt32}}
    const nbuf::Vector{Vector{UInt32}}
    const cbuf::Vector{Vector{UInt32}}
    const pm::Matrix{Int}
    const sex::Vector{Bool}
    const sires::Vector{Int}
    const dams::Vector{Int}
    const w::Vector{Float64}
    const tm::AliasTable
    const tf::AliasTable
    const srng::Vector{Xoshiro}
    const census::Vector{Int}
    const males::Vector{Int}
    sel::Vector{Vector{UInt32}}
    soff::Vector{Vector{UInt32}}
    ssubs::Vector{UInt32}
    const smbuf::Vector{Vector{UInt32}}
    const snbuf::Vector{Vector{UInt32}}
    const fmean::Vector{Float64}
    const flmean::Vector{Float64}
    const flvar::Vector{Float64}
    const gval::Vector{Float64}
    const pval::Vector{Float64}
    const gmean::Vector{Float64}
    const gvar::Vector{Float64}
    const pmean::Vector{Float64}
    const pvar::Vector{Float64}
    const fcnt::Vector{Int}
    const mcnt::Vector{Int}
    const fidx::Vector{Int}
    const bsires::Vector{Int}
    const bdams::Vector{Int}
    const pcnt::Vector{Int}
    const swt::Vector{Float64}
    const ppairs::Vector{Float64}
    const mpairs::Vector{Float64}
end

# State for the starting cohort `prt` (two haplotypes per individual), with
# buffers for up to `nmax` individuals; `prt` is extended in place with spare
# buffers, as is `sel`, the selected stream of the same cohort (none by
# default). The only draws from `rng` are the 2nmax stream seeds, taken in
# slot order.
function SimState(
    prt::Vector{Vector{UInt32}},
    substitutions::Vector{UInt32},
    nmax::Int,
    rng::AbstractRNG,
    sel::Vector{Vector{UInt32}}=[Vector{UInt32}() for _ in prt],
    ssubs::Vector{UInt32}=UInt32[],
)
    n = length(prt) ÷ 2
    (iseven(length(prt)) && 1 < n <= nmax) ||
        throw(ArgumentError("starting cohort of $(length(prt)) haplotypes does not fit nmax = $nmax"))
    length(sel) == length(prt) ||
        throw(ArgumentError("the selected stream must have one entry per haplotype"))
    nh = 2nmax
    bufs(k) = [Vector{UInt32}() for _ = 1:k]
    append!(prt, bufs(nh - 2n))
    append!(sel, bufs(nh - 2n))
    return SimState(
        n, prt, bufs(nh), substitutions, bufs(nh), bufs(nh), bufs(nmax),
        Matrix{Int}(undef, nmax, 2), sizehint!(Vector{Bool}(undef, n), nmax), Int[], Int[],
        ones(nmax), AliasTable(nmax), AliasTable(nmax),
        [Xoshiro(rand(rng, UInt64)) for _ = 1:nh], [n], Int[],
        sel, bufs(nh), ssubs, bufs(nh), bufs(nh),
        Float64[], Float64[], Float64[],
        zeros(nmax), zeros(nmax), Float64[], Float64[], Float64[], Float64[],
        sizehint!(Int[], nmax), sizehint!(Int[], nmax), sizehint!(Int[], nmax),
        sizehint!(Int[], nmax), sizehint!(Int[], nmax),
        zeros(Int, nmax), zeros(nmax), Float64[], Float64[],
    )
end

# Haplotypes of the current cohort.
_cohort(st::SimState) = view(st.prt, 1:2*st.n)

# Positions of a starting population: each vector strictly increasing and
# within the genome.
function _check_positions(v::Vector{UInt32}, tbp::UInt32, what)
    isempty(v) && return nothing
    @inbounds begin
        v[1] >= 1 || throw(ArgumentError("$what: positions are 1-based"))
        for k = 2:length(v)
            v[k-1] < v[k] ||
                throw(ArgumentError("$what: positions must be sorted and unique"))
        end
        v[end] <= tbp ||
            throw(ArgumentError("$what: position $(Int(v[end])) exceeds the genome length $(Int(tbp))"))
    end
    return nothing
end

function _check_haplotypes(haps::Vector{Vector{UInt32}}, nh::Int, tbp::UInt32)
    length(haps) == nh ||
        throw(ArgumentError("init has $(length(haps)) haplotypes; 2ne = $nh are required"))
    for (i, h) in pairs(haps)
        _check_positions(h, tbp, "init haplotype $i")
    end
    return nothing
end

# Starting haplotypes and substitutions. `init` is copied, never modified or
# aliased by the run.
_initial_state(::Nothing, nh::Int, cbp::Vector{UInt32}, result::Bool) =
    ([Vector{UInt32}() for _ = 1:nh], UInt32[])

function _initial_state(init::Vector{Vector{UInt32}}, nh::Int, cbp::Vector{UInt32}, result::Bool)
    _check_haplotypes(init, nh, cbp[end])
    return ([copy(h) for h in init], UInt32[])
end

function _initial_state(init::FisherWrightResult, nh::Int, cbp::Vector{UInt32}, result::Bool)
    init.chromosome_ends == cbp ||
        throw(ArgumentError("init.chromosome_ends does not match cumsum(chr)"))
    haps, subs = init.active_haplotypes, init.substitutions
    _check_haplotypes(haps, nh, cbp[end])
    _check_positions(subs, cbp[end], "init.substitutions")
    (length(init.census) == length(init.males) + 1 && init.census[end] == nh ÷ 2) ||
        throw(ArgumentError("init.census and init.males do not describe the $(nh ÷ 2) individuals of init"))
    all(v -> length(v) == length(init.males), (init.mean_fitness, init.mean_log_fitness,
        init.var_log_fitness, init.mean_genotypic, init.var_genotypic, init.mean_phenotypic,
        init.var_phenotypic, init.paternal_pairs, init.maternal_pairs)) ||
        throw(ArgumentError("the fitness and trait records of init must have one entry per generation"))
    # A substitution is carried by every haplotype: kept apart when the run
    # tracks substitutions, merged back into the haplotypes when it does not.
    return result ? ([_remove_sorted!(copy(h), subs) for h in haps], copy(subs)) :
           ([merge_sorted(h, subs) for h in haps], UInt32[])
end

# Census record of the result: that of this call, appended to the record of
# a FisherWrightResult the run was continued from.
_record(::Any, st::SimState) =
    (st.census, st.males, st.fmean, st.flmean, st.flvar, st.gmean, st.gvar, st.pmean, st.pvar,
        st.ppairs, st.mpairs)
_record(init::FisherWrightResult, st::SimState) = (
    vcat(init.census, view(st.census, 2:length(st.census))), vcat(init.males, st.males),
    vcat(init.mean_fitness, st.fmean), vcat(init.mean_log_fitness, st.flmean),
    vcat(init.var_log_fitness, st.flvar), vcat(init.mean_genotypic, st.gmean),
    vcat(init.var_genotypic, st.gvar), vcat(init.mean_phenotypic, st.pmean),
    vcat(init.var_phenotypic, st.pvar), vcat(init.paternal_pairs, st.ppairs),
    vcat(init.maternal_pairs, st.mpairs))

# Realized family sizes of one generation: the number of ordered pairs of
# its `no` offspring that share their sire, and that share their dam, from
# the mating table.
function _record_families!(st::SimState, no::Int)
    c, pm = st.pcnt, st.pm
    for (col, out) in ((1, st.ppairs), (2, st.mpairs))
        fill!(view(c, 1:st.n), 0)
        @inbounds for i = 1:no
            c[pm[i, col]] += 1
        end
        pairs = 0
        @inbounds for k = 1:st.n
            pairs += c[k] * (c[k] - 1)
        end
        push!(out, pairs)
    end
    return nothing
end

_fitness!(::Neutral, ::SimState) = nothing

# Fitness record of the parents of one generation. A model that evaluates
# fitness on a scale of its own records it in `_fitness!` and makes this a
# no-op; otherwise the generation is recorded as neutral.
function _record_fitness!(::AbstractFitness, st::SimState)
    push!(st.fmean, 1.0)
    push!(st.flmean, 0.0)
    push!(st.flvar, 0.0)
    return nothing
end

# Trait record of the parents of one generation: NaN unless the fitness
# model is a trait model, which records it in `_fitness!`.
function _record_trait!(::AbstractFitness, st::SimState)
    push!(st.gmean, NaN)
    push!(st.gvar, NaN)
    push!(st.pmean, NaN)
    push!(st.pvar, NaN)
    return nothing
end

# Fitness models that draw random numbers (environmental deviations) take
# the run's serial rng; the others ignore it.
_fitness!(f::AbstractFitness, st::SimState, ::AbstractRNG) = _fitness!(f, st)

# Errors in a fitness model that the run would otherwise hit midway.
_check_fitness(::AbstractFitness, ::UInt32) = nothing

# Combinations of mating scheme and fitness model that are not supported.
_check_selection(::AbstractMating, ::AbstractFitness) = nothing

# Sex of the parents, then sire and dam of each of `no` offspring.
# Conditioning on both sexes being present is exact rejection sampling;
# P(one sex only) = 2^(1-n), so a fixed retry cap fails for small n.
function _mate!(m::RandomMating, f::AbstractFitness, st::SimState, no::Int, rng::AbstractRNG)
    # Within the capacity reserved for the largest cohort: no allocation.
    resize!(st.sex, st.n)
    nm = _nmale(m, st.n)
    if nm == 0
        _draw_sex!(st.sex, st.sires, st.dams, rng, typemax(Int))
    else
        _draw_fixed_sex!(st.sex, st.sires, st.dams, nm, rng)
    end
    _parents!(f, st, no, rng)
    return nothing
end

# Uniform parents under neutrality; otherwise proportional to fitness within
# each sex.
_parents!(::Neutral, st::SimState, no::Int, rng::AbstractRNG) =
    _draw_parents!(st.pm, st.sires, st.dams, no, rng)
function _parents!(::AbstractFitness, st::SimState, no::Int, rng::AbstractRNG)
    w = st.w
    # Soft selection cannot rescue a sex without any fit member.
    (any(i -> w[i] > 0, st.sires) && any(i -> w[i] > 0, st.dams)) || error(
        "no parent with positive fitness in one sex: the population cannot reproduce")
    return _draw_parents!(st.pm, st.sires, st.dams, no, rng, w, st.tm, st.tf)
end

# One meiosis per parent for each of `no` offspring (threaded), written into
# `st.off`, which then becomes the current cohort of `no` individuals. Slot i
# draws its two crossover lists and strands from `st.srng[i]`, sire first.
function _meiosis!(st::SimState, recomb_map::RecombinationMap, no::Int)
    prt, off, pm, cbuf, srng = st.prt, st.off, st.pm, st.cbuf, st.srng
    Threads.@threads for i = 1:no
        s = pm[i, 1]
        d = pm[i, 2]
        # recombine empties its target, and cobp! its crossover buffer, so
        # both are safe to reuse across generations.
        co = cbuf[i]
        r = srng[i]
        recombine(prt[2s-1], prt[2s], off[2i-1], cobp!(co, recomb_map; rng=r); rng=r)
        recombine(prt[2d-1], prt[2d], off[2i], cobp!(co, recomb_map; rng=r); rng=r)
    end
    st.prt, st.off = off, prt
    st.n = no
    return nothing
end

# As `_meiosis!`, for a fitness model with selected alleles: each meiosis
# applies one crossover list and one starting strand to the neutral and the
# selected haplotypes of the parent. The draws from `st.srng[i]` are those
# of `_meiosis!`, in the same order.
function _meiosis_selected!(st::SimState, recomb_map::RecombinationMap, no::Int)
    prt, off, sel, soff = st.prt, st.off, st.sel, st.soff
    pm, cbuf, srng = st.pm, st.cbuf, st.srng
    Threads.@threads for i = 1:no
        co = cbuf[i]
        r = srng[i]
        for (p, k) in ((pm[i, 1], 2i - 1), (pm[i, 2], 2i))
            cobp!(co, recomb_map; rng=r)
            o = rand(r, Bool)
            _recombine!(prt[2p-1], prt[2p], off[k], co, o)
            _recombine!(sel[2p-1], sel[2p], soff[k], co, o)
        end
    end
    st.prt, st.off = off, prt
    st.sel, st.soff = soff, sel
    st.n = no
    return nothing
end

# New mutations on each haplotype of the current cohort (threaded).
# Haplotype i draws from `st.srng[i]`.
function _mutation!(st::SimState, p_mut::Poisson{Float64}, span::UnitRange{UInt32})
    prt, mbuf, nbuf, srng, subs = st.prt, st.mbuf, st.nbuf, st.srng, st.substitutions
    Threads.@threads for i = 1:2*st.n
        r = srng[i]
        nm = rand(r, p_mut)
        if nm > 0
            newm = nbuf[i]
            resize!(newm, nm)
            @inbounds for k = 1:nm
                newm[k] = rand(r, span)
            end
            nm <= 32 ? sort!(newm; alg=InsertionSort) : sort!(newm)
            # A substituted site is carried by every haplotype, so a new
            # mutation there is a no-op, as it is when result = false.
            _drop_present!(newm, subs)
            # Insert into the scratch buffer by block copies, then swap it
            # in: no allocation once both buffers have reached their
            # steady-state size.
            if !isempty(newm)
                _insert_sorted!(mbuf[i], prt[i], newm)
                prt[i], mbuf[i] = mbuf[i], prt[i]
            end
        end
    end
    return nothing
end

"""
    _simulate!(st, demog, mating, fitness, nt, recomb_map, p_mut, span; kwargs...) -> SimState

Advance `st` by `nt` generations and return it. Each generation evaluates
fitness on the adult cohort, mates (parents drawn in proportion to fitness
within each sex unless `fitness` is `Neutral`), runs meiosis for the `_popsize(demog, g)`
offspring and then mutation; with
`track`, fixed sites are moved to `st.substitutions` every
`fixation_interval` generations and after the last one.

Mating is drawn serially from `rng`; meiosis and mutation use the per-slot
streams in `st.srng`. The census size and the number of male parents of
every generation are appended to `st.census` and `st.males`.

With a fitness model that keeps selected alleles (`_selected(fitness)`),
meiosis carries the selected stream `st.sel` along, mutation also hits the
fitness loci, and fixed selected sites are moved to `st.ssubs`.
"""
function _simulate!(
    st::SimState,
    demog::AbstractDemography,
    mating::AbstractMating,
    fitness::AbstractFitness,
    nt::Int,
    recomb_map::RecombinationMap,
    p_mut::Poisson{Float64},
    span::UnitRange{UInt32};
    mutate::Bool,
    track::Bool,
    fixation_interval::Int,
    verbose::Bool,
    progress_interval::Int,
    rng::AbstractRNG,
)
    t0 = time()
    selected = _selected(fitness)
    p_sel = _sel_poisson(fitness)
    for g = 1:nt
        no = _popsize(demog, g)
        _fitness!(fitness, st, rng)
        _record_fitness!(fitness, st)
        _record_trait!(fitness, st)
        _mate!(mating, fitness, st, no, rng)
        push!(st.males, length(st.sires))
        _record_families!(st, no)
        selected ? _meiosis_selected!(st, recomb_map, no) : _meiosis!(st, recomb_map, no)
        push!(st.census, no)
        mutate && _mutation!(st, p_mut, span)
        mutate && p_sel !== nothing && _sel_mutation!(fitness, st, p_sel)
        if track && (g % fixation_interval == 0 || g == nt)
            st.substitutions = _fixation_step!(_cohort(st), st.substitutions)
            # A site fixed in the selected stream multiplies every fitness
            # by the same factor, so removing it leaves selection unchanged.
            selected && (st.ssubs = _fixation_step!(view(st.sel, 1:2*st.n), st.ssubs))
        end
        if verbose && (g % progress_interval == 0 || g == nt)
            _report_progress(stderr, g, nt, t0, _cohort(st), length(st.substitutions), track)
        end
    end
    return st
end
