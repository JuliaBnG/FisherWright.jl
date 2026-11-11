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
`_mate!(m, st, no, rng)`, which assigns sex to the `st.n` parents and writes
the sire and dam of each of `no` offspring to rows `1:no` of `st.pm`.
"""
abstract type AbstractMating end

"""
    AbstractFitness

Fitness of the adult cohort. A subtype implements `_fitness!(f, st)`, called
on the parents before mating.
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
- `srng`: one stream per meiosis/haplotype slot, seeded serially from the
  caller's rng, so a run does not depend on the thread count or scheduling.
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
    const srng::Vector{Xoshiro}
end

# State for the starting cohort `prt` (two haplotypes per individual), with
# buffers for up to `nmax` individuals; `prt` is extended in place with spare
# buffers. The only draws from `rng` are the 2nmax stream seeds, taken in
# slot order.
function SimState(
    prt::Vector{Vector{UInt32}},
    substitutions::Vector{UInt32},
    nmax::Int,
    rng::AbstractRNG,
)
    n = length(prt) ÷ 2
    (iseven(length(prt)) && 1 < n <= nmax) ||
        throw(ArgumentError("starting cohort of $(length(prt)) haplotypes does not fit nmax = $nmax"))
    nh = 2nmax
    bufs(k) = [Vector{UInt32}() for _ = 1:k]
    append!(prt, bufs(nh - 2n))
    return SimState(
        n, prt, bufs(nh), substitutions, bufs(nh), bufs(nh), bufs(nmax),
        Matrix{Int}(undef, nmax, 2), sizehint!(Vector{Bool}(undef, n), nmax), Int[], Int[],
        [Xoshiro(rand(rng, UInt64)) for _ = 1:nh],
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
    # A substitution is carried by every haplotype: kept apart when the run
    # tracks substitutions, merged back into the haplotypes when it does not.
    return result ? ([_remove_sorted!(copy(h), subs) for h in haps], copy(subs)) :
           ([merge_sorted(h, subs) for h in haps], UInt32[])
end

_fitness!(::Neutral, ::SimState) = nothing

# Conditioning on both sexes being present is exact rejection sampling;
# P(one sex only) = 2^(1-n), so a fixed retry cap fails for small n.
function _mate!(m::RandomMating, st::SimState, no::Int, rng::AbstractRNG)
    # Within the capacity reserved for the largest cohort: no allocation.
    resize!(st.sex, st.n)
    nm = _nmale(m, st.n)
    if nm == 0
        _random_mate!(st.pm, st.sex, st.sires, st.dams, no, rng, typemax(Int))
    else
        _fixed_sex_mate!(st.pm, st.sex, st.sires, st.dams, nm, no, rng)
    end
    return nothing
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
fitness on the adult cohort, mates, runs meiosis for the `_popsize(demog, g)`
offspring and then mutation; with
`track`, fixed sites are moved to `st.substitutions` every
`fixation_interval` generations and after the last one.

Mating is drawn serially from `rng`; meiosis and mutation use the per-slot
streams in `st.srng`.
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
    for g = 1:nt
        no = _popsize(demog, g)
        _fitness!(fitness, st)
        _mate!(mating, st, no, rng)
        _meiosis!(st, recomb_map, no)
        mutate && _mutation!(st, p_mut, span)
        if track && (g % fixation_interval == 0 || g == nt)
            st.substitutions = _fixation_step!(_cohort(st), st.substitutions)
        end
        if verbose && (g % progress_interval == 0 || g == nt)
            _report_progress(stderr, g, nt, t0, _cohort(st), length(st.substitutions), track)
        end
    end
    return st
end
