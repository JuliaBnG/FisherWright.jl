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

Census size per generation. A subtype implements `_popsize(d, g)`, the
number of offspring of generation `g`, and `_maxsize(d, nt)`, the largest
size over `nt` generations, by which the working buffers are allocated.
"""
abstract type AbstractDemography end

"""
    AbstractMating

Assignment of a sire and a dam to each offspring. A subtype implements
`_mate!(m, st, rng)`, which fills `st.pm`.
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

Constant census size of `n` diploid individuals: the only demography so far.
"""
struct ConstantSize <: AbstractDemography
    n::Int
end
_popsize(d::ConstantSize, ::Integer) = d.n
_maxsize(d::ConstantSize, ::Integer) = d.n

"""
    RandomMating()

Sex is redrawn each generation as Bernoulli(½) per individual, conditioned on
both sexes being present; each offspring draws its sire uniformly from the
males and its dam uniformly from the females, with replacement.
"""
struct RandomMating <: AbstractMating end

"""
    Neutral()

No selection: every parent has the same fitness, and nothing is evaluated.
"""
struct Neutral <: AbstractFitness end

"""
    SimState

Mutable state of one run. Individual `i` is haplotypes `2i - 1` and `2i`.

- `prt`: haplotypes of the current cohort; each sorted and unique.
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

# Buffers for up to `nmax` individuals. The only draws from `rng` are the
# 2nmax stream seeds, taken in slot order.
function SimState(
    prt::Vector{Vector{UInt32}},
    substitutions::Vector{UInt32},
    nmax::Int,
    rng::AbstractRNG,
)
    nh = 2nmax
    bufs(n) = [Vector{UInt32}() for _ = 1:n]
    return SimState(
        prt, bufs(nh), substitutions, bufs(nh), bufs(nh), bufs(nmax),
        Matrix{Int}(undef, nmax, 2), Vector{Bool}(undef, nmax), Int[], Int[],
        [Xoshiro(rand(rng, UInt64)) for _ = 1:nh],
    )
end

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
function _mate!(::RandomMating, st::SimState, rng::AbstractRNG)
    random_mate!(st.pm, st.sex, st.sires, st.dams; rng=rng, max_tries=typemax(Int))
    return nothing
end

# One meiosis per parent for each of `no` offspring (threaded), written into
# `st.off`, which then becomes the current cohort. Slot i draws its two
# crossover lists and strands from `st.srng[i]`, sire first.
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
    return nothing
end

# New mutations on each of the 2no haplotypes of the current cohort
# (threaded). Haplotype i draws from `st.srng[i]`.
function _mutation!(st::SimState, p_mut::Poisson{Float64}, span::UnitRange{UInt32}, no::Int)
    prt, mbuf, nbuf, srng, subs = st.prt, st.mbuf, st.nbuf, st.srng, st.substitutions
    Threads.@threads for i = 1:2no
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
fitness on the adult cohort, mates, runs meiosis and then mutation; with
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
        _mate!(mating, st, rng)
        _meiosis!(st, recomb_map, no)
        mutate && _mutation!(st, p_mut, span, no)
        if track && (g % fixation_interval == 0 || g == nt)
            st.substitutions = _fixation_step!(st.prt, st.substitutions)
        end
        if verbose && (g % progress_interval == 0 || g == nt)
            _report_progress(stderr, g, nt, t0, st.prt, length(st.substitutions), track)
        end
    end
    return st
end
