"""
    random_mate(
        n1::Integer,
        n2::Integer;
        rng::AbstractRNG = Random.default_rng(),
        max_tries::Int = 10,
    ) -> Tuple{Matrix{Int}, Vector{Bool}}

Randomly assign biological sexes to `n1` individuals and sample `n2` sire–dam mating pairs
with replacement.

# Arguments
- `n1::Integer`: Number of parental individuals (`n1 > 1`).
- `n2::Integer`: Number of offspring mating pairs to sample (`n2 > 0`).

# Keywords
- `rng::AbstractRNG = Random.default_rng()`: Random number generator.
- `max_tries::Int = 10`: Maximum attempts to generate a valid sex assignment where both sexes are present.

# Returns
- `pm::Matrix{Int}`: An `n2 × 2` matrix of 1-based individual IDs where column 1 contains sires and column 2 contains dams.
- `sex::Vector{Bool}`: A boolean vector of length `n1` indicating individual sexes (`true` = sire / male, `false` = dam / female).

# Errors
- Throws an `ErrorException` if `n1 <= 1` or `n2 <= 0`.
- Throws an `ErrorException` if both sexes cannot be sampled within `max_tries` iterations.

# See also
[`random_mate!`](@ref)

# Examples
```julia
using FisherWright

pm, sex = random_mate(20, 10)
size(pm) == (10, 2)
length(sex) == 20
```
"""
function random_mate(
    n1::Integer,
    n2::Integer;
    rng = Random.default_rng(),
    max_tries::Int = 10,
)
    (n1 > 1 && n2 > 0) || error("n1 must >1 and n2 >0")
    return random_mate!(
        Matrix{Int}(undef, n2, 2),
        Vector{Bool}(undef, n1),
        Int[],
        Int[];
        rng = rng,
        max_tries = max_tries,
    )
end
"""
    random_mate!(
        pm::Matrix{Int},
        sex::Vector{Bool},
        sires::Vector{Int},
        dams::Vector{Int};
        rng::AbstractRNG = Random.default_rng(),
        max_tries::Int = 10,
    ) -> Tuple{Matrix{Int}, Vector{Bool}}

In-place form of [`random_mate`](@ref), writing sex assignments into `sex` and sampled mating pairs
into `pm`, using `sires` and `dams` as pre-allocated index buffers.

Reusing these buffers across simulation generations avoids repeated heap allocations.

# Arguments
- `pm::Matrix{Int}`: Destination mating matrix of size `(n2, 2)`.
- `sex::Vector{Bool}`: Destination sex vector of length `n1` (`true` = sire, `false` = dam).
- `sires::Vector{Int}`: Scratch buffer for sire indices.
- `dams::Vector{Int}`: Scratch buffer for dam indices.

# Keywords
- `rng::AbstractRNG = Random.default_rng()`: Random number generator.
- `max_tries::Int = 10`: Maximum retries allowed to generate a population with both sexes present.

# Returns
- `(pm, sex)`: Tuple referencing the populated destination buffers.

# Errors
- Throws an `ErrorException` if `length(sex) <= 1` or `size(pm, 1) <= 0`.
- Throws an `ErrorException` if both sexes fail to appear within `max_tries` draws.
"""
function random_mate!(
    pm::Matrix{Int},
    sex::Vector{Bool},
    sires::Vector{Int},
    dams::Vector{Int};
    rng = Random.default_rng(),
    max_tries::Int = 10,
)
    return _random_mate!(pm, sex, sires, dams, size(pm, 1), rng, max_tries)
end

# Kernel of `random_mate!`: `length(sex)` parents, and matings written to rows
# `1:n2` of `pm`, which may have more rows. The draws from `rng` are the sex
# vector (redrawn until both sexes are present), then sire and dam of each
# offspring in turn.
function _random_mate!(
    pm::Matrix{Int},
    sex::Vector{Bool},
    sires::Vector{Int},
    dams::Vector{Int},
    n2::Int,
    rng,
    max_tries::Int,
)
    n1 = length(sex)
    (n1 > 1 && 0 < n2 <= size(pm, 1)) || error("n1 must >1 and n2 >0")
    _draw_sex!(sex, sires, dams, rng, max_tries)
    _draw_parents!(pm, sires, dams, n2, rng)
    return pm, sex
end

# Sex of each parent as Bernoulli(½), redrawn until both sexes are present;
# `sires` and `dams` receive the sorted indices of the males and females.
function _draw_sex!(sex::Vector{Bool}, sires::Vector{Int}, dams::Vector{Int}, rng, max_tries::Int)
    n1 = length(sex)
    tries = 0
    while true
        rand!(rng, sex)
        tries += 1
        (any(sex) && !all(sex)) && break
        tries >= max_tries &&
            error("Failed to generate both sexes in $max_tries tries")
    end
    # Build sire and dam index vectors once
    empty!(sires)
    empty!(dams)
    sizehint!(sires, div(n1, 2)+1)
    sizehint!(dams, div(n1, 2)+1)
    @inbounds for i = 1:n1
        if sex[i]
            push!(sires, i)
        else
            push!(dams, i)
        end
    end
    return sex
end

# Sire and dam of each of `n2` offspring, uniform with replacement.
function _draw_parents!(pm::Matrix{Int}, sires::Vector{Int}, dams::Vector{Int}, n2::Int, rng)
    @inbounds for i = 1:n2
        pm[i, 1] = rand(rng, sires)
        pm[i, 2] = rand(rng, dams)
    end
    return pm
end

# Sire and dam of each of `n2` offspring, with replacement and with
# probability proportional to the weight `w[i]` of parent `i` within its sex.
# `tm` and `tf` are rebuilt for the current sires and dams. Per offspring the
# draws from `rng` are two for the sire, then two for the dam.
function _draw_parents!(
    pm::Matrix{Int},
    sires::Vector{Int},
    dams::Vector{Int},
    n2::Int,
    rng,
    w::AbstractVector{Float64},
    tm::AliasTable,
    tf::AliasTable,
)
    _alias_build!(tm, w, sires)
    _alias_build!(tf, w, dams)
    @inbounds for i = 1:n2
        pm[i, 1] = sires[_alias_draw(tm, rng)]
        pm[i, 2] = dams[_alias_draw(tf, rng)]
    end
    return pm
end

# As `_random_mate!`, with exactly `nm` males among the `length(sex)` parents.
function _fixed_sex_mate!(
    pm::Matrix{Int},
    sex::Vector{Bool},
    sires::Vector{Int},
    dams::Vector{Int},
    nm::Int,
    n2::Int,
    rng,
)
    n1 = length(sex)
    (0 < nm < n1 && 0 < n2 <= size(pm, 1)) ||
        error("need 0 < nm < n1 and n2 > 0")
    _draw_fixed_sex!(sex, sires, dams, nm, rng)
    _draw_parents!(pm, sires, dams, n2, rng)
    return pm, sex
end

# Exactly `nm` males among the `length(sex)` parents, a uniformly random
# subset drawn by selection sampling (Knuth, TAOCP 2, Algorithm S): one
# uniform draw per parent that is still undecided.
function _draw_fixed_sex!(sex::Vector{Bool}, sires::Vector{Int}, dams::Vector{Int}, nm::Int, rng)
    n1 = length(sex)
    empty!(sires)
    empty!(dams)
    need = nm
    @inbounds for i = 1:n1
        left = n1 - i + 1
        # Male with probability need / left; no draw once it is 0 or 1.
        male = need > 0 && (need == left || rand(rng) * left < need)
        sex[i] = male
        if male
            need -= 1
            push!(sires, i)
        else
            push!(dams, i)
        end
    end
    return sex
end

"""
    random_mate(
        sex::Vector,
        n2::Integer;
        rng::AbstractRNG = Random.default_rng(),
    ) -> Matrix{Int}

Sample `n2` sire–dam mating pairs (with replacement) from an existing vector of sexes.

# Arguments
- `sex::Vector`: Vector indicating individual sexes, where `1` designates a sire and `0` designates a dam.
  Individuals with any other value (e.g. `-1` or `missing`) are excluded from mating selection.
- `n2::Integer`: Number of mating pairs to sample (`n2 > 0`).

# Keywords
- `rng::AbstractRNG = Random.default_rng()`: Random number generator.

# Returns
- `Matrix{Int}`: Sorted `n2 × 2` matrix of mating pairs (column 1: sire IDs, column 2: dam IDs).

# Errors
- Throws an `ErrorException` if `sex` contains no sires (`1`) or no dams (`0`).

# Examples
```julia
using FisherWright

sex = [1, 1, 0, 0, -1] # 2 sires, 2 dams, 1 excluded individual
pm = random_mate(sex, 5)
size(pm) == (5, 2)
```
"""
function random_mate(sex::Vector, n2::Integer; rng = Random.default_rng())
    sir_idx = findall(==(1), sex)
    dam_idx = findall(==(0), sex)
    length(sir_idx) > 0 || error("No sires in sex vector")
    length(dam_idx) > 0 || error("No dams in sex vector")
    pm = Matrix{Int}(undef, n2, 2)
    @inbounds for i = 1:n2
        pm[i, 1] = rand(rng, sir_idx)
        pm[i, 2] = rand(rng, dam_idx)
    end
    return sortslices(pm, dims=1)
end
