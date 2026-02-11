"""
    random_mate(n1::Integer, n2::Integer; rng=Random.default_rng())

Generate a random sex assignment for `n1` individuals (Bool: true = sire, false = dam),
then sample `n2` sire–dam pairs with replacement.

Returns (pm, sex) where:
- pm :: Matrix{Int} of size (n2, 2); columns: sire, dam
- sex :: Vector{Bool} length n1

Errors if a sex class is absent (resamples sex up to `max_tries` times).
"""
function random_mate(
    n1::Integer,
    n2::Integer;
    rng = Random.default_rng(),
    max_tries::Int = 10,
)
    (n1 > 1 && n2 > 0) || error("n1 must >1 and n2 >0")
    sex = rand(rng, Bool, n1)
    tries = 1
    while (all(sex) || all(.!sex)) && tries < max_tries
        sex .= rand(rng, Bool, n1)
        tries += 1
    end
    (all(sex) || all(.!sex)) && error("Failed to generate both sexes in $max_tries tries")
    # Build sire and dam index vectors once
    sir_idx = Vector{Int}()
    dam_idx = Vector{Int}()
    sizehint!(sir_idx, div(n1, 2)+1)
    sizehint!(dam_idx, div(n1, 2)+1)
    @inbounds for i = 1:n1
        if sex[i]
            push!(sir_idx, i)
        else
            push!(dam_idx, i)
        end
    end
    pm = Matrix{Int}(undef, n2, 2)
    @inbounds for i = 1:n2
        pm[i, 1] = rand(rng, sir_idx)
        pm[i, 2] = rand(rng, dam_idx)
    end
    return pm, sex
end

"""
    random_mate(sex::Vector, n2::Integer; rng=Random.default_rng())

Sample `n2` sire–dam pairs (with replacement) from an existing sex vector, where
sires are indicated by 1 and dams by 0. Returns pm :: Matrix{Int} (n2 × 2):
(sire, dam). An ID can be set as a value not equal to 0 or 1, e.g., -1, to
prevent it from being sampled as a sire or dam.
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
