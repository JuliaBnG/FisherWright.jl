"""
    function recombine(h₁::Vector{UInt32}, h₂::Vector{UInt32}, hₒ::Vector{UInt32}, cross_overs::Vector{UInt32})
- Recombine two haplotypes into an offspring haplotype.

The recombination is performed by iterating through the crossover points and
alternating between the two haplotypes. The crossover points are used to
determine the segments of the haplotypes that are inherited by the offspring.
The function returns the recombined haplotype in the output vector `hₒ`.
"""
function recombine(
    h₁::Vector{UInt32},
    h₂::Vector{UInt32},
    hₒ::Vector{UInt32},
    cross_overs::Vector{UInt32},
)
    # Initialize the output vector
    empty!(hₒ)

    # Randomly select which haplotype to start with
    i, j, o = 1, 1, rand(Bool)
    m, n = length(h₁), length(h₂)

    for co in cross_overs
        while m > 0 && h₁[i] < co && i < m
            o && push!(hₒ, h₁[i])
            i += 1
        end
        while n > 0 && h₂[j] < co && j < n
            o || push!(hₒ, h₂[j])
            j += 1
        end
        o = !o # flip haplotype
    end
    o && i <= m && append!(hₒ, h₁[i:m])
    !o && j <= n && append!(hₒ, h₂[j:n])
    return hₒ
end

"""
    function cobp(cbp::Vector{UInt32}, pᵣ::Vector{Poisson{Float64}})
- Generate crossover points for recombination.

The crossover points are generated based on the chromosome breakpoints and the
Poisson distribution of recombination events. The crossover points are generated
in the range of [1, bp-1] for each chromosome plus the accumulated breakpoints.
The crossover points are sorted and returned as a vector of UInt32.
"""
function cobp(cbp::Vector{UInt32}, pᵣ::Vector{Poisson{Float64}})
    result, sbp = UInt32[], 1
    for (i, bp) in enumerate(cbp)
        nr = rand(pᵣ[i]) # number of recombination events on this chromosome
        append!(result, sort(rand(sbp:(bp-1), nr)))
        rand() < 0.5 && push!(result, bp) # if crossover is between chromosomes
        sbp = bp + 1
    end

    return result
end
