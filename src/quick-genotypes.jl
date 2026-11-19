"""
    quickHap(
        nlc::Int,
        nid::Int;
        maf = 0.2,
        qd = Beta(0.75, 0.75),
        rng = Random.default_rng(),
        return_p::Bool = false,
    ) -> Haplotype | Tuple{Haplotype, Vector{Float64}}

Simulate synthetic haplotypes (0 or 1 allele calls) for `nid` diploid individuals
(`2 * nid` haplotypes) at `nlc` independent loci.

Allele frequencies ``p`` are sampled from `qd` and bounded by `maf < p < 1 - maf`.
Haploid alleles are then drawn from `Binomial(1, p)`.

# Arguments
- `nlc::Int`: Number of loci (matrix rows).
- `nid::Int`: Number of diploid individuals, yielding `2 * nid` haplotypes (matrix columns).

# Keywords
- `maf = 0.2`: Minor allele frequency threshold (`0 < maf < 0.5`).
- `qd = Beta(0.75, 0.75)`: Prior distribution used for sampling candidate allele frequencies.
- `rng = Random.default_rng()`: Random number generator for reproducibility.
- `return_p::Bool = false`: If `true`, returns `(haplotypes, freqs)` where `freqs` is the vector of accepted allele frequencies.

# Returns
- `Haplotype`: a `BnGStructs.Haplotype` of `nlc` loci and `2 * nid` haplotypes (individual `i` is haplotypes
  `2i - 1` and `2i`), `true` for the alternate allele. Before version 0.10 this was a `Matrix{Int8}` of 0 and 1;
  `Int8.(hp)` gives that matrix, with the same values for the same `rng`.
- (Optional) `Vector{Float64}`: Length-`nlc` vector of accepted allele frequencies when `return_p = true`.

# Errors
- Throws an `ArgumentError` if `maf <= 0` or `maf >= 0.5`.

# Examples
```julia
using FisherWright

hp = quickHap(100, 20)
size(hp) # (100, 40)
hp isa Haplotype # true
```
"""
function quickHap(
    nlc::Int,
    nid::Int;
    maf = 0.2,
    qd = Beta(0.75, 0.75),
    rng = Random.default_rng(),
    return_p::Bool = false,
)
    0 < maf < 0.5 || throw(ArgumentError("maf $maf not in (0, 0.5)"))
    nhp = 2 * nid
    hp = falses(nlc, nhp)
    freqs = return_p ? Vector{Float64}(undef, nlc) : nothing
    @inbounds for i = 1:nlc
        p = rand(rng, qd)
        while p ≤ maf || p ≥ 1 - maf
            p = rand(rng, qd)
        end
        freqs !== nothing && (freqs[i] = p)
        row = rand(rng, Binomial(1, p), nhp)
        @inbounds for j = 1:nhp
            hp[i, j] = row[j] == 1
        end
    end
    h = Haplotype(hp)
    return return_p ? (h, freqs) : h
end
