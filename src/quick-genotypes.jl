"""
    quickGT(
        nlc::Int,
        nid::Int;
        maf = 0.1,
        qd = Beta(0.75, 0.75),
        rng = Random.default_rng(),
        return_p::Bool = false,
    ) -> Matrix{Int8} | Tuple{Matrix{Int8}, Vector{Float64}}

Simulate synthetic diploid genotypes (0, 1, or 2 alternate allele copies) for `nid` individuals
at `nlc` independent biallelic loci.

Allele frequencies ``p`` are drawn from the prior distribution `qd` and accepted only when
`maf < p < 1 - maf`. Diploid genotypes for each locus are then drawn from `Binomial(2, p)`.

# Arguments
- `nlc::Int`: Number of loci (matrix rows).
- `nid::Int`: Number of diploid individuals (matrix columns).

# Keywords
- `maf = 0.1`: Minor allele frequency threshold (`0 < maf < 0.5`). Frequencies outside `(maf, 1 - maf)` are rejected.
- `qd = Beta(0.75, 0.75)`: Prior distribution used for sampling candidate allele frequencies.
- `rng = Random.default_rng()`: Random number generator for reproducibility.
- `return_p::Bool = false`: If `true`, returns a tuple `(genotypes, freqs)` containing the simulated genotypes
  and the vector of accepted allele frequencies `freqs`. If `false`, returns only `genotypes`.

# Returns
- `Matrix{Int8}`: Genotype matrix of shape `(nlc, nid)` with values in `{0, 1, 2}`.
- (Optional) `Vector{Float64}`: Length-`nlc` vector of accepted allele frequencies when `return_p = true`.

# Errors
- Throws an `ErrorException` if `maf <= 0` or `maf >= 0.5`.

# Examples
```julia
using FisherWright

gt = quickGT(100, 20)
size(gt) # (100, 20)

gt, freqs = quickGT(50, 10; maf = 0.05, return_p = true)
all(0.05 .< freqs .< 0.95) # true
```
"""
function quickGT(
    nlc::Int,
    nid::Int;
    maf = 0.1,
    qd = Beta(0.75, 0.75),
    rng = Random.default_rng(),
    return_p::Bool = false,
)
    (maf ≤ 0 || maf ≥ 0.5) && error("maf $maf not in (0, 0.5)")
    gt = Matrix{Int8}(undef, nlc, nid)
    freqs = return_p ? Vector{Float64}(undef, nlc) : nothing
    @inbounds for i = 1:nlc
        p = rand(rng, qd)
        while p ≤ maf || p ≥ 1 - maf
            p = rand(rng, qd)
        end
        freqs !== nothing && (freqs[i] = p)
        # Sample nid genotype counts for this locus
        row = rand(rng, Binomial(2, p), nid)
        @inbounds for j = 1:nid
            gt[i, j] = Int8(row[j])
        end
    end
    return return_p ? (gt, freqs) : gt
end

"""
    quickHap(
        nlc::Int,
        nid::Int;
        maf = 0.2,
        qd = Beta(0.75, 0.75),
        rng = Random.default_rng(),
        return_p::Bool = false,
    ) -> Matrix{Int8} | Tuple{Matrix{Int8}, Vector{Float64}}

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
- `Matrix{Int8}`: Haplotype matrix of shape `(nlc, 2 * nid)` with values in `{0, 1}`.
- (Optional) `Vector{Float64}`: Length-`nlc` vector of accepted allele frequencies when `return_p = true`.

# Errors
- Throws an `ErrorException` if `maf <= 0` or `maf >= 0.5`.

# Examples
```julia
using FisherWright

hp = quickHap(100, 20)
size(hp) # (100, 40)
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
    (maf ≤ 0 || maf ≥ 0.5) && error("maf $maf not in (0, 0.5)")
    nhp = 2 * nid
    hp = Matrix{Int8}(undef, nlc, nhp)
    freqs = return_p ? Vector{Float64}(undef, nlc) : nothing
    @inbounds for i = 1:nlc
        p = rand(rng, qd)
        while p ≤ maf || p ≥ 1 - maf
            p = rand(rng, qd)
        end
        freqs !== nothing && (freqs[i] = p)
        row = rand(rng, Binomial(1, p), nhp)
        @inbounds for j = 1:nhp
            hp[i, j] = Int8(row[j])
        end
    end
    return return_p ? (hp, freqs) : hp
end
