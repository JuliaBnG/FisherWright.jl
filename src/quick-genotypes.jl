"""
    quickGT(nlc::Int, nid::Int; maf=0.1, qd=Beta(0.75,0.75), rng=Random.default_rng(), return_p=false)

Simulate diploid genotypes (nlc loci × nid individuals) with allele frequencies
drawn from `qd` (default Beta(0.75,0.75)), rejecting values with p ≤ maf or
p ≥ 1 - maf. Returns Matrix{Int8} (values 0,1,2). If `return_p=true`, also
returns the Vector{Float64} of accepted allele frequencies.
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
    quickHap(nlc::Int, nid::Int; maf=0.2, qd=Beta(0.75,0.75), rng=Random.default_rng(), return_p=false)

Simulate haplotypes: returns Matrix{Int8} of size (nlc × 2nid) with 0/1 alleles.
Allele frequencies drawn as in `quickGT`. If `return_p=true`, also return
frequency vector.
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
