"""
    quickGT(nlc::Int, nid::Int; maf = 0.1, qd = Beta(0.75, 0.75))
A quick way to simulate SNP genotypes of `nlc` loci and `nid` individuals.
Allele frequencies are sampled from a default `Beta(.75, .75)`, conditioned on
`maf`. This distribution is U-shaped, which is like real situations. The loci,
however, are independently sampled, so there is no linkage disequilibrium among
loci.

Please refer [`fisher_wright`](@ref) for genotypes of Fisher-Wright
population.
"""
function quickGT(nlc::Int, nid::Int; maf = 0.1, qd = Beta(0.75, 0.75))
    (maf ≤ 0 || maf ≥ 0.5) && error("maf $maf not in (0, 0.5)")
    gt = zeros(Int8, nlc, nid)
    for iic = 1:nlc
        p = 0
        while p <= maf || p >= 1 - maf
            p = rand(qd)
        end
        rand!(Binomial(2, p), view(gt, iic, :))
    end
    gt
end

"""
    function quickHap(nlc::Int, nid::Int; maf = .2, bp = .75)
A quick way to simulate SNP genotype of `nlc` loci, and `2nid` haplotypes.
Allele Frequencies are sampled from a default distribution of `Beta(.75, .75)`,
conditioned on `maf`. This distribution is U-shaped, which is like real
situations. The loci, however, are independently sampled, so there is no linkage
disequilibrium among loci.

Please refer [`fisher_wright`](@ref) for genotypes of Fisher-Wright
population.
"""
function quickHap(nlc, nid; maf = 0.2, qd = Beta(0.75, 0.75))
    (maf ≤ 0 || maf ≥ 0.5) && error("maf $maf not in (0, 0.5)")
    nhp = 2nid
    hp = zeros(Int8, nlc, nhp)
    for iic = 1:nlc
        p = 0
        while p <= maf || p >= 1 - maf
            p = rand(qd)
        end
        rand!(Binomial(1, p), view(hp, iic, :))
    end
    hp
end
