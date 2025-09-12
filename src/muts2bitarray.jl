"""
    function muts2bitarray(muts::Vector{Vector{UInt32}}, cbp::Vector{UInt32})
Convert mutations in haplotypes to a `BitArray` of haplotypes and a DataFrame of
linkage map which contains the chromosome number, position, reference and
alternative alleles, and frequency of each mutation. When `flip` is `true`, the
allele names at randomly selected loci are swapped. The loci that are fixed are
 removed.
"""
function muts2bitarray(muts::Vector{Vector{UInt32}}, cbp::Vector{UInt32}; flip = false)
    all_mts = UInt32[]
    for t in muts
        all_mts = merge_sorted(all_mts, t)
    end
    idx = Dict(mut => i for (i, mut) in enumerate(all_mts))
    nlc, nhp = length(all_mts), length(muts)
    xy = BitArray(undef, nlc, nhp)
    xy .= 0

    Threads.@threads for i ∈ 1:nhp
        ix = [idx[mut] for mut in muts[i]]
        view(xy, ix, i) .= 1
    end
    if flip
        Threads.@threads for i ∈ (1:nlc)[rand(Bool, nlc)]
            xy[i, :] = .!xy[i, :]
        end
    end

    chr = map(all_mts) do mut
        searchsortedfirst(cbp, mut) # find the chromosome number
    end
    ref = rand("ACGT", nlc)
    alt = rand.(setdiff.("ACGT", ref))
    frq = vec(sum(xy, dims = 2))
    oo = 0 .< frq .< nhp # remove fixed loci
    frq /= nhp
    lmp = DataFrame(chr = Int8.(chr), pos = all_mts, ref = ref, alt = alt, frq = frq)
    return xy[oo, :], lmp[oo, :]
end

