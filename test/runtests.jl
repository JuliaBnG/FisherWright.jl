using Test
using Random
using FisherWright: fisher_wright, muts2bitarray, merge_sorted, recombine, cobp, FisherWrightResult, to_haplotype, _fixed_mutations, _drop_mutations, _fixation_step, _extract_fixed
using BnGStructs: Haplotype
using Distributions: Poisson

@testset "FisherWright basic" begin
    ne = 100
    nt = 200
    chr = [100_000_000, 100_000_000, 100_000_000]
    mr = 1.0
    mts, cbp = fisher_wright(ne, nt, chr, mr)
    @test length(cbp) == length(chr)
    @test length(mts) == 2ne
    # Basic invariant: haplotype mutation positions are sorted
    xy, lmp = muts2bitarray(mts, cbp)
    @test size(xy, 1) == lmp.pos |> length
    @test size(xy, 2) == 2ne
    @test all(0 .< lmp.frq .< 1)
end

@testset "Structured result boundary" begin
    ne = 12
    nt = 6
    chr = [10_000, 10_000]
    mr = 0.5
    result = fisher_wright(ne, nt, chr, mr; result=true, fixation_interval=1)

    @test result isa FisherWrightResult
    @test length(result.active_haplotypes) == 2ne
    @test result.chromosome_ends == UInt32.(cumsum(chr))

    manual = FisherWrightResult(
        [UInt32[1, 3], UInt32[2], UInt32[], UInt32[1, 2, 3]],
        UInt32[3],
        UInt32[5],
    )
    hap, loci = to_haplotype(manual)
    @test hap isa Haplotype
    @test hap.nlc == length(loci.pos)
    @test hap.nhp == 4
    @test hap.gt[1:hap.nlc, 1:hap.nhp] isa BitMatrix

    hap_fixed, loci_fixed = to_haplotype(manual; include_fixed=true)
    @test hap_fixed.nlc == 4
    @test length(loci_fixed.pos) == 4
    @test any(loci_fixed.pos .== UInt32(5))
end

@testset "Fixed mutation extraction" begin
    source = FisherWrightResult(
        [UInt32[1, 2, 5], UInt32[2, 5], UInt32[2, 5, 7], UInt32[2, 5]],
        UInt32[10],
        UInt32[11],
    )

    @test _fixed_mutations(source.active_haplotypes) == UInt32[2, 5]
    @test _drop_mutations(source.active_haplotypes, UInt32[2, 5]) ==
          [UInt32[1], UInt32[], UInt32[7], UInt32[]]

    cleaned, substitutions = _fixation_step(source.active_haplotypes, source.substitutions)
    @test cleaned == [UInt32[1], UInt32[], UInt32[7], UInt32[]]
    @test substitutions == UInt32[2, 5, 11]

    extracted = _extract_fixed(source)
    @test extracted.substitutions == UInt32[2, 5, 11]
    @test extracted.active_haplotypes == [UInt32[1], UInt32[], UInt32[7], UInt32[]]
    hap, loci = to_haplotype(extracted; include_fixed=true)
    @test hap.nlc == length(loci.pos)
    @test any(loci.pos .== UInt32(2))
    @test any(loci.pos .== UInt32(5))
end

@testset "Flat-vector invariants" begin
    @test merge_sorted(UInt32[1, 3, 5], UInt32[2, 3, 4]) ==
          UInt32[1, 2, 3, 4, 5]

    Random.seed!(1234)
    out = UInt32[]
    recombine(UInt32[1, 3, 5, 7], UInt32[2, 4, 6, 8], out, UInt32[4])
    @test issorted(out)
    @test length(out) == length(unique(out))

    parents = UInt32[1, 2, 3, 4]
    Random.seed!(1234)
    boundary = UInt32[]
    recombine(parents, parents, boundary, UInt32[3])
    @test boundary == parents

    Random.seed!(1234)
    cross = cobp(UInt32[3, 6], [Poisson(0.0), Poisson(0.0)])
    @test issorted(cross)
    @test all(1 .<= cross .<= UInt32(6))
end

@testset "Dense export boundary" begin
    muts = [
        UInt32[1, 3],
        UInt32[2],
        UInt32[],
        UInt32[1, 2, 3],
    ]
    cbp = UInt32[3]
    xy, lmp = muts2bitarray(muts, cbp)

    @test size(xy) == (3, 4)
    @test length(lmp.pos) == 3

    hp = Haplotype(xy)
    @test hp.nlc == 3
    @test hp.nhp == 4
    @test hp.gt[1:hp.nlc, 1:hp.nhp] == xy
end
