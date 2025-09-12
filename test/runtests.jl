using Test
using FisherWright

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
