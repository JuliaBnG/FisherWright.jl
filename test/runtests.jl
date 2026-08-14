using Test
using Random
using FisherWright: fisher_wright, muts2bitarray, merge_sorted, merge_sorted!, recombine,
    cobp, cobp!, FisherWrightResult, to_haplotype, RecombinationMap,
    uniform_recombination_map, random_mate, random_mate!, _fixed_mutations,
    _drop_mutations, _drop_mutations!, _fixation_step, _fixation_step!, _extract_fixed,
    _intersect_sorted!, _remove_sorted!
using BnGStructs: Haplotype, Species, Cattle, GenericSpecies
using Distributions: Poisson

"""
Reference implementation of meiosis: emit the positions of alternating parents
between successive crossovers. Deliberately naive, used as the correctness
oracle for `recombine`.
"""
function ref_recombine(h₁, h₂, cross_overs, start_from_h₁::Bool)
    out = UInt32[]
    o = start_from_h₁
    prev = UInt32(0)
    for b in vcat(cross_overs, typemax(UInt32))
        for v in (o ? h₁ : h₂)
            (prev <= v < b) && push!(out, v)
        end
        prev = b
        o = !o
    end
    return out
end

# `recombine` picks its starting haplotype internally, so a correct result is
# one of the two possible phases.
function recombine_matches_reference(h₁, h₂, cross_overs, out)
    return out == ref_recombine(h₁, h₂, cross_overs, true) ||
           out == ref_recombine(h₁, h₂, cross_overs, false)
end

@testset "FisherWright basic" begin
    ne = 100
    nt = 200
    chr = [100_000_000, 100_000_000, 100_000_000]
    mr = 1.0
    mts, cbp = fisher_wright(ne, nt, chr, mr)
    @test length(cbp) == length(chr)
    @test length(mts) == 2ne
    # Basic invariant: haplotype mutation positions are sorted and unique
    @test all(issorted, mts)
    @test all(h -> length(h) == length(unique(h)), mts)
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
    @test all(issorted, result.active_haplotypes)
    @test issorted(result.substitutions)
    # A fixed position must not remain in the active haplotypes.
    @test all(h -> isempty(intersect(h, result.substitutions)), result.active_haplotypes)

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
    # the non-mutating form must leave its input alone
    @test source.active_haplotypes[1] == UInt32[1, 2, 5]

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

    # in-place forms agree with the copying forms
    haps = [UInt32[1, 2, 5], UInt32[2, 5], UInt32[2, 5, 7], UInt32[2, 5]]
    updated = _fixation_step!(haps, UInt32[11])
    @test haps == [UInt32[1], UInt32[], UInt32[7], UInt32[]]
    @test updated == UInt32[2, 5, 11]

    haps2 = [UInt32[1, 2, 3], UInt32[2, 3, 4]]
    @test _drop_mutations!(haps2, UInt32[2]) === haps2
    @test haps2 == [UInt32[1, 3], UInt32[3, 4]]
    # nothing fixed: haplotypes and substitutions pass through untouched
    @test _fixation_step!([UInt32[1], UInt32[2]], UInt32[7]) == UInt32[7]
end

@testset "Sorted set operations" begin
    dest = UInt32[]
    @test _intersect_sorted!(dest, UInt32[1, 3, 5, 7], UInt32[3, 4, 5, 9]) == UInt32[3, 5]
    @test _intersect_sorted!(dest, UInt32[], UInt32[1, 2]) == UInt32[]
    @test _intersect_sorted!(dest, UInt32[1, 2], UInt32[]) == UInt32[]
    @test _intersect_sorted!(dest, UInt32[2, 4], UInt32[1, 3, 5]) == UInt32[]
    # exercise the binary-search branch (short candidate list, long haplotype)
    long = UInt32.(1:1000)
    @test _intersect_sorted!(dest, UInt32[7, 500, 1001], long) == UInt32[7, 500]

    @test _remove_sorted!(UInt32[1, 2, 3, 4], UInt32[2, 4]) == UInt32[1, 3]
    @test _remove_sorted!(UInt32[1, 2, 3], UInt32[]) == UInt32[1, 2, 3]
    @test _remove_sorted!(UInt32[], UInt32[1]) == UInt32[]
    @test _remove_sorted!(UInt32[1, 2], UInt32[1, 2]) == UInt32[]
    @test _remove_sorted!(UInt32[5, 6], UInt32[1, 2]) == UInt32[5, 6]

    # random cross-check against Base's hash-based set operations
    Random.seed!(99)
    for _ = 1:500
        a = sort(unique(rand(UInt32(1):UInt32(60), rand(0:25))))
        b = sort(unique(rand(UInt32(1):UInt32(60), rand(0:25))))
        @test _intersect_sorted!(dest, a, b) == sort(intersect(a, b))
        @test _remove_sorted!(copy(a), b) == sort(setdiff(a, b))
    end
end

@testset "merge_sorted" begin
    @test merge_sorted(UInt32[1, 3, 5], UInt32[2, 3, 4]) == UInt32[1, 2, 3, 4, 5]
    @test merge_sorted(UInt32[], UInt32[2, 2, 3]) == UInt32[2, 3]
    @test merge_sorted(UInt32[1, 1], UInt32[]) == UInt32[1]

    # in-place form matches, and reuses its buffer
    dest = UInt32[]
    @test merge_sorted!(dest, UInt32[1, 3, 5], UInt32[2, 3, 4]) == UInt32[1, 2, 3, 4, 5]
    @test merge_sorted!(dest, UInt32[9], UInt32[8]) == UInt32[8, 9]
    @test dest == UInt32[8, 9]
    v = UInt32[1, 2]
    @test_throws ArgumentError merge_sorted!(v, v, UInt32[3])

    Random.seed!(7)
    for _ = 1:500
        a = sort(unique(rand(UInt32(1):UInt32(40), rand(0:20))))
        b = sort(unique(rand(UInt32(1):UInt32(40), rand(0:20))))
        @test merge_sorted(a, b) == sort(union(a, b))
    end
end

@testset "recombine correctness" begin
    # Regression: the last position of each parent used to be reachable only
    # through the trailing append, which lost it or emitted it out of order.
    h₁, h₂, co = UInt32[5, 14, 18, 26], UInt32[10], UInt32[25]
    out = UInt32[]
    for _ = 1:50
        recombine(h₁, h₂, out, co)
        @test recombine_matches_reference(h₁, h₂, co, out)
        @test issorted(out)
    end

    # Randomised property test over both phases.
    Random.seed!(20240)
    out = UInt32[]
    for _ = 1:20_000
        h₁ = sort(unique(rand(UInt32(1):UInt32(50), rand(0:6))))
        h₂ = sort(unique(rand(UInt32(1):UInt32(50), rand(0:6))))
        co = sort(unique(rand(UInt32(1):UInt32(50), rand(0:3))))
        recombine(h₁, h₂, out, co)
        @test recombine_matches_reference(h₁, h₂, co, out)
        @test issorted(out)
        @test length(out) == length(unique(out))
        @test all(p -> p in h₁ || p in h₂, out)
    end

    # Both phases must actually occur, otherwise the property test above is
    # only ever checking one of them.
    Random.seed!(5)
    parents₁, parents₂ = UInt32[1, 2, 3], UInt32[4, 5, 6]
    seen = Set{Vector{UInt32}}()
    for _ = 1:100
        recombine(parents₁, parents₂, out, UInt32[])
        push!(seen, copy(out))
    end
    @test seen == Set([parents₁, parents₂])

    # Degenerate inputs.
    @test recombine(UInt32[], UInt32[], out, UInt32[3]) == UInt32[]
    @test recombine(UInt32[1, 2, 3], UInt32[1, 2, 3], out, UInt32[2]) == UInt32[1, 2, 3]

    # An explicit rng makes the phase reproducible.
    a = recombine(UInt32[1, 3], UInt32[2, 4], UInt32[], UInt32[2]; rng=MersenneTwister(1))
    b = recombine(UInt32[1, 3], UInt32[2, 4], UInt32[], UInt32[2]; rng=MersenneTwister(1))
    @test a == b
end

@testset "Flat-vector invariants" begin
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

@testset "Recombination map" begin
    map = uniform_recombination_map([3, 3]; M=1e8)
    @test map isa RecombinationMap
    @test map.cbp == UInt32[3, 6]
    @test map.interval_ends == [UInt32[3], UInt32[6]]
    @test map.interval_rates[1][1] == 3 / 1e8
    @test length(map.samplers) == 2
    @test_throws ArgumentError uniform_recombination_map([3, 3]; M=0)

    Random.seed!(1234)
    cross_map = cobp(map)
    @test issorted(cross_map)
    @test all(1 .<= cross_map .<= UInt32(6))

    custom = RecombinationMap(
        UInt32[6],
        [UInt32[3, 6]],
        [Float64[0.0, 0.0]],
    )
    Random.seed!(1234)
    cross_custom = cobp(custom)
    @test issorted(cross_custom)
    @test all(1 .<= cross_custom .<= UInt32(6))

    # cobp! reuses its buffer and agrees with cobp
    buf = UInt32[]
    @test cobp!(buf, map; rng=MersenneTwister(3)) == cobp(map; rng=MersenneTwister(3))
    @test cobp!(buf, map; rng=MersenneTwister(4)) == cobp(map; rng=MersenneTwister(4))
    @test buf === cobp!(buf, map)

    # Chromosome ends segregate independently: each is emitted about half the
    # time, and crossovers stay sorted with a dense map.
    dense = uniform_recombination_map([1_000_000, 1_000_000]; M=1e6)
    rng = MersenneTwister(11)
    boundary_hits = 0
    for _ = 1:2000
        cobp!(buf, dense; rng=rng)
        @test issorted(buf)
        UInt32(1_000_000) in buf && (boundary_hits += 1)
    end
    @test 800 < boundary_hits < 1200
end

@testset "random_mate" begin
    Random.seed!(3)
    pm, sex = random_mate(20, 30)
    @test size(pm) == (30, 2)
    @test length(sex) == 20
    @test all(i -> sex[i], pm[:, 1])
    @test all(i -> !sex[i], pm[:, 2])

    # in-place form fills the supplied buffers
    pm2 = Matrix{Int}(undef, 30, 2)
    sex2 = Vector{Bool}(undef, 20)
    out_pm, out_sex = random_mate!(pm2, sex2, Int[], Int[]; rng=MersenneTwister(2))
    @test out_pm === pm2
    @test out_sex === sex2
    @test all(i -> sex2[i], pm2[:, 1])
    @test all(i -> !sex2[i], pm2[:, 2])
end

@testset "Rate parameters are independent" begin
    # Sized so the carried-mutation counts concentrate (~20 000 copies): the
    # ratios below are then stable to a few percent across seeds and thread
    # counts, which a smaller run is far too noisy to give.
    chr = [100_000_000, 100_000_000]
    Random.seed!(21)
    a, _ = fisher_wright(100, 50, chr, 1.0; M=1e8)
    Random.seed!(21)
    b, _ = fisher_wright(100, 50, chr, 1.0; M=1e6)
    Random.seed!(21)
    c, _ = fisher_wright(100, 50, chr, 1.0; mut_base=1e7)
    na, nb, nc = sum(length, a), sum(length, b), sum(length, c)

    # M drives recombination only: a hundredfold denser crossover map must
    # leave the mutation rate alone.
    @test 0.75 < na / nb < 1.3
    # mut_base drives mutation only: a tenfold smaller basis is tenfold the
    # rate.
    @test 7.5 < nc / na < 13.5

    @test_throws ArgumentError fisher_wright(30, 5, chr, 1.0; mut_base=0.0)
end

@testset "Simulation invariants at scale" begin
    Random.seed!(31)
    res = fisher_wright(60, 80, [5_000_000, 5_000_000, 5_000_000], 4.0;
        result=true, fixation_interval=5)
    @test all(issorted, res.active_haplotypes)
    @test all(h -> length(h) == length(unique(h)), res.active_haplotypes)
    @test issorted(res.substitutions)
    @test length(res.substitutions) == length(unique(res.substitutions))
    @test all(h -> isempty(intersect(h, res.substitutions)), res.active_haplotypes)
    @test all(h -> all(p -> 1 <= p <= 15_000_000, h), res.active_haplotypes)
end

@testset "Silent by default" begin
    chr = [100_000, 100_000]

    function captured_stdout(f)
        pipe = Pipe()
        Base.link_pipe!(pipe; reader_supports_async=true, writer_supports_async=true)
        reader = @async read(pipe, String)
        try
            redirect_stdout(f, pipe)
        finally
            close(pipe.in)
        end
        return fetch(reader)
    end

    @test isempty(captured_stdout(() -> fisher_wright(10, 120, chr, 1.0)))
    @test occursin("Generation", captured_stdout(
        () -> fisher_wright(10, 120, chr, 1.0; verbose=true)))
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

@testset "Species-based simulation" begin
    # Test using a small generic species
    sp = GenericSpecies("Mini", Int32(20), UInt32[100_000, 100_000], UInt32(50_000_000))
    res = fisher_wright(sp, 30, 1.0; result=true)
    @test res isa FisherWrightResult
    @test length(res.active_haplotypes) == 40
    @test res.chromosome_ends == UInt32[100_000, 200_000]

    # Test default tuple return and default mr=1.0
    mts, cbp = fisher_wright(sp, 10)
    @test length(mts) == 40
    @test cbp == UInt32[100_000, 200_000]

    # Test Cattle instance
    cattle = Cattle(15)
    res_cattle = fisher_wright(cattle, 10, 1.0; result=true)
    @test res_cattle isa FisherWrightResult
    @test length(res_cattle.active_haplotypes) == 30
    @test length(res_cattle.chromosome_ends) == 29

    # Recombination map from Species
    rmap = uniform_recombination_map(cattle)
    @test rmap isa RecombinationMap
    @test length(rmap.cbp) == 29
end
