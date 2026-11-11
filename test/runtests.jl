using Test
using Random
using FisherWright: fisher_wright, muts2bitarray, extract_chip_bitarray, merge_sorted, merge_sorted!, recombine, _recombine!, _insert_sorted!, _isearchfirst,
    cobp, cobp!, FisherWrightResult, to_haplotype, RecombinationMap,
    uniform_recombination_map, random_mate, random_mate!, _fixed_sex_mate!, _nmale, _fixed_mutations,
    _drop_mutations, _drop_mutations!, _fixation_step, _fixation_step!, _extract_fixed,
    _intersect_sorted!, _remove_sorted!, _with_fixed, _simulate!, _initial_state,
    SimState, ConstantSize, SizeTrajectory, RandomMating, Neutral
using BnGStructs
using Distributions: Poisson
using Statistics: mean, var, cor

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

@testset "Zero heap allocations in meiosis inner loop operations" begin
    chr = [10_000_000, 10_000_000]
    rmap = uniform_recombination_map(chr)
    cbuf = zeros(UInt32, 100)
    empty!(cbuf)
    cobp!(cbuf, rmap) # warm up
    @test @allocated(cobp!(cbuf, rmap)) == 0

    h1 = sort(rand(UInt32(1):UInt32(20_000_000), 500))
    h2 = sort(rand(UInt32(1):UInt32(20_000_000), 500))
    ho = zeros(UInt32, 2000)
    empty!(ho)
    recombine(h1, h2, ho, cbuf) # warm up
    @test @allocated(recombine(h1, h2, ho, cbuf)) == 0
    _recombine!(h1, h2, ho, cbuf, true) # warm up
    @test @allocated(_recombine!(h1, h2, ho, cbuf, true)) == 0

    muts = [UInt32(15), UInt32(30)]
    buf = zeros(UInt32, 2000)
    empty!(buf)
    merge_sorted!(buf, ho, muts) # warm up
    @test @allocated(merge_sorted!(buf, ho, muts)) == 0
    _insert_sorted!(buf, ho, muts) # warm up
    @test @allocated(_insert_sorted!(buf, ho, muts)) == 0
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
        UInt32[5],
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
        UInt32[11],
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

@testset "_isearchfirst matches searchsortedfirst" begin
    ref(v, x, lo, hi) = searchsortedfirst(v, x, lo, hi, Base.Order.Forward)
    @test _isearchfirst(UInt32[], UInt32(5), 1, 0) == 1
    v = UInt32[10, 20, 30]
    @test _isearchfirst(v, UInt32(5), 1, 3) == 1
    @test _isearchfirst(v, UInt32(30), 1, 3) == 3
    @test _isearchfirst(v, UInt32(31), 1, 3) == 4
    @test _isearchfirst(v, UInt32(25), 3, 2) == 3
    Random.seed!(5)
    G = UInt32(1_000_000_000)
    for _ = 1:300
        # Uniform, clustered (most positions in a narrow band) and with repeats.
        L = rand(0:5000)
        v = sort!(rand(UInt32(1):G, L))
        w = sort!(vcat(rand(UInt32(1):UInt32(1000), L), rand(UInt32(1):G, 5)))
        u = sort!(rand(UInt32(1):UInt32(50), L))
        for h in (v, w, u), _ = 1:50
            n = length(h)
            lo = rand(1:n+1)
            hi = rand(lo-1:n)
            x = isempty(h) || rand(Bool) ? rand(UInt32(0):G) : h[rand(1:n)]
            @test _isearchfirst(h, x, lo, hi) == ref(h, x, lo, hi)
        end
    end
    # Measured inside a function: at testset scope, Julia 1.10 counts the
    # boxing of the returned Int against the call.
    alloc(v, x) = @allocated _isearchfirst(v, x, 1, length(v))
    h = sort!(rand(UInt32(1):G, 10_000))
    alloc(h, UInt32(500_000_000)) # warm up
    @test alloc(h, UInt32(500_000_000)) == 0
end

@testset "_insert_sorted! matches merge_sorted!" begin
    ins! = _insert_sorted!
    @test ins!(UInt32[], UInt32[], UInt32[]) == UInt32[]
    @test ins!(UInt32[], UInt32[], UInt32[2, 2, 5]) == UInt32[2, 5]
    @test ins!(UInt32[], UInt32[3, 7], UInt32[]) == UInt32[3, 7]
    @test ins!(UInt32[], UInt32[3, 7], UInt32[1, 3, 3, 7, 9, 9]) == UInt32[1, 3, 7, 9]
    v = UInt32[1, 2]
    @test_throws ArgumentError ins!(v, v, UInt32[3])

    Random.seed!(11)
    dest, ref = UInt32[], UInt32[]
    for _ = 1:2000
        v = sort!(unique!(rand(UInt32(1):UInt32(60), rand(0:40))))
        # Short, sorted, with repeats and elements of v.
        ins = sort!(vcat(rand(UInt32(1):UInt32(60), rand(0:6)), isempty(v) ? UInt32[] : rand(v, rand(0:2))))
        @test ins!(dest, v, ins) == merge_sorted!(ref, v, ins)
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

    # The kernel with an explicit strand reproduces the reference phase exactly,
    # and `recombine` is the kernel after exactly one `rand(rng, Bool)` draw,
    # so the simulation's RNG draw order is unchanged by the split.
    Random.seed!(20241)
    for _ = 1:5_000
        h₁ = sort(unique(rand(UInt32(1):UInt32(50), rand(0:6))))
        h₂ = sort(unique(rand(UInt32(1):UInt32(50), rand(0:6))))
        co = sort(unique(rand(UInt32(1):UInt32(50), rand(0:3))))
        for o in (true, false)
            @test _recombine!(h₁, h₂, out, co, o) == ref_recombine(h₁, h₂, co, o)
        end
        seed = rand(UInt64)
        r₁, r₂ = Xoshiro(seed), Xoshiro(seed)
        @test recombine(h₁, h₂, UInt32[], co; rng=r₁) == _recombine!(h₁, h₂, out, co, rand(r₂, Bool))
        @test rand(r₁, UInt64) == rand(r₂, UInt64)
    end
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
        UInt32(1_000_001) in buf && (boundary_hits += 1)
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

@testset "Mutation switch" begin
    chr = [100_000, 200_000]
    args = (12, 40, chr, 5.0)
    kw = (M=5e4, mut_base=1e5)

    # Without mutation a monomorphic start stays monomorphic, in both return
    # modes and through the Species method.
    haps, cbp = fisher_wright(args...; kw..., mutate=false, rng=Xoshiro(1))
    @test length(haps) == 24 && all(isempty, haps)
    @test cbp == UInt32.(cumsum(chr))
    res = fisher_wright(args...; kw..., mutate=false, result=true, rng=Xoshiro(1))
    @test all(isempty, res.active_haplotypes) && isempty(res.substitutions)
    sp = GenericSpecies("NoMut", Int32(8), UInt32.(chr), UInt32(50_000))
    @test all(isempty, fisher_wright(sp, 10; mutate=false, rng=Xoshiro(2))[1])

    # mutate = true is the default and the switch does not alter it.
    @test fisher_wright(args...; kw..., rng=Xoshiro(3))[1] ==
          fisher_wright(args...; kw..., mutate=true, rng=Xoshiro(3))[1]
    @test !all(isempty, fisher_wright(args...; kw..., rng=Xoshiro(3))[1])

    # mr is still validated when it is not used.
    @test_throws ArgumentError fisher_wright(12, 4, chr, 0.0; mutate=false)
    @test_throws ArgumentError fisher_wright(12, 4, chr, 25.0; mutate=false)
end

@testset "Starting population (init)" begin
    chr = [100_000, 200_000]
    ne, tbp = 12, UInt32(300_000)
    kw = (M=5e4, mut_base=1e5)
    burn = fisher_wright(ne, 300, chr, 5.0; kw..., result=true, rng=Xoshiro(1))
    @test !isempty(burn.substitutions) && !all(isempty, burn.active_haplotypes)
    full = _with_fixed(burn)

    # The same continuation whether the state is given as a result or as plain
    # haplotypes, in both return modes and for any fixation interval.
    a = fisher_wright(ne, 40, chr, 5.0; kw..., init=burn, result=true, rng=Xoshiro(2))
    b, cbp = fisher_wright(ne, 40, chr, 5.0; kw..., init=burn, rng=Xoshiro(2))
    c, _ = fisher_wright(ne, 40, chr, 5.0; kw..., init=full, rng=Xoshiro(2))
    d = fisher_wright(ne, 40, chr, 5.0; kw..., init=full, result=true, fixation_interval=7, rng=Xoshiro(2))
    @test cbp == burn.chromosome_ends
    @test _with_fixed(a) == b == c == _with_fixed(d)
    @test a.substitutions == d.substitutions
    @test issubset(burn.substitutions, a.substitutions)
    @test length(a.substitutions) > length(burn.substitutions)
    @test all(h -> isempty(intersect(h, a.substitutions)), a.active_haplotypes)

    # init is neither modified nor aliased by the run.
    snap = deepcopy(burn)
    r = fisher_wright(ne, 3, chr, 5.0; kw..., init=burn, result=true, rng=Xoshiro(3))
    @test burn.active_haplotypes == snap.active_haplotypes
    @test burn.substitutions == snap.substitutions
    @test r.substitutions !== burn.substitutions
    @test !any(x -> any(y -> x === y, burn.active_haplotypes), r.active_haplotypes)
    snap = deepcopy(full)
    h, _ = fisher_wright(ne, 3, chr, 5.0; kw..., init=full, rng=Xoshiro(3))
    @test full == snap
    @test !any(x -> any(y -> x === y, full), h)

    # Drift and recombination only: no new site appears, every site carried
    # by all haplotypes stays, and the substitutions are those of init plus
    # sites lost from segregation by fixation.
    nomut = fisher_wright(ne, 60, chr, 5.0; kw..., init=burn, mutate=false, result=true, rng=Xoshiro(4))
    seg₀ = reduce(union, burn.active_haplotypes)
    @test issubset(reduce(union, nomut.active_haplotypes), seg₀)
    @test issubset(burn.substitutions, nomut.substitutions)
    @test issubset(setdiff(nomut.substitutions, burn.substitutions), seg₀)
    @test length(reduce(union, nomut.active_haplotypes)) < length(seg₀)

    # Without crossovers (one chromosome, 1e-7 Morgan) each offspring
    # haplotype is a copy of one parental haplotype.
    start = [sort!(unique!(rand(Xoshiro(100 + i), UInt32(1):UInt32(1000), 20))) for i = 1:2ne]
    kids, _ = fisher_wright(ne, 1, [1000], 1.0; M=1e10, init=start, mutate=false, rng=Xoshiro(5))
    @test all(in(start), kids)
    # Individual i is haplotypes 2i-1 and 2i: its two gametes come from a
    # sire and a dam, never from the same parent.
    owner(k) = (findfirst(==(k), start) + 1) ÷ 2
    @test allunique(start)
    @test all(i -> owner(kids[2i-1]) != owner(kids[2i]), 1:ne)

    # A site carried by every starting haplotype is a substitution.
    fixedstart = [UInt32[7, 50 + i] for i = 1:2ne]
    res = fisher_wright(ne, 1, [1000], 1.0; init=fixedstart, mutate=false, result=true, rng=Xoshiro(6))
    @test res.substitutions == UInt32[7]
    @test all(h -> !(UInt32(7) in h), res.active_haplotypes)

    # Species method.
    sp = GenericSpecies("Init", Int32(ne), UInt32.(chr), UInt32(50_000))
    @test fisher_wright(sp, 40, 5.0; mut_base=1e5, init=burn, rng=Xoshiro(2))[1] == b

    # Invalid starting populations.
    empty_start() = [UInt32[] for _ = 1:2ne]
    bad(f) = (v = empty_start(); f(v); v)
    run(init; n=ne, c=chr) = fisher_wright(n, 2, c, 5.0; init=init)
    @test_throws ArgumentError run(empty_start()[1:end-1])
    @test_throws ArgumentError run(burn; n=ne + 1)
    @test_throws ArgumentError run(bad(v -> v[3] = UInt32[5, 4]))
    @test_throws ArgumentError run(bad(v -> v[3] = UInt32[5, 5]))
    @test_throws ArgumentError run(bad(v -> v[3] = UInt32[0, 5]))
    @test_throws ArgumentError run(bad(v -> v[end] = UInt32[5, tbp + 1]))
    @test run(bad(v -> v[end] = UInt32[1, tbp]))[1] isa Vector{Vector{UInt32}}
    @test_throws ArgumentError run(burn; c=[200_000, 100_000])
    @test_throws ArgumentError run(FisherWrightResult(empty_start(), UInt32[100_000, 300_000], UInt32[9, 3]))
end

@testset "Simulation engine" begin
    chr = [1_000_000, 1_000_000]
    ne, cbp = 50, UInt32[1_000_000, 2_000_000]
    kw = (M=1e6, mut_base=1e6)
    rmap = uniform_recombination_map(chr; M=1e6)
    span = UInt32(1):UInt32(2_000_000)
    engine!(st, rng, nt; mutate=true, track=false) = _simulate!(
        st, ConstantSize(ne), RandomMating(), Neutral(), nt, rmap, Poisson(10.0), span;
        mutate=mutate, track=track, fixation_interval=1, verbose=false,
        progress_interval=100, rng=rng)

    # The state takes exactly one seed per haplotype slot from rng, in order.
    rng, ref = Xoshiro(5), Xoshiro(5)
    st = SimState(_initial_state(nothing, 2ne, cbp, false)..., ne, rng)
    seeds = [rand(ref, UInt64) for _ = 1:2ne]
    @test st.srng == Xoshiro.(seeds)
    @test rand(rng, UInt64) == rand(ref, UInt64)
    @test length(st.prt) == length(st.off) == length(st.mbuf) == length(st.nbuf) == 2ne
    @test length(st.cbuf) == length(st.sex) == size(st.pm, 1) == ne

    # fisher_wright is the engine behind validation: same seed, same output.
    rng = Xoshiro(9)
    st = SimState(_initial_state(nothing, 2ne, cbp, true)..., ne, rng)
    @test engine!(st, rng, 120; track=true) === st
    res = fisher_wright(ne, 120, chr, 5.0; kw..., result=true, rng=Xoshiro(9))
    @test st.prt == res.active_haplotypes && st.substitutions == res.substitutions

    # Inference: the engine returns its state, and fisher_wright one of its
    # two documented return types, with no abstract element type.
    @test (@inferred engine!(st, rng, 1)) === st
    f(r) = fisher_wright(ne, 2, chr, 5.0; result=r)
    @test only(Base.return_types(f, (Bool,))) ==
          Union{FisherWrightResult,Tuple{Vector{Vector{UInt32}},Vector{UInt32}}}

    # Steady state: a generation allocates task overhead only, far less than
    # one copy of the haplotypes (over 1 MB here), with or without mutation.
    rng = Xoshiro(11)
    st = SimState(_initial_state(nothing, 2ne, cbp, false)..., ne, rng)
    engine!(st, rng, 400)
    @test sum(length, st.prt) > 250_000
    for mutate in (true, false)
        engine!(st, rng, 50; mutate=mutate)
        @test @allocated(engine!(st, rng, 50; mutate=mutate)) / 50 < 2^16
    end

    # The same holds while the census size changes: buffers are sized for the
    # largest cohort, and those of a smaller one are kept.
    demog = SizeTrajectory(ne, repeat([20, 35, 50, 35], 25))
    cycle!(st, rng) = _simulate!(
        st, demog, RandomMating(), Neutral(), 100, rmap, Poisson(10.0), span;
        mutate=false, track=false, fixation_interval=1, verbose=false,
        progress_interval=100, rng=rng)
    cycle!(st, rng)
    @test st.n == 35 && length(st.prt) == 2ne
    @test @allocated(cycle!(st, rng)) / 100 < 2^16

    # ... and with a fixed number of males.
    fixed!(st, rng) = _simulate!(
        st, demog, RandomMating(nmale=7), Neutral(), 100, rmap, Poisson(10.0), span;
        mutate=false, track=false, fixation_interval=1, verbose=false,
        progress_interval=100, rng=rng)
    fixed!(st, rng)
    @test sum(st.sex) == 7 && length(st.sex) == 50   # parents of the last generation
    @test @allocated(fixed!(st, rng)) / 100 < 2^16
end

@testset "Changing population size" begin
    chr = [100_000, 200_000]
    args = (20, 60, chr, 5.0)
    kw = (M=5e4, mut_base=1e5)
    valid(h, G) = issorted(h) && allunique(h) && all(p -> 1 <= p <= G, h)

    # A trajectory that stays at ne is the run without popsize, as a vector
    # and as a function.
    a = fisher_wright(args...; kw..., rng=Xoshiro(3))
    @test fisher_wright(args...; kw..., popsize=fill(20, 60), rng=Xoshiro(3)) == a
    @test fisher_wright(args...; kw..., popsize=g -> 20, rng=Xoshiro(3)) == a
    @test fisher_wright(args...; kw..., popsize=fill(Int32(20), 60), rng=Xoshiro(3)) == a

    # Bottleneck, expansion and an odd final size: 2N_nt haplotypes come
    # back, both return modes agree, and a vector equals the same function.
    sizes = [fill(5, 20); fill(45, 39); 13]
    res = fisher_wright(args...; kw..., popsize=sizes, result=true, rng=Xoshiro(4))
    haps, cbp = fisher_wright(args...; kw..., popsize=sizes, rng=Xoshiro(4))
    @test length(haps) == length(res.active_haplotypes) == 26
    @test cbp == res.chromosome_ends == UInt32[100_000, 300_000]
    @test _with_fixed(res) == haps
    @test fisher_wright(args...; kw..., popsize=g -> sizes[g], rng=Xoshiro(4))[1] == haps
    @test fisher_wright(args...; kw..., popsize=sizes, result=true, fixation_interval=7,
        rng=Xoshiro(4)).substitutions == res.substitutions
    @test !isempty(res.substitutions) && !all(isempty, res.active_haplotypes)
    @test all(h -> valid(h, 300_000), res.active_haplotypes)
    @test isempty(_fixed_mutations(res.active_haplotypes))
    @test all(h -> isempty(intersect(h, res.substitutions)), res.active_haplotypes)

    # The smallest population, and growth from it.
    tiny = fisher_wright(2, 80, chr, 5.0; kw..., popsize=fill(2, 80), result=true, rng=Xoshiro(5))
    @test length(tiny.active_haplotypes) == 4
    grown, _ = fisher_wright(2, 6, chr, 5.0; kw..., popsize=g -> 2^(g + 1), rng=Xoshiro(5))
    @test length(grown) == 2 * 2^7 && all(h -> valid(h, 300_000), grown)

    # With init, ne is the size of the starting population; without mutation
    # an expansion only redistributes its variants.
    start = fisher_wright(8, 100, chr, 5.0; kw..., result=true, rng=Xoshiro(6))
    big = fisher_wright(8, 5, chr, 5.0; kw..., init=start, popsize=fill(64, 5),
        mutate=false, result=true, rng=Xoshiro(7))
    @test length(big.active_haplotypes) == 128
    @test issubset(reduce(union, big.active_haplotypes), reduce(union, start.active_haplotypes))
    @test issubset(start.substitutions, big.substitutions)
    @test start.active_haplotypes == fisher_wright(8, 100, chr, 5.0; kw..., result=true,
        rng=Xoshiro(6)).active_haplotypes   # init left untouched

    # Species method.
    sp = GenericSpecies("Size", Int32(20), UInt32.(chr), UInt32(50_000))
    @test fisher_wright(sp, 60, 5.0; mut_base=1e5, popsize=sizes, rng=Xoshiro(4))[1] == haps

    # Invalid trajectories.
    run(p; nt=10) = fisher_wright(20, nt, chr, 5.0; popsize=p)
    @test_throws ArgumentError run(fill(20, 9))
    @test_throws ArgumentError run(fill(20, 11))
    @test_throws ArgumentError run([fill(20, 9); 1])
    @test_throws ArgumentError run([0; fill(20, 9)])
    @test_throws ArgumentError run(g -> 20.0)
    @test_throws ArgumentError run(g -> g == 4 ? 1 : 20)
    @test_throws ArgumentError run(20)
    @test_throws "generation 4" run(g -> g == 4 ? 1 : 20)
end

@testset "Mating with fixed sex numbers" begin
    # Constructor.
    @test RandomMating() == RandomMating(nmale=nothing, sexratio=nothing)
    @test_throws ArgumentError RandomMating(nmale=3, sexratio=0.5)
    @test_throws ArgumentError RandomMating(nmale=0)
    @test_throws ArgumentError RandomMating(sexratio=0.0)
    @test_throws ArgumentError RandomMating(sexratio=1.0)
    @test _nmale(RandomMating(), 20) == 0                  # binomial sex
    @test _nmale(RandomMating(nmale=3), 20) == 3
    @test [_nmale(RandomMating(sexratio=0.5), n) for n in (2, 3, 20, 21)] == [1, 2, 10, 10]
    @test [_nmale(RandomMating(sexratio=0.1), n) for n in (2, 14, 16, 100)] == [1, 1, 2, 10]
    @test _nmale(RandomMating(sexratio=0.99), 5) == 4      # one female is kept

    # Kernel: exactly nm males, sires are males and dams females, and rows
    # beyond n2 are left alone.
    rng = Xoshiro(1)
    n, nm, n2 = 12, 3, 40
    pm, sex, sires, dams = fill(-1, 50, 2), Vector{Bool}(undef, n), Int[], Int[]
    _fixed_sex_mate!(pm, sex, sires, dams, nm, n2, rng)
    @test sum(sex) == nm && sires == findall(sex) && dams == findall(!, sex)
    @test all(sex[pm[1:n2, 1]]) && !any(sex[pm[1:n2, 2]])
    @test all(==(-1), pm[n2+1:end, :])
    for (n, nm) in ((2, 1), (5, 4), (5, 1))   # the extremes
        s = Vector{Bool}(undef, n)
        _fixed_sex_mate!(pm, s, sires, dams, nm, 10, rng)
        @test sum(s) == nm
    end
    @test_throws ErrorException _fixed_sex_mate!(pm, sex, sires, dams, 0, n2, rng)
    @test_throws ErrorException _fixed_sex_mate!(pm, sex, sires, dams, n, n2, rng)
    @test_throws ErrorException _fixed_sex_mate!(pm, sex, sires, dams, nm, 51, rng)

    # The males are a uniform random subset: over R draws each of the 12
    # individuals is male with probability 1/4 and each pair with probability
    # 3·2/(12·11) = 1/22, and within a draw the 3 sires and the 9 dams are
    # used uniformly. Pearson χ² against the uniform expectation; the three
    # statistics are checked at the 3.5σ-equivalent quantile of their
    # (approximate) χ² laws, Wilson–Hilferty.
    R = 20_000
    male, pairs, siruse, damuse = zeros(Int, n), zeros(Int, n, n), zeros(Int, nm), zeros(Int, n - nm)
    for _ = 1:R
        _fixed_sex_mate!(pm, sex, sires, dams, nm, n2, rng)
        male .+= sex
        for a in sires, b in sires
            a < b && (pairs[a, b] += 1)
        end
        for i = 1:n2
            siruse[searchsortedfirst(sires, pm[i, 1])] += 1
            damuse[searchsortedfirst(dams, pm[i, 2])] += 1
        end
    end
    χ²(obs, e) = sum(o -> (o - e)^2 / e, obs)
    crit(df) = df * (1 - 2 / 9df + 3.5 * sqrt(2 / 9df))^3
    pr = [pairs[a, b] for a = 1:n for b = a+1:n]
    @test χ²(pr, R / 22) < crit(length(pr) - 1)
    # Marginal counts are negatively correlated beyond the multinomial case
    # (each draw adds exactly 3), which only lowers the statistic.
    @test χ²(male, R / 4) < crit(n - 1)
    @test χ²(siruse, R * n2 / nm) < crit(nm - 1)
    @test χ²(damuse, R * n2 / (n - nm)) < crit(n - nm - 1)

    # Through fisher_wright: valid output, reproducible, both return modes
    # agree, and the default scheme given explicitly is the default run.
    chr, kw = [100_000, 200_000], (M=5e4, mut_base=1e5)
    valid(h) = issorted(h) && allunique(h) && all(p -> 1 <= p <= 300_000, h)
    a = fisher_wright(20, 60, chr, 5.0; kw..., rng=Xoshiro(3))
    @test fisher_wright(20, 60, chr, 5.0; kw..., mating=RandomMating(), rng=Xoshiro(3)) == a
    for m in (RandomMating(nmale=1), RandomMating(nmale=19), RandomMating(sexratio=0.5), RandomMating(sexratio=0.2))
        res = fisher_wright(20, 60, chr, 5.0; kw..., mating=m, result=true, rng=Xoshiro(3))
        haps, _ = fisher_wright(20, 60, chr, 5.0; kw..., mating=m, rng=Xoshiro(3))
        @test length(haps) == 40 && all(valid, haps) && _with_fixed(res) == haps
        @test haps != a[1]
    end
    # One sire: every paternal haplotype descends from one individual per
    # generation, so fixation is much faster than under the default.
    one = fisher_wright(20, 60, chr, 5.0; kw..., mating=RandomMating(nmale=1), result=true, rng=Xoshiro(3))
    dft = fisher_wright(20, 60, chr, 5.0; kw..., result=true, rng=Xoshiro(3))
    @test length(one.substitutions) > length(dft.substitutions)

    # With a changing census size: nmale must fit every parental cohort, the
    # last cohort is not a parent, and sexratio follows the size.
    sizes = [fill(4, 5); fill(30, 4); 3]
    run(m; ne=20, p=sizes) = fisher_wright(ne, 10, chr, 5.0; kw..., popsize=p, mating=m, rng=Xoshiro(8))
    @test length(run(RandomMating(nmale=3))[1]) == 6
    @test length(run(RandomMating(sexratio=0.9))[1]) == 6
    @test_throws ArgumentError run(RandomMating(nmale=4))
    @test_throws "parents of generation 2" run(RandomMating(nmale=4))
    @test_throws "parents of generation 1" run(RandomMating(nmale=3); ne=3)
    @test_throws ArgumentError fisher_wright(20, 10, chr, 5.0; mating=RandomMating(nmale=20))
    @test_throws TypeError fisher_wright(20, 10, chr, 5.0; mating=:random)

    # Species method.
    sp = GenericSpecies("Sex", Int32(20), UInt32.(chr), UInt32(50_000))
    m = RandomMating(nmale=4)
    @test fisher_wright(sp, 60, 5.0; mut_base=1e5, mating=m, rng=Xoshiro(4)) ==
          fisher_wright(20, 60, chr, 5.0; kw..., mating=m, rng=Xoshiro(4))
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

    function captured(redirect, f)
        pipe = Pipe()
        Base.link_pipe!(pipe; reader_supports_async=true, writer_supports_async=true)
        reader = @async read(pipe, String)
        try
            redirect(f, pipe)
        finally
            close(pipe.in)
        end
        return fetch(reader)
    end

    @test isempty(captured(redirect_stdout, () -> fisher_wright(10, 120, chr, 1.0)))
    @test isempty(captured(redirect_stderr, () -> fisher_wright(10, 120, chr, 1.0)))
    @test isempty(captured(redirect_stdout,
        () -> redirect_stderr(devnull) do
            fisher_wright(10, 120, chr, 1.0; verbose=true)
        end))
    log = captured(redirect_stderr,
        () -> fisher_wright(10, 120, chr, 1.0; verbose=true, progress_interval=50,
                            result=true))
    # Non-TTY stderr: one line per report at g = 50, 100 and the final g = 120.
    @test count("Generation", log) == 3
    @test occursin("Generation 120 / 120", log)
    @test occursin("substitutions", log)
    @test occursin("ETA", log)
    @test_throws ArgumentError fisher_wright(10, 5, chr, 1.0; progress_interval=0)
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

@testset "Direct chip extraction and subsetting" begin
    muts = [
        UInt32[10, 30, 50],
        UInt32[20, 30],
        UInt32[10, 50],
        UInt32[20, 40, 50],
    ]
    cbp = UInt32[50]
    res = FisherWrightResult(muts, cbp, UInt32[])

    # Extract specific 3 chip markers: [10, 30, 50]
    chip_positions = UInt32[10, 30, 50]
    chip_bit = extract_chip_bitarray(muts, chip_positions)
    @test size(chip_bit) == (3, 4)
    # Marker 1 (pos 10): present in hap 1 and hap 3
    @test chip_bit[1, :] == [true, false, true, false]
    # Marker 2 (pos 30): present in hap 1 and hap 2
    @test chip_bit[2, :] == [true, true, false, false]
    # Marker 3 (pos 50): present in hap 1, 3, 4
    @test chip_bit[3, :] == [true, false, true, true]
    @test_throws ArgumentError extract_chip_bitarray(muts, UInt32[30, 10])
    @test_throws ArgumentError extract_chip_bitarray(muts, UInt32[10, 10])

    # Test via to_haplotype with chip positions
    hap_chip = to_haplotype(res, chip_positions)
    @test hap_chip isa BnGStructs.Haplotype
    @test size(hap_chip) == (3, 4)
    @test hap_chip[1, 1] == true
    @test hap_chip[1, 2] == false

    # Test via to_haplotype with LocusSet
    lset = BnGStructs.LocusSet("TestChip", [1, 3]) # indices 1 and 3 into [10, 30, 50] -> pos 10, 50
    hap_lset = to_haplotype(res, lset, chip_positions)
    @test size(hap_lset) == (2, 4)
    @test hap_lset[1, :] == [true, false, true, false]
    @test hap_lset[2, :] == [true, false, true, true]
end

@testset "Thread safety and BitMatrix chunk alignment in export" begin
    # Test non-64-aligned locus count (nlc % 64 != 0) where adjacent columns share UInt64 chunks
    nlc = 65
    nhp = 16 * max(4, Threads.nthreads())
    all_mts = collect(UInt32(1):UInt32(nlc))
    muts = [copy(all_mts) for _ in 1:nhp]
    cbp = [UInt32(nlc)]

    # Serial reference extraction
    ref_xy = falses(nlc, nhp)
    for i in 1:nhp, m in muts[i]
        ref_xy[m, i] = true
    end

    # Stress test across multiple repetitions to ensure no lost bits from data races
    for _ in 1:500
        xy_chip = extract_chip_bitarray(muts, all_mts)
        @test xy_chip == ref_xy
        @test count(xy_chip) == nlc * nhp

        xy_muts, _ = muts2bitarray(muts, cbp; flip=false, include_fixed=true)
        @test xy_muts == ref_xy
        @test count(xy_muts) == nlc * nhp
    end

    # Test random sparse mutations across threads
    rng = MersenneTwister(1234)
    for _ in 1:50
        rand_muts = [sort(unique(rand(rng, UInt32(1):UInt32(nlc), rand(rng, 5:25)))) for _ in 1:nhp]
        ref_sparse = falses(nlc, nhp)
        for i in 1:nhp, m in rand_muts[i]
            ref_sparse[m, i] = true
        end

        xy_chip = extract_chip_bitarray(rand_muts, all_mts)
        @test xy_chip == ref_sparse

        xy_muts, _ = muts2bitarray(rand_muts, cbp; flip=false, include_fixed=true)
        @test xy_muts == ref_sparse
    end
end

# BEGIN population genetics validation (included in docs/src/manual/validation.md)

"""
Pooled ``σ²_d = ΣD² / Σ p₁q₁p₂q₂`` for SNP pairs on the same chromosome, per
distance bin, using SNPs with minor allele frequency ≥ `maf`. Returns the
numerator and denominator sums so replicates can be pooled.
"""
function ld_sums(bm, chr, pos, bins; maf = 0.1)
    n = size(bm, 2)
    p = vec(sum(bm; dims = 2)) ./ n
    num, den = zeros(length(bins)), zeros(length(bins))
    for c in unique(chr)
        k = findall(i -> chr[i] == c && maf <= p[i] <= 1 - maf, eachindex(p))
        G = Float64.(bm[k, :])
        q = p[k]
        D = G * G' ./ n .- q * q'
        V = (q .* (1 .- q)) * (q .* (1 .- q))'
        d = Int.(pos[k])' .- Int.(pos[k])
        for (j, (lo, hi)) in enumerate(bins)
            m = (lo .<= d) .& (d .< hi)
            num[j] += sum(abs2, D[m])
            den[j] += sum(V[m])
        end
    end
    return num, den
end

@testset "Neutral Wright-Fisher theory" begin
    # One set of equilibrium runs (10×2N generations) feeds every check below.
    # Each bound was set from replicate runs: it passes for the correct model
    # with several seeds, and fails if the simulated Ne is halved or the
    # recombination rate doubled (and, for S, on the pre-v0.3.5 loop order).
    Random.seed!(2026)
    ne, chr, reps, nsub = 100, fill(10_000_000, 10), 20, 20
    θ = 4ne * 1e-8 * sum(chr)
    groups = [1:1, 2:2, 3:4, 5:9, 10:19]          # derived-allele counts in a sample of nsub
    bins = [(300_000, 1_000_000), (1_000_000, 3_000_000)]
    S = H = 0.0
    ξ = zeros(nsub - 1)
    num, den = zeros(length(bins)), zeros(length(bins))
    for _ = 1:reps
        res = fisher_wright(ne, 20ne, chr, 1.0; result = true)
        haps = res.active_haplotypes
        n = length(haps)
        cnt = Dict{UInt32,Int}()
        for h in haps, p in h
            cnt[p] = get(cnt, p, 0) + 1
        end
        for c in values(cnt)
            0 < c < n || continue
            S += 1
            H += 2 * (c / n) * (1 - c / n)
        end
        # Site frequency spectrum of a random sample of nsub ≪ 2N haplotypes,
        # where the coalescent expectation E[ξₖ] = θ/k applies.
        sub = Dict{UInt32,Int}()
        for h in haps[randperm(n)[1:nsub]], p in h
            sub[p] = get(sub, p, 0) + 1
        end
        for c in values(sub)
            c < nsub && (ξ[c] += 1)
        end
        bm, lmp = muts2bitarray(haps, res.chromosome_ends)
        a, b = ld_sums(bm, lmp.chr, lmp.pos, bins)
        num .+= a
        den .+= b
    end

    @testset "Segregating sites and heterozygosity" begin
        # S should match Watterson's θ·aₙ and Σ2pq should match θ = 4Nμ·L. The
        # old mutate-then-mate order gave S ≈ 0.92·θ·aₙ. SE is about 0.7% for S
        # and 1.3% for Σ2pq.
        aₙ = sum(1 / i for i = 1:2ne-1)
        @test 0.97 < S / reps / (θ * aₙ) < 1.06
        @test 0.90 < H / reps / θ < 1.10
    end

    @testset "Site frequency spectrum" begin
        # Observed/expected for each group of classes; SE is 2-3% per group.
        # Halving the simulated Ne gives about 0.5 in every group.
        for g in groups
            @test 0.90 < sum(ξ[g]) / reps / (θ * sum(1 / k for k in g)) < 1.10
        end
    end

    @testset "Linkage disequilibrium decay" begin
        # Reference σ²_d from msprime DTWF, same settings (N = 100, whole
        # population, MAF ≥ 0.1, μ = r = 1e-8, 100 × 10 chromosomes of 10 Mb):
        # 0.2869 ± 0.0024 at 0.3-1 Mb and 0.1398 ± 0.0015 at 1-3 Mb, from
        # bench/validate-popgen-msprime.py --reps 100 --seed 1 --summary. SE here is
        # about 2-3%. Doubling the recombination rate gives ratios of about 0.66
        # and 0.55; halving Ne gives about 1.4 and 1.6. Sved's 1/(1 + 4Nc) is not
        # used: it overestimates r² at short distances and ignores the MAF filter.
        ref = [0.2869, 0.1398]
        for j in eachindex(bins)
            @test 0.85 < num[j] / den[j] / ref[j] < 1.15
        end
    end
end

# Expected Σ2pq / L over all haplotypes of each generation, starting from
# independent alleles at frequency ½ in `n0` individuals; `sizes[g]` is the
# number of offspring of generation g. Exact for the model: sire and dam are
# drawn with replacement, and the parents' sex is Binomial(N, ½) conditioned
# on both sexes or, with `nmale`, exactly `nmale(N)` males among N parents.
# hB and hW are the heterozygosities of two genes in different individuals
# and in one individual (docs/src/manual/validation.md, "Effective size").
# The sex counts are independent of the genotypes and across generations, so
# the expectation uses a = E[1/(4Nm) + 1/(4Nf)] over the N parents of each
# generation, whatever the number of offspring.
function expected_heterozygosity(n0, sizes; nmale=nothing)
    function coalescence(N)
        if nmale !== nothing
            m = nmale(N)
            return 1 / 4m + 1 / 4(N - m)
        end
        w = [binomial(big(N), m) for m = 1:N-1]
        return Float64(sum(w[m] * (1 // (4m) + 1 // (4(N - m))) for m = 1:N-1) / sum(w))
    end
    hB = hW = 0.5
    H = Float64[]
    N = n0
    for n in sizes
        a = coalescence(N)
        hB, hW = (1 - a) * hB + a / 2 * hW, hB
        push!(H, hW / 2n + (1 - 1 / n) * hB)
        N = n
    end
    return H
end

# Mean Σ2pq / L over the haplotypes `haps`, for sites 1:L.
function site_heterozygosity(haps, L)
    n = length(haps)
    c = zeros(Int, L)
    for h in haps, p in h
        c[p] += 1
    end
    return sum(x -> 2 * (x / n) * (1 - x / n), c) / L
end

@testset "Heterozygosity decay without mutation" begin
    # Drift alone, from a polymorphic start given by init: Σ2pq must follow
    # the exact recursion of the two-sex model, whose leading eigenvalue is
    # 1 - 1/(2Ne) with Ne = 19.38 for N = 20. 400 replicates of 1,000 sites on
    # 20 chromosomes of 1 Morgan; the second leg continues the first through
    # a FisherWrightResult. z-scores against the replicate standard error;
    # over four seeds they ranged from -1.63 to 2.53.
    N, chr, reps, t₁, t₂ = 20, fill(50, 20), 400, 10, 30
    L = sum(chr)
    rng = Xoshiro(20261110)
    het(haps) = site_heterozygosity(haps, L)
    H₁, H₂ = zeros(reps), zeros(reps)
    for r = 1:reps
        start = [UInt32.(findall(rand(rng, Bool, L))) for _ = 1:2N]
        mid = fisher_wright(N, t₁, chr, 1.0; M=50, init=start, mutate=false, result=true, rng=rng)
        last, _ = fisher_wright(N, t₂, chr, 1.0; M=50, init=mid, mutate=false, rng=rng)
        H₁[r] = het(_with_fixed(mid))
        H₂[r] = het(last)
    end
    E = expected_heterozygosity(N, fill(N, t₁ + t₂))
    z₁ = (mean(H₁) - E[t₁]) / sqrt(var(H₁) / reps)
    z₂ = (mean(H₂) - E[t₁+t₂]) / sqrt(var(H₂) / reps)
    @test abs(z₁) < 3.5
    @test abs(z₂) < 3.5
    # The test separates Ne from N: an ideal monoecious population of N = 20,
    # ½(1 - 1/2N)^(t+1) (the extra factor is the realized Σ2pq of the start),
    # is 8 to 11 standard errors away.
    ideal = 0.5 * (1 - 1 / 2N)^(t₁ + t₂ + 1)
    @test (ideal - mean(H₂)) / sqrt(var(H₂) / reps) > 3.5
end

@testset "Heterozygosity through a bottleneck and an expansion" begin
    # Drift under a changing census size. 30 individuals are reduced to 6
    # for 8 generations and then expand to 60 for 12; the second leg starts
    # from the first through init. Σ2pq must follow the exact recursion, in
    # which the drift of generation g is set by its N_{g-1} parents. 400
    # replicates of 1,000 sites on 20 chromosomes of 1 Morgan; z-scores
    # against the replicate standard error; over four seeds |z| was at most
    # 1.34, and the alternative below 12 to 15 standard errors away.
    N₀, Nb, Nx, chr, reps, t₁, t₂ = 30, 6, 60, fill(50, 20), 400, 8, 12
    L = sum(chr)
    rng = Xoshiro(20261111)
    H₁, H₂ = zeros(reps), zeros(reps)
    for r = 1:reps
        start = [UInt32.(findall(rand(rng, Bool, L))) for _ = 1:2N₀]
        mid = fisher_wright(N₀, t₁, chr, 1.0; M=50, init=start, popsize=fill(Nb, t₁),
            mutate=false, result=true, rng=rng)
        last, _ = fisher_wright(Nb, t₂, chr, 1.0; M=50, init=mid, popsize=g -> Nx,
            mutate=false, rng=rng)
        H₁[r] = site_heterozygosity(_with_fixed(mid), L)
        H₂[r] = site_heterozygosity(last, L)
    end
    sizes = [fill(Nb, t₁); fill(Nx, t₂)]
    E = expected_heterozygosity(N₀, sizes)
    se₁, se₂ = sqrt(var(H₁) / reps), sqrt(var(H₂) / reps)
    @test abs(mean(H₁) - E[t₁]) / se₁ < 3.5
    @test abs(mean(H₂) - E[end]) / se₂ < 3.5
    # The test can tell parents from offspring: if the drift of generation g
    # followed its own size N_g, the first generation would drift as 6
    # individuals instead of 30.
    wrong = expected_heterozygosity(Nb, sizes)
    @test (mean(H₁) - wrong[t₁]) / se₁ > 3.5
    @test (mean(H₂) - wrong[end]) / se₂ > 3.5
end

@testset "Heterozygosity decay under a fixed sex ratio" begin
    # Drift alone with fixed numbers of each sex, against the exact
    # recursion with a = 1/(4Nm) + 1/(4Nf). Three schemes over the same 20
    # generations, 400 replicates of 1,000 sites on 20 chromosomes of
    # 1 Morgan each: 3 sires and 17 dams among N = 20 (Ne ≈ 4·3·17/20 = 10.2);
    # equal sexes, 10 and 10 (Ne ≈ 20); and a proportion of 0.2 males through
    # a change in size, 20 → 10 for 10 generations → 40 (4, 2 and 8 males).
    # z-scores against the replicate standard error.
    chr, reps, t = fill(50, 20), 400, 20
    L = sum(chr)
    rng = Xoshiro(20261112)
    step = [fill(10, 10); fill(40, 10)]
    cases = (
        (RandomMating(nmale=3), nothing, N -> 3),
        (RandomMating(sexratio=0.5), nothing, N -> N ÷ 2),
        (RandomMating(sexratio=0.2), step, N -> N ÷ 5),
    )
    for (mating, popsize, nmale) in cases
        H = zeros(reps)
        for r = 1:reps
            start = [UInt32.(findall(rand(rng, Bool, L))) for _ = 1:40]
            haps, _ = fisher_wright(20, t, chr, 1.0; M=50, init=start, mutate=false,
                popsize=popsize, mating=mating, rng=rng)
            H[r] = site_heterozygosity(haps, L)
        end
        sizes = something(popsize, fill(20, t))
        se = sqrt(var(H) / reps)
        @test abs(mean(H) - expected_heterozygosity(20, sizes; nmale=nmale)[end]) / se < 3.5
        # The test separates the fixed scheme from binomial sex.
        @test abs(mean(H) - expected_heterozygosity(20, sizes)[end]) / se > 3.5
    end
end

@testset "Crossovers per meiosis matches Poisson and assortment" begin
    Random.seed!(42)
    # Chromosome 1: 50 Mb (expected 0.50 crossovers at M=1e8)
    # Chromosome 2: 100 Mb (expected 1.00 crossovers at M=1e8)
    chr = [50_000_000, 100_000_000]
    M = 1e8
    rmap = uniform_recombination_map(chr; M=M)
    cbuf = UInt32[]
    K = 10_000
    counts_chr1 = zeros(Int, K)
    counts_chr2 = zeros(Int, K)
    assortment = zeros(Int, K)

    for k in 1:K
        cobp!(cbuf, rmap)
        counts_chr1[k] = count(<=(50_000_000), cbuf)
        counts_chr2[k] = count(x -> 50_000_001 < x <= 150_000_000, cbuf)
        assortment[k] = count(==(50_000_001), cbuf)
    end

    # Sample means must conform to theoretical Poisson expectations
    @test abs(mean(counts_chr1) - 0.50) < 0.025
    @test abs(mean(counts_chr2) - 1.00) < 0.035

    # Dispersion index (variance/mean) for Poisson process must be ~1.0
    @test 0.90 < var(counts_chr1) / mean(counts_chr1) < 1.10
    @test 0.90 < var(counts_chr2) / mean(counts_chr2) < 1.10

    # Independent assortment boundary crossover probability is 0.50
    @test abs(mean(assortment) - 0.50) < 0.020
end

@testset "Chromosome boundary assortment segregates terminal loci" begin
    # Chr 1: 1..1000, Chr 2: 1001..2000
    # Assortment breakpoint is at 1001. Locus 1000 is on Chr 1, locus 1001 is on Chr 2.
    h1 = UInt32[1000, 1001]
    h2 = UInt32[]
    ho = UInt32[]
    cross_overs = UInt32[1001]
    # Under recombination with cross_overs=[1001], one segment is < 1001 (chr 1)
    # and the second segment is >= 1001 (chr 2).
    # The offspring receives either [1000] or [1001], never both together and never neither.
    seen_1000 = false
    seen_1001 = false
    for _ in 1:100
        recombine(h1, h2, ho, cross_overs)
        @test (ho == UInt32[1000]) || (ho == UInt32[1001])
        ho == UInt32[1000] && (seen_1000 = true)
        ho == UInt32[1001] && (seen_1001 = true)
    end
    @test seen_1000 && seen_1001
end

# END population genetics validation

include("properties.jl")
include("golden.jl")
