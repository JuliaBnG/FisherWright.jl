using Test
using Random
using FisherWright: fisher_wright, muts2bitarray, extract_chip_bitarray, merge_sorted, merge_sorted!, recombine, _recombine!, _insert_sorted!, _isearchfirst,
    cobp, cobp!, FisherWrightResult, effective_size, to_haplotype, RecombinationMap,
    uniform_recombination_map, random_mate, random_mate!, _fixed_sex_mate!, _nmale, _fixed_mutations,
    _drop_mutations, _drop_mutations!, _fixation_step, _fixation_step!, _extract_fixed,
    _intersect_sorted!, _remove_sorted!, _with_fixed, _simulate!, _initial_state,
    SimState, ConstantSize, SizeTrajectory, RandomMating, Neutral, bottleneck, growth,
    AbstractFitness, AliasTable, _alias_build!, _alias_draw, _alias_probabilities,
    _draw_parents!, _mate!, FitnessLoci, fitness_loci, SiteFitness, _log_fitness, _selected,
    _meiosis!, _meiosis_selected!, _sel_mutation!, TraitLoci, trait_loci, TraitFitness, Stabilizing,
    Directional, _genotypic_value, FamilySize, Hierarchical, Monoecious, AvoidRelatives, CustomMating,
    _record_families!, _related
using Distributions: Gamma, Exponential, Normal, Chisq, quantile
import FisherWright: _fitness!
using BnGStructs
using DataFrames: nrow, names
using Distributions: Poisson
using Statistics: mean, var, cor, std

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
    @test st.census == [ne] && isempty(st.males)

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

@testset "Trajectory helpers" begin
    # bottleneck: nb in generations start to start + duration - 1.
    b = bottleneck(100, 10, 21, 5)
    @test b.(1:30) == [fill(100, 20); fill(10, 5); fill(100, 5)]
    @test bottleneck(100, 10, 1, 1).(1:3) == [10, 100, 100]
    @test bottleneck(100, 10, 3, 2; recovery=40).(1:6) == [100, 100, 10, 10, 40, 40]
    @test bottleneck(Int32(8), 0x02, big(2), 1)(2) === 2
    @test (@inferred b(3)) === 100
    for bad in ((1, 10, 2, 2), (100, 1, 2, 2), (100, 10, 0, 2), (100, 10, 2, 0))
        @test_throws ArgumentError bottleneck(bad...)
    end
    @test_throws ArgumentError bottleneck(100, 10, 2, 2; recovery=1)

    # growth: geometric, rounded, optionally delayed and limited.
    @test growth(50, 0.05).(1:20) == [round(Int, 50 * 1.05^g) for g = 1:20]
    @test growth(50, 0.0).(1:5) == fill(50, 5)
    @test growth(50, 1; start=3).(1:6) == [50, 50, 50, 100, 200, 400]
    @test growth(50, 0.5; limit=120).(1:4) == [75, 112, 120, 120]
    @test growth(50, -0.5; limit=10).(1:4) == [25, 12, 10, 10]
    @test growth(50, -0.1; start=20, limit=10).([20, 21, 60]) == [50, 45, 10]
    @test (@inferred growth(50, 0.05; limit=500)(7)) isa Int
    @test_throws ArgumentError growth(50, 1)(100)          # 50·2^100 does not fit
    @test_throws ArgumentError growth(1, 0.1)
    @test_throws ArgumentError growth(50, -1)
    @test_throws ArgumentError growth(50, Inf)
    @test_throws ArgumentError growth(50, 0.1; start=-1)
    @test_throws ArgumentError growth(50, 0.1; limit=1)
    @test_throws ArgumentError growth(50, 0.1; limit=40)   # ceiling below the start
    @test_throws ArgumentError growth(50, -0.1; limit=60)  # floor above the start

    # As popsize: the census record is the trajectory, and the run is the
    # one given by the same sizes as a vector.
    chr, kw = [100_000, 200_000], (M=5e4, mut_base=1e5, result=true)
    for (ne, f) in ((30, bottleneck(30, 4, 6, 5)), (10, growth(10, 0.1; limit=40)),
        (30, growth(30, -0.15; start=5, limit=6)))
        res = fisher_wright(ne, 25, chr, 5.0; kw..., popsize=f, rng=Xoshiro(2))
        vec = fisher_wright(ne, 25, chr, 5.0; kw..., popsize=f.(1:25), rng=Xoshiro(2))
        @test res.census == [ne; f.(1:25)]
        @test res.active_haplotypes == vec.active_haplotypes && res.males == vec.males
    end
    # The bottleneck cohorts are the parents of the following generations.
    res = fisher_wright(30, 25, chr, 5.0; kw..., popsize=bottleneck(30, 4, 6, 5), rng=Xoshiro(2))
    @test all(<(4), res.males[7:11]) && effective_size(res, 7:11; from=:sexes) < 4 < effective_size(res, 1:6; from=:sexes)
    @test effective_size(res, 7:11) < effective_size(res, 1:6)
    # A decline without a floor is stopped by fisher_wright at the first
    # generation whose size is below 2: 20·0.5^4 = 1.25.
    @test_throws "generation 4" fisher_wright(20, 10, chr, 5.0; popsize=growth(20, -0.5))
end

@testset "Census record and effective size" begin
    chr, kw = [100_000, 200_000], (M=5e4, mut_base=1e5, result=true)

    # Constant size, default mating: N_0, …, N_nt and one sex count per
    # generation, with both sexes present.
    res = fisher_wright(20, 60, chr, 5.0; kw..., rng=Xoshiro(3))
    @test res.census == fill(20, 61) && length(res.males) == 60
    @test all(m -> 1 <= m <= 19, res.males) && length(unique(res.males)) > 1
    # The record does not change the realization or depend on the scan.
    @test _with_fixed(res) == fisher_wright(20, 60, chr, 5.0; M=5e4, mut_base=1e5, rng=Xoshiro(3))[1]
    scan = fisher_wright(20, 60, chr, 5.0; kw..., fixation_interval=7, rng=Xoshiro(3))
    @test (scan.census, scan.males) == (res.census, res.males)

    # Binomial sex: the male count of N = 50 parents has mean 25 and variance
    # 12.5 (truncation at 0 and 50 is below 1e-14). z-test of the mean over
    # 4,000 independent generations, and the variance within 10%.
    long = fisher_wright(50, 4000, [1_000], 1.0; result=true, mutate=false, rng=Xoshiro(12))
    @test abs(mean(long.males) - 25) / sqrt(12.5 / 4000) < 3.5
    @test 0.9 < var(long.males) / 12.5 < 1.1
    # Ne from the record: Jensen's inequality puts it below 4·25·25/50 = 50,
    # at E[1/4Nm + 1/4Nf]⁻¹ = 48.96 for N = 50.
    @test 48.7 < effective_size(long; from=:sexes) < 49.3

    # Fixed sexes and a changing size: the record is the request.
    sizes = [fill(10, 5); fill(40, 4); 7]
    fix = fisher_wright(20, 10, chr, 5.0; kw..., popsize=sizes, mating=RandomMating(nmale=3), rng=Xoshiro(4))
    @test fix.census == [20; sizes] && fix.males == fill(3, 10)
    @test length(fix.active_haplotypes) == 2fix.census[end]
    rat = fisher_wright(20, 10, chr, 5.0; kw..., popsize=sizes, mating=RandomMating(sexratio=0.2), rng=Xoshiro(4))
    @test rat.males == [4; fill(2, 5); fill(8, 4)]   # parents: 20, 10 × 5, 40 × 4

    # effective_size(; from = :sexes): Wright's formula per generation, harmonic mean over
    # several, parents only (the last cohort of 7 does not enter).
    ne(m, n) = 4m * (n - m) / n
    @test effective_size(fix, 1; from=:sexes) ≈ ne(3, 20) ≈ 10.2
    @test effective_size(fix, 2; from=:sexes) ≈ ne(3, 10)
    @test effective_size(fix, 7:10; from=:sexes) ≈ ne(3, 40)
    @test effective_size(fix; from=:sexes) ≈ 10 / (1 / ne(3, 20) + 5 / ne(3, 10) + 4 / ne(3, 40))
    @test effective_size(fix, [1, 10]; from=:sexes) ≈ 2 / (1 / ne(3, 20) + 1 / ne(3, 40))
    @test effective_size(rat, 2:6; from=:sexes) ≈ ne(2, 10) && effective_size(rat, 10; from=:sexes) ≈ ne(8, 40)
    @test effective_size(fisher_wright(20, 5, chr, 5.0; kw..., mating=RandomMating(sexratio=0.5)); from=:sexes) ≈ 20
    @test (@inferred effective_size(fix)) isa Float64
    @test (@inferred effective_size(fix; from=:sexes)) isa Float64
    @test_throws BoundsError effective_size(fix, 11)
    @test_throws BoundsError effective_size(fix, 0)
    @test_throws ArgumentError effective_size(fix, 1:0)
    @test_throws ArgumentError effective_size(fix; from=:pedigree)

    # The realized family sizes: the recorded pairs are those of the mating
    # table, and under random mating their expectation is the sex-count
    # formula, so the two agree over a long run.
    @test length(fix.paternal_pairs) == length(fix.maternal_pairs) == 10
    rng = Xoshiro(13)
    rmap = uniform_recombination_map(chr; M=5e4)
    st = SimState([UInt32[] for _ = 1:40], UInt32[], 20, rng)
    for _ = 1:20
        _simulate!(st, ConstantSize(20), RandomMating(), Neutral(), 1, rmap, Poisson(1.0), UInt32(1):UInt32(300_000);
            mutate=false, track=false, fixation_interval=1, verbose=false, progress_interval=100, rng=rng)
        pairs(col) = (c = [count(==(k), st.pm[1:20, col]) for k = 1:20]; sum(c .* (c .- 1)))
        @test (st.ppairs[end], st.mpairs[end]) == (pairs(1), pairs(2))
    end
    @test abs(effective_size(long) / effective_size(long; from=:sexes) - 1) < 0.02

    # A run continued from a result appends to its record; from plain
    # haplotypes, or in a result built by hand, the record starts anew.
    cont = fisher_wright(7, 3, chr, 5.0; kw..., init=fix, popsize=[9, 9, 12], rng=Xoshiro(5))
    @test cont.census == [20; sizes; 9; 9; 12] && cont.males[1:10] == fix.males
    @test length(cont.males) == 13 && all(m -> 1 <= m <= 8, cont.males[12:13])
    @test effective_size(cont, 1:10) == effective_size(fix)
    @test cont.paternal_pairs[1:10] == fix.paternal_pairs && length(cont.maternal_pairs) == 13
    @test fix.census == [20; sizes] && fix.males == fill(3, 10)   # init left untouched
    plain = fisher_wright(7, 3, chr, 5.0; kw..., init=_with_fixed(fix), rng=Xoshiro(5))
    @test plain.census == fill(7, 4) && length(plain.males) == 3
    bare = FisherWrightResult(fix.active_haplotypes, fix.chromosome_ends, fix.substitutions)
    @test fisher_wright(7, 3, chr, 5.0; kw..., init=bare, rng=Xoshiro(5)).census == fill(7, 4)
    # A record that does not end at the individuals of init is rejected.
    bad(c, m) = FisherWrightResult(fix.active_haplotypes, fix.chromosome_ends, fix.substitutions, c, m)
    @test_throws ArgumentError fisher_wright(7, 3, chr, 5.0; init=bad([20, 8], [3]))
    @test_throws ArgumentError fisher_wright(7, 3, chr, 5.0; init=bad([20, 7], [3, 3]))

    # A hand-built result has no record.
    manual = FisherWrightResult([UInt32[1, 3], UInt32[2], UInt32[], UInt32[1, 2, 3]], UInt32[5], UInt32[5])
    @test manual.census == [2] && isempty(manual.males)
    @test_throws ArgumentError effective_size(manual)
    five = FisherWrightResult(manual.active_haplotypes, UInt32[5], UInt32[5], [2, 2], [1])
    @test all(isnan, five.paternal_pairs) && effective_size(five; from=:sexes) ≈ 2   # 1 male, 1 female
    @test_throws "not recorded" effective_size(five)
    kept = _extract_fixed(FisherWrightResult(copy.(fix.active_haplotypes), fix.chromosome_ends,
        fix.substitutions, fix.census, fix.males))
    @test (kept.census, kept.males) == (fix.census, fix.males)

    # Species method.
    sp = GenericSpecies("Census", Int32(20), UInt32.(chr), UInt32(50_000))
    @test fisher_wright(sp, 10, 5.0; mut_base=1e5, result=true, popsize=sizes).census == [20; sizes]
end

# Fitness by position in the cohort, for tests of weighted mating: not
# heritable, so it changes the variance of family size and nothing else.
struct IndexFitness <: AbstractFitness
    w::Vector{Float64}
end
_fitness!(f::IndexFitness, st::SimState) = (copyto!(st.w, 1, f.w, 1, st.n); nothing)

# (s, h, u, SE) from docs/src/validation/fixation-reference.txt, N = 30.
const FIXATION_REFERENCE = [
    (0.2, 0.5, 0.160166, 0.000259),
    (0.2, 0.0, 0.055675, 0.000162),
    (0.2, 1.0, 0.278524, 0.000317),
    (0.0, 0.5, 0.016693, 0.000091),
]

const χ²crit(df) = df * (1 - 2 / 9df + 3.5 * sqrt(2 / 9df))^3   # 3.5σ-equivalent, Wilson–Hilferty

@testset "Alias sampler" begin
    rng = Xoshiro(41)
    t = AliasTable(8)
    # The table is exact: the probabilities it implies are the normalised
    # weights, for weights of very different size, zeros and a subset.
    for (w, idx) in (
        ([1.0, 2.0, 3.0, 4.0], [1, 2, 3, 4]),
        ([5.0], [1]),
        ([1e-9, 1.0, 1e9, 0.0, 3.0], [1, 2, 3, 4, 5]),
        ([0.0, 0.0, 7.0], [1, 2, 3]),
        (fill(0.3, 7), collect(1:7)),
        (rand(rng, 200) .^ 8, collect(1:200)),
        ([9.0, 1.0, 9.0, 2.0, 9.0, 5.0], [2, 4, 6]),     # outcomes 1:3 ↔ individuals 2, 4, 6
    )
        _alias_build!(t, w, idx)
        @test length(t) == length(idx)
        @test _alias_probabilities(t) ≈ w[idx] ./ sum(w[idx]) rtol = 1e-12
        @test all(p -> 0 <= p <= 1, t.prob) && all(in(1:length(idx)), t.alias)
        @test all(_alias_draw(t, rng) in 1:length(idx) for _ = 1:100)
    end
    # Draws follow the weights (Pearson χ²), and a zero weight is never drawn.
    w = [4.0, 0.0, 1.0, 2.5, 0.5, 2.0]
    _alias_build!(t, w, collect(1:6))
    R = 200_000
    c = zeros(Int, 6)
    for _ = 1:R
        c[_alias_draw(t, rng)] += 1
    end
    @test c[2] == 0
    keep = w .> 0
    @test sum((c[keep] .- R .* w[keep] ./ 10) .^ 2 ./ (R .* w[keep] ./ 10)) < χ²crit(4)
    # Rebuilding within the reserved size does not allocate, nor does a draw.
    idx = collect(1:6)
    rebuild(t, w, idx) = @allocated _alias_build!(t, w, idx)
    draw(t, rng) = @allocated _alias_draw(t, rng)
    rebuild(t, w, idx); draw(t, rng)
    @test rebuild(t, w, idx) == 0 && draw(t, rng) == 0
    # Invalid weights.
    for bad in ([1.0, -0.1], [1.0, NaN], [Inf, 1.0], [0.0, 0.0])
        @test_throws ArgumentError _alias_build!(t, bad, [1, 2])
    end
    @test_throws ArgumentError _alias_build!(t, [1.0], Int[])
end

@testset "Weighted parents" begin
    rng = Xoshiro(42)
    # Parents are used in proportion to their weight within each sex.
    w = [1.0, 5.0, 2.0, 0.0, 3.0, 1.0, 4.0, 2.0, 0.5, 1.5]
    sires, dams = [2, 3, 4, 9], [1, 5, 6, 7, 8, 10]
    pm, tm, tf = zeros(Int, 60, 2), AliasTable(), AliasTable()
    R = 5_000
    use = zeros(Int, 10, 2)
    for _ = 1:R
        _draw_parents!(pm, sires, dams, 50, rng, w, tm, tf)
        for i = 1:50, s = 1:2
            use[pm[i, s], s] += 1
        end
    end
    @test all(==(0), pm[51:60, :])                       # rows beyond n2 untouched
    @test sum(use[dams, 1]) == sum(use[sires, 2]) == 0   # sires are males, dams females
    @test use[4, 1] == 0                                 # zero fitness: never a parent
    χ²(obs, p) = sum((obs .- 50R .* p) .^ 2 ./ (50R .* p))
    ps, pd = w[[2, 3, 9]] ./ sum(w[sires]), w[dams] ./ sum(w[dams])
    @test χ²(use[[2, 3, 9], 1], ps) < χ²crit(2)
    @test χ²(use[dams, 2], pd) < χ²crit(5)
    # A sex whose members all have zero fitness cannot reproduce.
    @test_throws ArgumentError _draw_parents!(pm, [4], dams, 50, rng, w, tm, tf)

    # Through the engine: the state holds one weight per individual of the
    # largest cohort, Neutral leaves it alone, and a fitness model fills it.
    chr, cbp = [1_000, 1_000], UInt32[1_000, 2_000]
    rmap = uniform_recombination_map(chr; M=1e3)
    run!(st, f, nt, r; demog=ConstantSize(st.n)) = _simulate!(
        st, demog, RandomMating(nmale=4), f, nt, rmap, Poisson(1.0), UInt32(1):UInt32(2_000);
        mutate=true, track=false, fixation_interval=1, verbose=false, progress_interval=100, rng=r)
    st = SimState(_initial_state(nothing, 20, cbp, false)..., 10, rng)
    @test st.w == ones(10)
    run!(st, Neutral(), 5, rng)
    @test st.w == ones(10)
    f = IndexFitness(w)
    run!(st, f, 5, rng)
    @test st.w == w && length(st.tm) == 4 && length(st.tf) == 6
    @test (@inferred run!(st, f, 1, rng)) === st
    run!(st, f, 50, rng)
    @test @allocated(run!(st, f, 50, rng)) / 50 < 2^16
    # Individual 4 has no offspring: with four males of which it may be one,
    # the weighted draw must still never pick it.
    for _ = 1:200
        _mate!(RandomMating(nmale=4), f, st, 10, rng)
        @test !(4 in view(st.pm, 1:10, :))
    end
    # Equal fitness is random mating in distribution (not the same draws).
    flat = IndexFitness(ones(10))
    _mate!(RandomMating(), flat, st, 10, rng)
    @test all(st.sex[st.pm[1:10, 1]]) && !any(st.sex[st.pm[1:10, 2]])
end

@testset "Fitness loci" begin
    # Explicit table: scalars are expanded, vectors kept.
    l = FitnessLoci([10, 20, 30], 0.1)
    @test l.pos == UInt32[10, 20, 30] && l.s == fill(0.1, 3) && l.h == fill(0.5, 3) && l.μ == 0.0
    @test length(l) == 3
    l = FitnessLoci(UInt32[5, 9], [-1.0, 0.3]; h=[0.0, 2.0], μ=1e-4)
    @test l.s == [-1.0, 0.3] && l.h == [0.0, 2.0] && l.μ == 1e-4
    @test length(FitnessLoci(Int[], 0.1)) == 0
    for bad in (
        () -> FitnessLoci([20, 10], 0.1),            # not sorted
        () -> FitnessLoci([10, 10], 0.1),            # not unique
        () -> FitnessLoci([0, 10], 0.1),             # not 1-based
        () -> FitnessLoci([10, 20], [0.1]),          # lengths
        () -> FitnessLoci([10], -1.5),               # s < -1
        () -> FitnessLoci([10], NaN),
        () -> FitnessLoci([10], -0.5; h=3.0),        # 1 + hs < 0
        () -> FitnessLoci([10], 0.1; μ=-1e-3),
    )
        @test_throws ArgumentError bad()
    end

    # Drawn table: distinct sites within the genome, effects from the DFE,
    # reproducible from the rng.
    chr = [1_000, 2_000]
    l = fitness_loci(chr, 200; s=-Gamma(0.5, 0.1) + 0.0, h=x -> 0.25, μ=1e-3, rng=Xoshiro(1))
    @test length(l) == 200 && issorted(l.pos) && allunique(l.pos) && all(p -> 1 <= p <= 3_000, l.pos)
    @test all(<(0), l.s) && all(==(0.25), l.h) && l.μ == 1e-3
    l2 = fitness_loci(chr, 200; s=-Gamma(0.5, 0.1) + 0.0, h=x -> 0.25, μ=1e-3, rng=Xoshiro(1))
    @test (l2.pos, l2.s) == (l.pos, l.s)
    @test fitness_loci(chr, 3_000; s=0.01).pos == 1:3_000        # every site
    @test length(fitness_loci(chr, 2_500; s=0.01).pos) == 2_500   # the dense branch
    @test length(fitness_loci(chr, 0; s=0.01)) == 0
    @test_throws ArgumentError fitness_loci(chr, 3_001; s=0.01)
    @test_throws ArgumentError fitness_loci(chr, 10; s="strong")
    @test_throws ArgumentError fitness_loci([0], 1; s=0.1)

    # Fitness of an individual from the merge-walk equals the product over
    # sites of the dense genotype: random tables and haplotypes.
    rng = Xoshiro(2)
    for _ = 1:300
        K = rand(rng, 1:40)
        pos = sort!(unique(rand(rng, UInt32(1):UInt32(500), K)))
        K = length(pos)
        l = FitnessLoci(pos, rand(rng, K) .* 1.5 .- 0.5; h=rand(rng, K) .* 1.2)
        f = SiteFitness(l)
        a, b = pos[rand(rng, Bool, K)], pos[rand(rng, Bool, K)]
        dense = prod(1:K; init=1.0) do k
            d = (pos[k] in a) + (pos[k] in b)
            d == 0 ? 1.0 : d == 1 ? 1 + l.h[k] * l.s[k] : 1 + l.s[k]
        end
        @test exp(_log_fitness(f, a, b)) ≈ dense
        @test _log_fitness(f, a, b) == _log_fitness(f, b, a)
    end
    f = SiteFitness(FitnessLoci([7, 9], [-1.0, 0.5]; h=[0.5, 0.5]))
    @test _log_fitness(f, UInt32[], UInt32[]) == 0.0
    @test _log_fitness(f, UInt32[7], UInt32[7]) == -Inf        # lethal homozygote
    @test exp(_log_fitness(f, UInt32[7, 9], UInt32[9])) ≈ 0.5 * 1.5
    alloc(f, a, b) = @allocated _log_fitness(f, a, b)
    alloc(f, UInt32[7], UInt32[9])
    @test alloc(f, UInt32[7], UInt32[9]) == 0
    @test _selected(f) && !_selected(Neutral()) && !_selected(IndexFitness([1.0]))
end

@testset "Selected stream" begin
    chr, cbp = [5_000, 5_000], UInt32[5_000, 10_000]
    rmap = uniform_recombination_map(chr; M=2_000)
    N = 20
    loci = FitnessLoci(collect(50:50:10_000), 0.05; μ=2e-3)
    f = SiteFitness(loci)
    rng = Xoshiro(3)
    rand_sel() = loci.pos[rand(rng, length(loci)) .< 0.2]
    rand_neu() = UInt32.(filter(p -> p % 50 != 0, sort!(unique(rand(rng, 1:10_000, 60)))))

    # Meiosis on both streams: the neutral stream is what `_meiosis!` gives
    # for the same matings and streams, and the two streams together are the
    # recombinant of their union under the same crossovers and strand.
    function state(seed)
        st = SimState([rand_neu() for _ = 1:2N], UInt32[], N, Xoshiro(seed), [rand_sel() for _ = 1:2N], UInt32[])
        st.pm[:, 1] .= rand(Xoshiro(seed + 1), 1:N, N)
        st.pm[:, 2] .= rand(Xoshiro(seed + 2), 1:N, N)
        return st
    end
    rng = Xoshiro(3); a = state(5)
    rng = Xoshiro(3); b = state(5)
    parents = [merge_sorted(a.prt[k], a.sel[k]) for k = 1:2N]
    streams = copy.(b.srng)
    _meiosis!(a, rmap, N)
    _meiosis_selected!(b, rmap, N)
    @test a.prt[1:2N] == b.prt[1:2N]
    @test a.srng == b.srng                        # same draws from every stream
    co, out = UInt32[], UInt32[]
    for i = 1:N, (p, k) in ((b.pm[i, 1], 2i - 1), (b.pm[i, 2], 2i))
        r = streams[i]
        cobp!(co, rmap; rng=r)
        _recombine!(parents[2p-1], parents[2p], out, co, rand(r, Bool))
        @test merge_sorted(b.prt[k], b.sel[k]) == out
        @test issubset(b.sel[k], loci.pos) && isdisjoint(b.prt[k], loci.pos)
    end

    # Mutation at the fitness loci: only sites of the table, sorted and
    # unique, about K·μ = 0.4 new alleles per haplotype per generation (fewer
    # by the hits on sites already carried), and none at substituted sites.
    st = SimState([UInt32[] for _ = 1:2N], UInt32[], N, Xoshiro(8))
    p_sel = Poisson(length(loci) * loci.μ)
    gained = 0
    for _ = 1:200
        before = sum(length, st.sel)
        _sel_mutation!(f, st, p_sel)
        gained += sum(length, st.sel) - before
        foreach(empty!, st.sel)
    end
    @test abs(gained / (200 * 2N) - 0.4) < 3.5 * sqrt(0.4 / (200 * 2N))
    st.ssubs = loci.pos[1:100]
    for _ = 1:50
        _sel_mutation!(f, st, p_sel)
    end
    @test all(h -> issorted(h) && allunique(h) && issubset(h, loci.pos[101:end]), st.sel[1:2N])
    @test sum(length, st.sel) > 100

    # Through fisher_wright: both streams valid, kept apart, reproducible,
    # and independent of the fixation interval.
    kw = (M=2_000, mut_base=1e4, result=true, fitness=f)
    res = fisher_wright(N, 300, chr, 2.0; kw..., rng=Xoshiro(4))
    @test length(res.selected_haplotypes) == length(res.active_haplotypes) == 2N
    @test all(h -> issorted(h) && allunique(h) && issubset(h, loci.pos), res.selected_haplotypes)
    @test issorted(res.selected_substitutions) && issubset(res.selected_substitutions, loci.pos)
    @test !isempty(res.selected_substitutions) && any(!isempty, res.selected_haplotypes)
    @test all(h -> isdisjoint(h, res.selected_substitutions), res.selected_haplotypes)
    @test isempty(_fixed_mutations(res.selected_haplotypes))
    again = fisher_wright(N, 300, chr, 2.0; kw..., fixation_interval=13, rng=Xoshiro(4))
    @test again.selected_haplotypes == res.selected_haplotypes
    @test again.selected_substitutions == res.selected_substitutions
    @test again.active_haplotypes == res.active_haplotypes && again.substitutions == res.substitutions
    # Advantageous alleles accumulate: more selected substitutions than
    # under s = 0 at the same sites and mutation rate (110 against 72 here;
    # over 20 seeds 117 ± 7 against 74 ± 7, mean ± SD).
    flat = SiteFitness(FitnessLoci(loci.pos, 0.0; μ=loci.μ))
    neu = fisher_wright(N, 300, chr, 2.0; kw..., fitness=flat, rng=Xoshiro(4))
    @test length(res.selected_substitutions) > length(neu.selected_substitutions) + 20

    # A neutral run has an empty selected stream and is the run without the
    # keyword.
    plain = fisher_wright(N, 50, chr, 2.0; M=2_000, mut_base=1e4, result=true, rng=Xoshiro(6))
    @test all(isempty, plain.selected_haplotypes) && isempty(plain.selected_substitutions)
    @test length(plain.selected_haplotypes) == 2N
    same = fisher_wright(N, 50, chr, 2.0; M=2_000, mut_base=1e4, result=true, fitness=Neutral(), rng=Xoshiro(6))
    @test same.active_haplotypes == plain.active_haplotypes

    # Continuing: the selected alleles of init are carried on, with their
    # substitutions; mutate = false adds none; init is left untouched.
    keep = deepcopy(res)
    cont = fisher_wright(N, 40, chr, 2.0; kw..., init=res, mutate=false, rng=Xoshiro(7))
    @test issubset(res.selected_substitutions, cont.selected_substitutions)
    @test issubset(reduce(union, cont.selected_haplotypes), reduce(union, res.selected_haplotypes))
    @test res.selected_haplotypes == keep.selected_haplotypes
    @test length(cont.census) == 341
    # A starting population with one copy of a selected allele.
    sel = [UInt32[] for _ = 1:2N]
    push!(sel[3], loci.pos[7])
    start = FisherWrightResult([UInt32[] for _ = 1:2N], cbp, UInt32[]; selected_haplotypes=sel)
    one = fisher_wright(N, 1, chr, 2.0; kw..., init=start, mutate=false, rng=Xoshiro(9))
    @test all(h -> h == UInt32[] || h == [loci.pos[7]], one.selected_haplotypes)
    @test isempty(sel[1]) && sel[3] == [loci.pos[7]]

    # Selection adds variance in family size: the realized effective size is
    # below the sex-count value, which ignores it.
    strong = SiteFitness(FitnessLoci(loci.pos, 0.5; μ=2e-3))
    sres = fisher_wright(N, 200, chr, 2.0; kw..., fitness=strong, rng=Xoshiro(12))
    @test effective_size(sres, 50:200) < 0.9 * effective_size(sres, 50:200; from=:sexes)

    # Errors.
    @test_throws ArgumentError fisher_wright(N, 5, chr, 2.0; fitness=f)                 # result = false
    @test_throws ArgumentError fisher_wright(N, 5, chr, 2.0; init=res, result=true)     # selected alleles, no model
    @test_throws ArgumentError fisher_wright(N, 5, [5_000], 2.0; fitness=f, result=true) # loci beyond the genome
    other = SiteFitness(FitnessLoci([1, 2, 3], 0.1))
    @test_throws ArgumentError fisher_wright(N, 5, chr, 2.0; fitness=other, init=res, result=true)
    short = FisherWrightResult([UInt32[] for _ = 1:2N], cbp, UInt32[]; selected_haplotypes=sel[1:4])
    @test_throws ArgumentError fisher_wright(N, 5, chr, 2.0; kw..., init=short)
    unsorted = FisherWrightResult([UInt32[] for _ = 1:2N], cbp, UInt32[];
        selected_haplotypes=[[[loci.pos[2], loci.pos[1]]]; sel[2:end]])
    @test_throws ArgumentError fisher_wright(N, 5, chr, 2.0; kw..., init=unsorted)

    # Species method, and a lethal recessive never becomes homozygous... in
    # an adult that reproduces: no selected substitution can occur.
    sp = GenericSpecies("Sel", Int32(N), UInt32.(chr), UInt32(2_000))
    @test fisher_wright(sp, 300, 2.0; mut_base=1e4, result=true, fitness=f, rng=Xoshiro(4)).selected_haplotypes ==
          res.selected_haplotypes
    lethal = SiteFitness(FitnessLoci(collect(50:50:10_000), -1.0; h=0.0, μ=1e-4))
    dead = fisher_wright(N, 200, chr, 2.0; kw..., fitness=lethal, mating=RandomMating(sexratio=0.5),
        rng=Xoshiro(10))
    @test isempty(dead.selected_substitutions) && any(!isempty, dead.selected_haplotypes)
    # A sex in which nobody has positive fitness ends the run with an error.
    doomed = FisherWrightResult([UInt32[] for _ = 1:2N], cbp, UInt32[];
        selected_haplotypes=[UInt32[50] for _ = 1:2N])
    @test_throws "cannot reproduce" fisher_wright(N, 5, chr, 2.0; kw..., fitness=lethal, init=doomed)
end

@testset "Export of selected alleles" begin
    chr, cbp, N = [5_000, 5_000], UInt32[5_000, 10_000], 20
    loci = FitnessLoci(collect(50:50:10_000), [0.02k / 200 for k = 1:200]; h=0.3, μ=2e-3)
    res = fisher_wright(N, 300, chr, 2.0; M=2_000, mut_base=1e4, result=true, fitness=SiteFitness(loci),
        rng=Xoshiro(4))
    sel, ssubs = res.selected_haplotypes, res.selected_substitutions
    @test !isempty(ssubs) && any(!isempty, sel) && !isempty(res.substitutions)

    # The fitness loci with their effects: one row per site of the table,
    # in its order, against a brute-force lookup.
    qtl, map = to_haplotype(res, loci)
    @test qtl isa Haplotype && size(qtl) == (200, 2N)
    @test propertynames(map) == [:chr, :pos, :s, :h, :frq]
    @test map.pos == loci.pos && map.s == loci.s && map.h == loci.h
    @test map.chr == [p <= 5_000 ? 1 : 2 for p in loci.pos] && eltype(map.chr) == Int8
    dense = [(loci.pos[k] in sel[j]) || (loci.pos[k] in ssubs) for k = 1:200, j = 1:2N]
    @test qtl == dense
    @test map.frq ≈ vec(sum(dense; dims=2)) ./ 2N
    @test all(==(1.0), map.frq[indexin(ssubs, loci.pos)])       # substituted sites are fixed
    # Without the substitutions they read as absent.
    qtl0, map0 = to_haplotype(res, loci; include_fixed=false)
    @test qtl0 == [loci.pos[k] in sel[j] for k = 1:200, j = 1:2N]
    @test all(==(0.0), map0.frq[indexin(ssubs, loci.pos)])
    # The fitness that the run used is recovered from the export.
    f = SiteFitness(loci)
    for i in (1, 7, N)
        w = prod(1:200) do k
            d = qtl[k, 2i-1] + qtl[k, 2i]
            d == 0 ? 1.0 : d == 1 ? 1 + loci.h[k] * loci.s[k] : 1 + loci.s[k]
        end
        fixed = prod(1 + loci.s[k] for k in indexin(ssubs, loci.pos))
        @test w ≈ exp(_log_fitness(f, sel[2i-1], sel[2i])) * fixed
    end

    # Neutral and selected alleles in one matrix.
    hap, lmap = to_haplotype(res)
    both, bmap = to_haplotype(res; selected=true)
    @test !isdisjoint(lmap.pos, loci.pos)       # neutral mutations also fall on coordinates of fitness loci
    segsel = [p for p in loci.pos if 0 < count(h -> p in h, sel) < 2N]
    @test sort(bmap.pos) == sort!(union(lmap.pos, segsel))
    for (r, p) in enumerate(bmap.pos)
        @test collect(both[r, :]) == [(p in res.active_haplotypes[j]) || (p in sel[j]) for j = 1:2N]
    end
    allb, amap = to_haplotype(res; selected=true, include_fixed=true)
    @test issubset(ssubs, amap.pos) && issubset(res.substitutions, amap.pos)
    @test all(==(1.0f0), amap.frq[indexin(ssubs, amap.pos)])
    # Chip extraction: a selected site is carried only with selected = true.
    chip = UInt32[loci.pos[1], segsel[1], ssubs[1], 9_999]
    sort!(unique!(chip))
    plain = to_haplotype(res, chip)
    withsel = to_haplotype(res, chip; selected=true, include_fixed=true)
    k = findfirst(==(segsel[1]), chip)
    @test !any(plain[k, :]) && collect(withsel[k, :]) == [segsel[1] in h for h in sel]
    @test all(withsel[findfirst(==(ssubs[1]), chip), :])
    @test to_haplotype(res, chip; selected=true)[findfirst(==(ssubs[1]), chip), :] == falses(2N)

    # A coordinate carried in both streams is one locus.
    clash = FisherWrightResult([UInt32[100], UInt32[], UInt32[300], UInt32[]], UInt32[1_000], UInt32[];
        selected_haplotypes=[UInt32[], UInt32[100], UInt32[], UInt32[]])
    m, cm = to_haplotype(clash; selected=true)
    @test cm.pos == UInt32[100, 300] && collect(m[1, :]) == [true, true, false, false]

    # A neutral result exports as before, and has nothing at fitness loci.
    neutral = fisher_wright(N, 50, chr, 2.0; M=2_000, mut_base=1e4, result=true, rng=Xoshiro(6))
    @test to_haplotype(neutral; selected=true)[2].pos == to_haplotype(neutral)[2].pos
    q, mq = to_haplotype(neutral, loci)
    @test !any(q) && all(==(0.0), mq.frq)

    # Errors.
    @test_throws ArgumentError to_haplotype(res, FitnessLoci(Int[], 0.1))
    @test_throws ArgumentError to_haplotype(res, FitnessLoci([20_000], 0.1))      # beyond the genome
    @test_throws ArgumentError to_haplotype(res, FitnessLoci([50, 100], 0.1))     # not the run's table
end

@testset "Fitness record" begin
    chr, cbp, N = [5_000, 5_000], UInt32[5_000, 10_000], 20
    kw = (M=2_000, mut_base=1e4, result=true)
    record(r) = (r.mean_fitness, r.mean_log_fitness, r.var_log_fitness)

    # A neutral run records fitness 1 for the parents of every generation.
    neutral = fisher_wright(N, 30, chr, 2.0; kw..., rng=Xoshiro(1))
    @test record(neutral) == (ones(30), zeros(30), zeros(30))

    # With selection, entry g is the fitness of the parents of generation
    # g: one generation from a known population, against the dense product.
    loci = FitnessLoci(collect(50:50:10_000), [0.03 * (-1)^k for k = 1:200]; h=0.3, μ=2e-3)
    f = SiteFitness(loci)
    res = fisher_wright(N, 200, chr, 2.0; kw..., fitness=f, rng=Xoshiro(2))
    @test all(v -> length(v) == length(res.males) == 200, record(res))
    @test !isempty(res.selected_substitutions)
    @test all(>(0), res.var_log_fitness[50:end])
    qtl, _ = to_haplotype(res, loci)          # substituted sites included
    lw = [sum(1:200) do k
        d = qtl[k, 2i-1] + qtl[k, 2i]
        d == 0 ? 0.0 : d == 1 ? log1p(loci.h[k] * loci.s[k]) : log1p(loci.s[k])
    end for i = 1:N]
    next = fisher_wright(N, 1, chr, 2.0; kw..., fitness=f, init=res, rng=Xoshiro(3))
    @test next.mean_log_fitness[end] ≈ mean(lw)
    @test next.var_log_fitness[end] ≈ var(lw)
    @test next.mean_fitness[end] ≈ mean(exp, lw)
    # A continued run appends to the record of init, which is left alone.
    @test all(v -> length(v) == 201, record(next))
    @test next.mean_fitness[1:200] == res.mean_fitness && length(res.mean_fitness) == 200
    # The record does not depend on the fixation interval, although the
    # sites move between the haplotypes and the substitutions.
    again = fisher_wright(N, 200, chr, 2.0; kw..., fitness=f, fixation_interval=17, rng=Xoshiro(2))
    @test again.mean_fitness ≈ res.mean_fitness && again.mean_log_fitness ≈ res.mean_log_fitness
    @test again.var_log_fitness ≈ res.var_log_fitness atol = 1e-12

    # A site carried by everyone: fitness 1 + s for all, before and after it
    # is moved to the substitutions.
    one = SiteFitness(FitnessLoci([100], 0.25))
    full = FisherWrightResult([UInt32[] for _ = 1:2N], cbp, UInt32[];
        selected_haplotypes=[UInt32[100] for _ = 1:2N])
    fixed = fisher_wright(N, 3, chr, 2.0; kw..., fitness=one, init=full, mutate=false)
    @test fixed.selected_substitutions == UInt32[100] && all(isempty, fixed.selected_haplotypes)
    @test fixed.mean_fitness ≈ fill(1.25, 3) && fixed.mean_log_fitness ≈ fill(log(1.25), 3)
    @test fixed.var_log_fitness == zeros(3)

    # Directional selection raises mean fitness; purifying selection keeps
    # it below 1.
    up = fisher_wright(N, 300, chr, 2.0; kw..., fitness=SiteFitness(FitnessLoci(loci.pos, 0.05; μ=2e-3)),
        rng=Xoshiro(4))
    @test up.mean_fitness[1] == 1.0 && up.mean_fitness[end] > 2 && issorted(up.mean_fitness[[1, 100, 200, 300]])
    down = fisher_wright(N, 300, chr, 2.0; kw..., fitness=SiteFitness(FitnessLoci(loci.pos, -0.05; μ=2e-3)),
        rng=Xoshiro(4))
    @test all(<=(1), down.mean_fitness) && down.mean_fitness[end] < 0.9

    # A lethal homozygote among the parents: mean log fitness is -Inf.
    lethal = SiteFitness(FitnessLoci([100], -1.0; h=0.0))
    mixed = FisherWrightResult([UInt32[] for _ = 1:2N], cbp, UInt32[];
        selected_haplotypes=[i <= 2 ? UInt32[100] : UInt32[] for i = 1:2N])
    dead = fisher_wright(N, 1, chr, 2.0; kw..., fitness=lethal, init=mixed, mutate=false,
        mating=RandomMating(sexratio=0.5))
    @test dead.mean_fitness == [(N - 1) / N] && dead.mean_log_fitness == [-Inf] && isnan(dead.var_log_fitness[1])

    # Results built by hand: a neutral record of the right length, and a
    # record of the wrong length is rejected as a start.
    bare = FisherWrightResult(res.active_haplotypes, cbp, res.substitutions)
    @test record(bare) == (Float64[], Float64[], Float64[])
    five = FisherWrightResult(res.active_haplotypes, cbp, res.substitutions, fill(N, 4), fill(10, 3))
    @test record(five) == (ones(3), zeros(3), zeros(3))
    wrong = FisherWrightResult(res.active_haplotypes, cbp, res.substitutions, fill(N, 4), fill(10, 3),
        res.selected_haplotypes, res.selected_substitutions, ones(2), zeros(3), zeros(3))
    @test_throws ArgumentError fisher_wright(N, 1, chr, 2.0; kw..., fitness=f, init=wrong)
    kept = _extract_fixed(res)
    @test record(kept) == record(res)
end

@testset "Trait loci and genotypic value" begin
    l = TraitLoci([10, 20, 30], [-2.0, 0.5, 3.0]; h=[0.5, 0.0, 1.0], μ=1e-4)
    @test l.pos == UInt32[10, 20, 30] && l.a == [-2.0, 0.5, 3.0] && l.h == [0.5, 0.0, 1.0] && l.μ == 1e-4
    @test TraitLoci([5], 1.0).h == [0.5] && length(TraitLoci(Int[], 1.0)) == 0
    for bad in (() -> TraitLoci([20, 10], 1.0), () -> TraitLoci([0], 1.0), () -> TraitLoci([1, 2], [1.0]),
        () -> TraitLoci([1], Inf), () -> TraitLoci([1], 1.0; h=NaN), () -> TraitLoci([1], 1.0; μ=-1.0))
        @test_throws ArgumentError bad()
    end
    t = trait_loci([1_000, 2_000], 100; a=Normal(0, 1), h=x -> x > 0 ? 0.2 : 0.8, μ=1e-3, rng=Xoshiro(1))
    @test length(t) == 100 && issorted(t.pos) && allunique(t.pos) && all(<=(3_000), t.pos)
    @test any(<(0), t.a) && any(>(0), t.a) && t.h == [x > 0 ? 0.2 : 0.8 for x in t.a]
    @test trait_loci([1_000, 2_000], 100; a=Normal(0, 1), h=x -> x > 0 ? 0.2 : 0.8, μ=1e-3, rng=Xoshiro(1)).a == t.a
    @test_throws ArgumentError trait_loci([1_000], 5; a="big")
    for bad in (() -> Stabilizing(0, 0), () -> Stabilizing(Inf, 1), () -> Directional(NaN),
        () -> TraitFitness(t; σe=-1, fitness=Directional(0.1)))
        @test_throws ArgumentError bad()
    end

    # The merge-walk equals the sum over the dense genotype.
    rng = Xoshiro(2)
    for _ = 1:300
        K = rand(rng, 1:40)
        pos = sort!(unique(rand(rng, UInt32(1):UInt32(500), K)))
        K = length(pos)
        l = TraitLoci(pos, randn(rng, K); h=rand(rng, K) .* 1.5 .- 0.25)
        x, y = pos[rand(rng, Bool, K)], pos[rand(rng, Bool, K)]
        dense = sum(1:K; init=0.0) do k
            d = (pos[k] in x) + (pos[k] in y)
            d == 0 ? 0.0 : d == 1 ? l.h[k] * l.a[k] : l.a[k]
        end
        @test _genotypic_value(l, x, y) ≈ dense
        @test _genotypic_value(l, x, y) ≈ _genotypic_value(l, y, x)
    end
    alloc(l, x, y) = @allocated _genotypic_value(l, x, y)
    l = TraitLoci([7, 9], [1.0, 2.0])
    alloc(l, UInt32[7], UInt32[9])
    @test alloc(l, UInt32[7], UInt32[9]) == 0
end

@testset "Selection on a trait" begin
    chr, cbp = [5_000, 5_000], UInt32[5_000, 10_000]
    kw = (M=2_000, mut_base=1e4, result=true)
    N = 30
    loci = TraitLoci(collect(100:100:10_000), [0.2 * (-1)^k for k = 1:100]; μ=2e-3)

    # The phenotype is G + e: with no derived alleles, G = 0 and the
    # variance of P is σe². 40 generations × 30 parents, χ² on 39·40 df.
    f = TraitFitness(TraitLoci(loci.pos, loci.a; μ=0.0); σe=2.0, fitness=Directional(0.3))
    r = fisher_wright(N, 40, chr, 2.0; kw..., fitness=f, rng=Xoshiro(1))
    @test r.mean_genotypic == zeros(40) && r.var_genotypic == zeros(40)
    df = 40 * (N - 1)
    @test abs(sum(r.var_phenotypic) * (N - 1) / 4 - df) / sqrt(2df) < 3.5
    @test abs(mean(r.mean_phenotypic) / sqrt(4 / (40N))) < 3.5
    @test r.mean_log_fitness ≈ 0.3 .* r.mean_phenotypic
    @test r.var_log_fitness ≈ 0.09 .* r.var_phenotypic

    # One generation of selection is unbiased: for additive loci, the
    # expected mean G of the offspring, given the parents and the fitness of
    # each, is the fitness-weighted mean G of the parents within each sex,
    # averaged over the sexes. Repeated from the same parents, with new
    # environmental deviations, matings and meioses each time.
    rmap = uniform_recombination_map(chr; M=2_000)
    rng = Xoshiro(3)
    lt = TraitLoci(collect(100:100:10_000), randn(rng, 100) .* 0.3)
    start = [lt.pos[rand(rng, 100) .< 0.3] for _ = 1:2N]
    for fit in (Directional(0.8), Stabilizing(1.0, 0.5))
        tf = TraitFitness(lt; σe=0.7, fitness=fit)
        R = 3_000
        d, shift = zeros(R), zeros(R)
        for r = 1:R
            st = SimState([UInt32[] for _ = 1:2N], UInt32[], N, rng, copy.(start), UInt32[])
            _simulate!(st, ConstantSize(N), RandomMating(), tf, 1, rmap, Poisson(1.0), UInt32(1):UInt32(10_000);
                mutate=false, track=false, fixation_interval=1, verbose=false, progress_interval=100, rng=rng)
            G, w, sex = st.gval[1:N], st.w[1:N], st.sex[1:N]
            if r == 1
                # Fitness is the function of the phenotype, up to the
                # scaling of the fittest to 1.
                lw = [fit isa Directional ? fit.β * p : -(p - fit.θ)^2 / (2fit.ω²) for p in st.pval[1:N]]
                @test log.(w) .- lw ≈ fill(-maximum(lw), N) atol = 1e-9
                @test maximum(w) == 1.0
            end
            expect = (sum(G[sex] .* w[sex]) / sum(w[sex]) + sum(G[.!sex] .* w[.!sex]) / sum(w[.!sex])) / 2
            got = mean(_genotypic_value(lt, st.sel[2i-1], st.sel[2i]) for i = 1:N)
            d[r] = got - expect
            shift[r] = got - mean(G)
        end
        @test abs(mean(d)) / (std(d) / sqrt(R)) < 3.5
        # ... and selection does act: the offspring differ from the mean of
        # their parents.
        @test abs(mean(shift)) / (std(shift) / sqrt(R)) > 10
    end

    # Directional selection moves the mean in its direction; stabilizing
    # selection holds it near the optimum.
    up = fisher_wright(N, 150, chr, 2.0; kw..., fitness=TraitFitness(loci; σe=0.5, fitness=Directional(1.0)),
        rng=Xoshiro(4))
    down = fisher_wright(N, 150, chr, 2.0; kw..., fitness=TraitFitness(loci; σe=0.5, fitness=Directional(-1.0)),
        rng=Xoshiro(4))
    @test up.mean_genotypic[end] > 2 && down.mean_genotypic[end] < -2
    held = fisher_wright(N, 300, chr, 2.0; kw..., fitness=TraitFitness(loci; σe=0.5, fitness=Stabilizing(1.5, 0.5)),
        rng=Xoshiro(5))
    @test abs(mean(held.mean_genotypic[201:300]) - 1.5) < 0.5
    @test !isempty(up.selected_substitutions) && !isempty(held.selected_substitutions)

    # Substituted sites keep their value: under stabilizing selection, which
    # is not multiplicative over sites, the realization does not depend on
    # how often fixed sites are moved out of the haplotypes.
    ff = TraitFitness(loci; σe=0.5, fitness=Stabilizing(1.5, 0.5))
    a1 = fisher_wright(N, 300, chr, 2.0; kw..., fitness=ff, rng=Xoshiro(5))
    a2 = fisher_wright(N, 300, chr, 2.0; kw..., fitness=ff, fixation_interval=17, rng=Xoshiro(5))
    @test a1.selected_haplotypes == held.selected_haplotypes == a2.selected_haplotypes
    @test a1.selected_substitutions == a2.selected_substitutions
    @test a1.mean_genotypic ≈ a2.mean_genotypic && a1.mean_phenotypic ≈ a2.mean_phenotypic

    # The export gives the genotypic values: one more generation without
    # noise records the mean G of the exported population.
    hap, map = to_haplotype(held, loci)
    @test propertynames(map) == [:chr, :pos, :a, :h, :frq] && map.a == loci.a
    G = [sum(k -> (d = hap[k, 2i-1] + hap[k, 2i]; d == 0 ? 0.0 : d == 1 ? loci.h[k] * loci.a[k] : loci.a[k]),
        1:length(loci)) for i = 1:N]
    nxt = fisher_wright(N, 1, chr, 2.0; kw..., init=held, fitness=TraitFitness(loci; fitness=Stabilizing(1.5, 0.5)),
        rng=Xoshiro(6))
    @test nxt.mean_genotypic[end] ≈ mean(G) && nxt.var_genotypic[end] ≈ var(G)
    @test nxt.mean_phenotypic[end] == nxt.mean_genotypic[end]          # σe = 0
    @test length(nxt.mean_genotypic) == 301 && nxt.mean_genotypic[1:300] == held.mean_genotypic

    # Without a trait model the trait record is NaN, and continues as such.
    neutral = fisher_wright(N, 5, chr, 2.0; kw..., rng=Xoshiro(7))
    @test all(isnan, neutral.mean_genotypic) && length(neutral.var_phenotypic) == 5
    site = fisher_wright(N, 5, chr, 2.0; kw..., fitness=SiteFitness(FitnessLoci([50], 0.1; μ=1e-3)), rng=Xoshiro(7))
    @test all(isnan, site.mean_phenotypic) && length(site.mean_phenotypic) == 5
    bare = FisherWrightResult(neutral.active_haplotypes, cbp, UInt32[])
    @test isempty(bare.mean_genotypic)

    # A user function of the phenotype: truncation on P (the best half by
    # a threshold) gives the same as the equivalent built-in where they
    # coincide, and an invalid value is an error.
    tr = TraitFitness(loci; σe=0.5, fitness=P -> P > 0 ? 1.0 : 0.0)
    @test fisher_wright(N, 50, chr, 2.0; kw..., fitness=tr, rng=Xoshiro(8)).mean_genotypic[end] > 0.5
    flat = TraitFitness(loci; σe=0.5, fitness=P -> 1.0)
    dirz = TraitFitness(loci; σe=0.5, fitness=Directional(0.0))
    @test fisher_wright(N, 20, chr, 2.0; kw..., fitness=flat, rng=Xoshiro(9)).selected_haplotypes ==
          fisher_wright(N, 20, chr, 2.0; kw..., fitness=dirz, rng=Xoshiro(9)).selected_haplotypes
    @test_throws "non-negative" fisher_wright(N, 5, chr, 2.0; kw..., fitness=TraitFitness(loci; fitness=P -> -1.0))
    @test_throws ArgumentError fisher_wright(N, 5, [5_000], 2.0; kw..., fitness=TraitFitness(loci; fitness=Directional(1.0)))
    @test_throws ArgumentError fisher_wright(N, 5, chr, 2.0; fitness=TraitFitness(loci; fitness=Directional(1.0)))
    other = TraitFitness(TraitLoci([1, 2], 1.0); fitness=Directional(1.0))
    @test_throws ArgumentError fisher_wright(N, 5, chr, 2.0; kw..., init=held, fitness=other)

    # The environmental deviations do not depend on the number of threads.
    other_t = Threads.nthreads() == 1 ? 3 : 1
    code = "using FisherWright, Random; l = TraitLoci(collect(100:100:10_000), 0.2; μ = 2e-3); " *
           "r = fisher_wright(30, 60, [5_000, 5_000], 2.0; M = 2_000, mut_base = 1e4, result = true, " *
           "fitness = TraitFitness(l; σe = 1.0, fitness = Stabilizing(1.0, 1.0)), rng = Xoshiro(11)); " *
           "print(hash((r.selected_haplotypes, r.mean_phenotypic)))"
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -t $other_t -e $code`
    l = TraitLoci(collect(100:100:10_000), 0.2; μ=2e-3)
    here = fisher_wright(30, 60, [5_000, 5_000], 2.0; M=2_000, mut_base=1e4, result=true,
        fitness=TraitFitness(l; σe=1.0, fitness=Stabilizing(1.0, 1.0)), rng=Xoshiro(11))
    @test parse(UInt, read(cmd, String)) == hash((here.selected_haplotypes, here.mean_phenotypic))
end

# Probability that two distinct items of `total` fall to the same one of `n`
# parents under the splits of FamilySize: E[Σ c(c-1)] / (total(total-1)).
function same_parent(mode, k, total, n)
    total < 2 && return 0.0
    if mode === :equal
        q, r = divrem(total, n)
        return (r * (q + 1) * q + (n - r) * q * (q - 1)) / (total * (total - 1))
    end
    return mode === :random ? 1 / n : (k + 1) / (n * k + 1)
end

# Coalescence term a = (P_f + P_m)/4 of FamilySize for N parents with Nm
# males and `no` offspring.
function family_a(m::FamilySize, N, Nm, no)
    Nf = N - Nm
    # Breeders: random subsets of the given sizes, if any.
    m.nsires > 0 && (Nm = m.nsires)
    m.ndams > 0 && (Nf = m.ndams)
    Pf = same_parent(m.dams, m.kd, no, Nf)
    π = same_parent(m.sires, m.km, Nf, Nm)
    return (Pf + (Pf + (1 - Pf) * π)) / 4
end

@testset "Mating with family sizes" begin
    @test_throws ArgumentError FamilySize(dams=:poisson)
    @test_throws ArgumentError FamilySize(sires=0)
    @test_throws ArgumentError FamilySize(dams=Inf)
    @test_throws ArgumentError FamilySize(nmale=2, sexratio=0.5)
    m = FamilySize(dams=0.5, sires=:equal, nmale=4)
    @test (m.dams, m.kd, m.sires, m.km, m.sex.nmale) == (:gamma, 0.5, :equal, 0.0, 4)

    # The mating table: every dam has one sire, and families add up.
    rng = Xoshiro(1)
    N = 40
    st = SimState([UInt32[] for _ = 1:2N], UInt32[], 60, rng)
    families(st, no) = (pm = st.pm[1:no, :]; [count(==(d), pm[:, 2]) for d in st.dams])
    for m in (FamilySize(dams=:equal, sires=:equal, nmale=7), FamilySize(dams=0.3, sires=0.8),
        FamilySize(dams=:random, sires=:equal, sexratio=0.5)), no in (40, 57, 13)
        _mate!(m, Neutral(), st, no, rng)
        pm = st.pm[1:no, :]
        @test all(st.sex[pm[:, 1]]) && !any(st.sex[pm[:, 2]])
        @test all(length(unique(pm[pm[:, 2] .== d, 1])) <= 1 for d in st.dams)   # one sire per dam
        nd = families(st, no)
        @test sum(nd) == no
        m.dams === :equal && @test maximum(nd) - minimum(nd) <= 1
        if m.sires === :equal
            dps = [length(unique(pm[pm[:, 1] .== s, 2])) for s in st.sires]
            # dams with offspring per sire; with all dams assigned, at most
            # one apart (dams without offspring are not visible here)
            m.dams === :equal && no >= length(st.dams) && @test maximum(dps) - minimum(dps) <= 1
        end
    end
    @test_throws ArgumentError _mate!(FamilySize(dams=:equal), SiteFitness(FitnessLoci([5], 0.1)), st, 10, rng)

    # Moments of dam family sizes: binomial under :random, and
    # Dirichlet-multinomial under Gamma(k): Var = μ(1 - 1/Nf)(no + Nf k)/(1 + Nf k).
    for (m, k) in ((FamilySize(dams=:random, nmale=10), Inf), (FamilySize(dams=0.5, nmale=10), 0.5))
        R, no, Nf = 4_000, 60, 30
        v = Float64[]
        for _ = 1:R
            _mate!(m, Neutral(), st, no, rng)
            append!(v, families(st, no))
        end
        μ = no / Nf
        expect = isinf(k) ? no * (1 / Nf) * (1 - 1 / Nf) : μ * (1 - 1 / Nf) * (no + Nf * k) / (1 + Nf * k)
        # Mean square about the known mean no/Nf (the counts sum to no), per
        # replicate; its expectation is Var(n_d). SE over replicates.
        per = [mean(abs2, v[(r-1)*Nf+1:r*Nf] .- μ) for r = 1:R]
        @test abs(mean(per) - expect) / (std(per) / sqrt(R)) < 3.5
    end

    # Within the reserved capacity, a generation allocates task overhead only.
    # Without mutation, so that growing haplotypes do not reallocate their
    # buffers (whose growth depends on drift, and so on the scheme).
    chr, cbp = [1_000, 1_000], UInt32[1_000, 2_000]
    rmap = uniform_recombination_map(chr; M=1e3)
    run!(st, m, nt) = _simulate!(st, ConstantSize(40), m, Neutral(), nt, rmap, Poisson(1.0), UInt32(1):UInt32(2_000);
        mutate=false, track=false, fixation_interval=1, verbose=false, progress_interval=100, rng=rng)
    st = SimState([UInt32.(findall(rand(rng, Bool, 2_000))) for _ = 1:80], UInt32[], 40, rng)
    for m in (FamilySize(dams=:equal, sires=:equal, nmale=8), FamilySize(dams=0.5, sires=2.0))
        run!(st, m, 50)
        @test @allocated(run!(st, m, 50)) / 50 < 2^14
    end

    # Breeder subsets: exactly nsires sires and ndams dams breed, each sire
    # with dams_per_sire dams, equal families, and the sexes of all parents
    # are kept for the record.
    @test_throws ArgumentError FamilySize(nsires=0)
    @test_throws ArgumentError Hierarchical(0, 3)
    h = Hierarchical(4, 3; sexratio=0.5)
    @test (h.dams, h.sires, h.nsires, h.ndams, h.sex.sexratio) == (:equal, :equal, 4, 12, 0.5)
    st = SimState([UInt32[] for _ = 1:2N], UInt32[], 60, rng)
    for m in (Hierarchical(4, 3), Hierarchical(4, 3; nmale=10), FamilySize(dams=0.7, sires=:random, nsires=5, ndams=9)),
        no in (40, 50)
        _mate!(m, Neutral(), st, no, rng)
        pm = st.pm[1:no, :]
        used_s, used_d = unique(pm[:, 1]), unique(pm[:, 2])
        @test length(used_s) <= m.nsires && length(used_d) <= m.ndams
        @test issubset(used_s, st.sires) && issubset(used_d, st.dams)
        @test all(length(unique(pm[pm[:, 2] .== d, 1])) == 1 for d in used_d)
        if m.dams === :equal
            @test length(used_s) == 4 && length(used_d) == 12
            nd = [count(==(d), pm[:, 2]) for d in used_d]
            @test maximum(nd) - minimum(nd) <= 1
            @test all(length(unique(pm[pm[:, 1] .== s, 2])) == 3 for s in used_s)
        end
        @test length(st.sires) + length(st.dams) == N          # every parent has a sex
        m.sex.nmale > 0 && @test length(st.sires) == 10
    end
    # The sires are a uniform subset: 2 of the 10 males per generation, and
    # the males are a uniform subset of the 40 parents, so every parent is a
    # sire with probability 1/20 per generation. Pearson χ² over the 40
    # parents, 5,000 generations.
    m = Hierarchical(2, 3; nmale=10)
    use = zeros(Int, N)
    for _ = 1:5_000
        _mate!(m, Neutral(), st, 40, rng)
        for s in unique(st.pm[1:40, 1])
            use[s] += 1
        end
    end
    @test sum(use) == 2 * 5_000
    e = sum(use) / N
    @test sum((use .- e) .^ 2 ./ e) < N - 1 + 3.5 * sqrt(2(N - 1))

    # effective_size from the realized family sizes is exact for equal
    # splits, which are deterministic: Hierarchical(4, 5) in 100, with any
    # sexes.
    r = fisher_wright(100, 30, [100_000], 1.0; result=true, mating=Hierarchical(4, 5), rng=Xoshiro(4))
    @test effective_size(r) ≈ 1 / family_a(Hierarchical(4, 5), 100, 50, 100)
    @test effective_size(r) < 15 < 90 < effective_size(r; from=:sexes)
    # Under FamilySize with Gamma weights it varies around the closed form.
    m = FamilySize(dams=0.5, sires=:random, nmale=10)
    r = fisher_wright(40, 2_000, [1_000], 1.0; result=true, mating=m, mutate=false, rng=Xoshiro(5))
    @test abs(effective_size(r) * family_a(m, 40, 10, 40) - 1) < 0.03

    # Binomial sex is redrawn until the breeders fit: 12 parents, 2 sires
    # and 8 dams.
    r = fisher_wright(12, 200, [100_000], 1.0; result=true, mating=Hierarchical(2, 4), rng=Xoshiro(3))
    @test all(m -> 2 <= m <= 4, r.males) && length(unique(r.males)) == 3
    # Infeasible designs are rejected before the run.
    @test_throws "generation 1" fisher_wright(20, 5, [1_000], 1.0; mating=Hierarchical(3, 4; nmale=2))
    @test_throws "generation 1" fisher_wright(20, 5, [1_000], 1.0; mating=Hierarchical(2, 5; nmale=12))
    @test_throws "generation 3" fisher_wright(20, 5, [1_000], 1.0; popsize=[20, 9, 20, 20, 20],
        mating=Hierarchical(3, 3))
    @test_throws ArgumentError fisher_wright(20, 5, [1_000], 1.0; mating=FamilySize(nsires=3, sexratio=0.1))

    # Through fisher_wright: valid output, census and males recorded,
    # nmale checked against every parental cohort, no selection.
    sizes = [fill(12, 5); fill(30, 5)]
    r = fisher_wright(20, 10, [100_000], 1.0; result=true, popsize=sizes,
        mating=FamilySize(dams=:equal, sires=:equal, nmale=5), rng=Xoshiro(2))
    @test r.census == [20; sizes] && r.males == fill(5, 10) && length(r.active_haplotypes) == 60
    @test_throws ArgumentError fisher_wright(20, 10, [100_000], 1.0; popsize=sizes, mating=FamilySize(nmale=12))
    for m in (Hierarchical(2, 3), FamilySize(sires=:equal), FamilySize(dams=:equal))
        @test_throws "nothing to select on" fisher_wright(20, 10, [100_000], 1.0; result=true, mating=m,
            fitness=SiteFitness(FitnessLoci([5], 0.1)))
    end
end

@testset "Family sizes under selection" begin
    # Parent usage: with fitness weights w and :random splits, each dam has
    # no·w/Σw offspring on average among the dams, and each sire no·w/Σw
    # among the sires. 12 parents, 4 of them males (a new random subset in
    # every replicate), 30 offspring, 20,000 replicates; z-score of the mean
    # deviation of every individual in each role.
    N, no, R = 12, 30, 20_000
    w = [0.2, 1.0, 3.0, 0.5, 2.0, 1.0, 0.1, 4.0, 1.5, 0.7, 2.5, 1.2]
    f = IndexFitness(w)
    rng = Xoshiro(21)
    st = SimState([UInt32[] for _ = 1:2N], UInt32[], 40, rng)
    st.w[1:N] .= w
    m = FamilySize(nmale=4)
    dev_d, dev_s = zeros(R, N), zeros(R, N)      # deviation from no·w/Σw
    uni_d = zeros(R, N)                          # deviation from uniform usage
    onesire = true
    for r = 1:R
        _mate!(m, f, st, no, rng)
        pm = st.pm[1:no, :]
        sires, dams = findall(st.sex[1:N]), findall(.!st.sex[1:N])
        for (dev, members, col) in ((dev_d, dams, 2), (dev_s, sires, 1))
            W = sum(w[members])
            for i in members
                dev[r, i] = count(==(i), pm[:, col]) - no * w[i] / W
            end
        end
        for i in dams
            uni_d[r, i] = count(==(i), pm[:, 2]) - no / length(dams)
        end
        onesire &= all(length(unique(pm[pm[:, 2] .== d, 1])) == 1 for d in unique(pm[:, 2]))
    end
    @test onesire
    zmax(dev) = maximum(abs(mean(c)) / (std(c) / sqrt(R)) for c in eachcol(dev))
    @test zmax(dev_d) < 4 && zmax(dev_s) < 4      # 12 individuals per role
    # Usage does follow fitness: against uniform usage the deviations are
    # far from zero.
    @test zmax(uni_d) > 20
    # Under :equal splits selection is refused here as well.
    @test_throws "nothing to select on" _mate!(FamilySize(dams=:equal), f, st, no, rng)

    # One generation of selection on a trait is unbiased under :random
    # splits, as under random mating: the offspring mean G given the parents
    # is the fitness-weighted mean within each sex, averaged over sexes.
    chr = [5_000, 5_000]
    rmap = uniform_recombination_map(chr; M=2_000)
    N = 30
    lt = TraitLoci(collect(100:100:10_000), randn(rng, 100) .* 0.3)
    start = [lt.pos[rand(rng, 100) .< 0.3] for _ = 1:2N]
    tf = TraitFitness(lt; σe=0.7, fitness=Directional(0.8))
    R = 3_000
    d, shift = zeros(R), zeros(R)
    for r = 1:R
        st = SimState([UInt32[] for _ = 1:2N], UInt32[], N, rng, copy.(start), UInt32[])
        _simulate!(st, ConstantSize(N), FamilySize(), tf, 1, rmap, Poisson(1.0), UInt32(1):UInt32(10_000);
            mutate=false, track=false, fixation_interval=1, verbose=false, progress_interval=100, rng=rng)
        G, wt, sex = st.gval[1:N], st.w[1:N], st.sex[1:N]
        expect = (sum(G[sex] .* wt[sex]) / sum(wt[sex]) + sum(G[.!sex] .* wt[.!sex]) / sum(wt[.!sex])) / 2
        got = mean(_genotypic_value(lt, st.sel[2i-1], st.sel[2i]) for i = 1:N)
        d[r] = got - expect
        shift[r] = got - mean(G)
    end
    @test abs(mean(d)) / (std(d) / sqrt(R)) < 3.5
    @test mean(shift) / (std(shift) / sqrt(R)) > 10      # selection acts

    # Through fisher_wright: selection with family sizes runs, responds, and
    # lowers the realized effective size below that of the neutral scheme.
    kw = (M=2_000, mut_base=1e4, result=true)
    loci = TraitLoci(collect(100:100:10_000), [0.2 * (-1)^k for k = 1:100]; μ=2e-3)
    m = FamilySize(dams=1.0, sires=:random, nmale=10)
    up = fisher_wright(N, 150, chr, 2.0; kw..., mating=m, fitness=TraitFitness(loci; σe=0.5, fitness=Directional(1.0)),
        rng=Xoshiro(22))
    flat = fisher_wright(N, 150, chr, 2.0; kw..., mating=m, rng=Xoshiro(22))
    @test up.mean_genotypic[end] > 2
    @test effective_size(up, 50:150) < effective_size(flat, 50:150)
    @test_throws "nothing to select on" fisher_wright(N, 5, chr, 2.0; kw..., mating=Hierarchical(3, 3),
        fitness=TraitFitness(loci; fitness=Directional(1.0)))
end

@testset "Monoecious mating" begin
    @test_throws ArgumentError Monoecious(selfing=-0.1)
    @test_throws ArgumentError Monoecious(selfing=1.5)
    @test sprint(show, Monoecious()) == "Monoecious()"
    @test sprint(show, Monoecious(selfing=0.25)) == "Monoecious(selfing = 0.25)"

    # The mating table: parents from the whole cohort, no sexes, and the
    # selfed offspring at the requested rate (binomial z-score).
    rng = Xoshiro(31)
    N = 40
    st = SimState([UInt32[] for _ = 1:2N], UInt32[], 60, rng)
    for (s, no) in ((nothing, 40), (0.0, 57), (1.0, 13), (0.3, 60))
        selfed = 0
        inrange = true
        R = 2_000
        for _ = 1:R
            _mate!(Monoecious(selfing=s), Neutral(), st, no, rng)
            pm = st.pm[1:no, :]
            inrange &= all(1 .<= pm .<= N)
            selfed += count(pm[:, 1] .== pm[:, 2])
        end
        @test inrange
        @test isempty(st.sires) && isempty(st.dams)
        p = s === nothing ? 1 / N : s
        if 0 < p < 1
            @test abs(selfed / (R * no) - p) / sqrt(p * (1 - p) / (R * no)) < 3.5
        else
            @test selfed == p * R * no
        end
    end

    # Pairs of offspring sharing a parent, against the brute-force count
    # over ordered pairs (i, j), i ≠ j, of the roles of i's parent in j:
    # the cross term vanishes with separate sexes.
    for m in (Monoecious(), Monoecious(selfing=0.5), RandomMating()), no in (30, 45)
        _mate!(m, Neutral(), st, no, rng)
        empty!(st.ppairs), empty!(st.mpairs), empty!(st.selfing)
        _record_families!(st, no)
        pm = st.pm[1:no, :]
        pp = sum((pm[i, 1] == pm[j, 1]) + (pm[i, 1] == pm[j, 2]) for i = 1:no, j = 1:no if i != j)
        mp = sum((pm[i, 2] == pm[j, 1]) + (pm[i, 2] == pm[j, 2]) for i = 1:no, j = 1:no if i != j)
        @test (st.ppairs[1], st.mpairs[1]) == (pp, mp)
        @test st.selfing[1] == count(pm[:, 1] .== pm[:, 2]) / no
        if m isa RandomMating
            @test pp == sum(c * (c - 1) for c in (count(==(k), pm[:, 1]) for k = 1:N))
        end
    end

    # Under selection: each parent is first parent no·w/Σw times on average,
    # and a parent of weight zero is never used. 12 parents, 30 offspring,
    # 20,000 replicates.
    Nw, no, R = 12, 30, 20_000
    w = [0.2, 1.0, 3.0, 0.5, 2.0, 1.0, 0.0, 4.0, 1.5, 0.7, 2.5, 1.2]
    f = IndexFitness(w)
    st = SimState([UInt32[] for _ = 1:2Nw], UInt32[], 40, rng)
    st.w[1:Nw] .= w
    for s in (nothing, 0.0, 0.4)
        firsts = zeros(R, Nw)
        used0, outcrossed = false, true
        for r = 1:R
            _mate!(Monoecious(selfing=s), f, st, no, rng)
            pm = st.pm[1:no, :]
            for k = 1:Nw
                firsts[r, k] = count(==(k), pm[:, 1])
            end
            used0 |= any(==(7), pm)
            s == 0 && (outcrossed &= all(pm[:, 1] .!= pm[:, 2]))
        end
        @test !used0 && outcrossed
        dev = firsts .- no .* (w ./ sum(w))'
        z = [mean(dev[:, k]) / (std(dev[:, k]) / sqrt(R)) for k = 1:Nw if w[k] > 0]
        @test maximum(abs, z) < 4
    end
    lone = zeros(Nw)
    lone[3] = 1.0
    st.w[1:Nw] .= lone
    @test_throws ErrorException _mate!(Monoecious(selfing=0.5), IndexFitness(lone), st, no, rng)
    _mate!(Monoecious(selfing=1), IndexFitness(lone), st, no, rng)
    @test all(==(3), st.pm[1:no, :])

    # A run: no males, the selfing record, hermaphrodites in the pedigree,
    # and a continuation that appends to the record.
    chr = [100_000]
    res = fisher_wright(30, 20, chr, 1.0; mating=Monoecious(selfing=0.5), result=true, pedigree=3, rng=rng)
    @test res.males == zeros(Int, 20)
    @test length(res.selfing) == 20 && all(0 .<= res.selfing .<= 1)
    @test all(ismissing, res.pedigree.sex)
    kids = res.pedigree[res.pedigree.generation .== 20, :]
    @test count(kids.sire .== kids.dam) / 30 == res.selfing[20]
    @test_throws ArgumentError effective_size(res; from=:sexes)
    more = fisher_wright(30, 5, chr, 1.0; init=res, mating=RandomMating(), result=true, rng=rng)
    @test length(more.selfing) == 25 && all(iszero, more.selfing[21:25])
    @test effective_size(more, 21:25; from=:sexes) > 0

    # E[a_g] = 1/N for every selfing model (the four parent slots of two
    # offspring are uniform and independent); 1,000 generations, z-score.
    for s in (nothing, 0.0, 0.5, 1.0)
        r = fisher_wright(30, 1_000, [1_000], 1.0; mut_base=1e12, mating=Monoecious(selfing=s), result=true, rng=rng)
        a = (r.paternal_pairs .+ r.maternal_pairs) ./ (4 * 30 * 29)
        @test abs(mean(a) - 1 / 30) / (std(a) / sqrt(length(a))) < 3.5
        σ = s === nothing ? 1 / 30 : s
        # Realized Ne about N(2 - σ)/2.
        @test isapprox(effective_size(r), 30 * (2 - σ) / 2; rtol=0.05)
    end

    # Within the reserved capacity, a generation allocates task overhead only.
    rmap = uniform_recombination_map([1_000, 1_000]; M=1e3)
    st = SimState([UInt32.(findall(rand(rng, Bool, 2_000))) for _ = 1:80], UInt32[], 40, rng)
    run!(st, m, nt) = _simulate!(st, ConstantSize(40), m, Neutral(), nt, rmap, Poisson(1.0), UInt32(1):UInt32(2_000);
        mutate=false, track=false, fixation_interval=1, verbose=false, progress_interval=100, rng=rng)
    for m in (Monoecious(), Monoecious(selfing=0.3))
        run!(st, m, 50)
        @test @allocated(run!(st, m, 50)) / 50 < 2^14
    end
end

# Whether the parents `x` and `y` (IDs) are related at `level` (1 full sibs,
# 2 sibs, 3 sibs or first cousins), from a pedigree DataFrame; 0 is an
# unknown parent.
function related_by_pedigree(par, x, y, level)
    px = filter(!=(0), collect(get(par, x, (0, 0))))
    py = filter(!=(0), collect(get(par, y, (0, 0))))
    level == 1 && return length(px) == 2 && sort(px) == sort(py)
    isempty(intersect(px, py)) || return true
    level == 2 && return false
    gp(ps) = filter(!=(0), reduce(vcat, [collect(get(par, p, (0, 0))) for p in ps]; init=Int[]))
    return !isempty(intersect(gp(px), gp(py)))
end

@testset "Avoidance of relatives" begin
    @test_throws ArgumentError AvoidRelatives(FamilySize())
    @test_throws ArgumentError AvoidRelatives(; level=:sibs)
    @test_throws ArgumentError AvoidRelatives(; max_tries=0)
    @test sprint(show, AvoidRelatives(RandomMating(nmale=3); level=:cousin)) ==
          "AvoidRelatives(RandomMating(nmale = 3); level = :cousin, max_tries = 100)"

    # The relationship test on a hand-made record: columns are individuals,
    # rows the two parents and four grandparents (0 unknown).
    rel = zeros(Int, 6, 6)
    rel[:, 1] = [1, 2, 5, 6, 7, 8]
    rel[:, 2] = [2, 1, 7, 8, 5, 6]     # full sib of 1, roles swapped
    rel[:, 3] = [1, 3, 5, 6, 9, 9]     # half sib of 1
    rel[:, 4] = [4, 9, 6, 0, 0, 0]     # first cousin of 1 through grandparent 6
    rel[:, 5] = [0, 0, 0, 0, 0, 0]     # unknown
    rel[:, 6] = [4, 9, 6, 0, 0, 0]     # full sib of 4
    @test [_related(rel, 1, y, l) for y = 2:5, l = 1:3] ==
          Bool[1 1 1; 0 1 1; 0 0 1; 0 0 0]
    @test _related(rel, 4, 6, 1) && !_related(rel, 5, 5, 1)

    # The second parent, given the first: uniform among the unrelated others,
    # or in proportion to fitness among them. Monoecious without selfing,
    # 8 individuals, individual 1 a half sib of 2 and 3; 40,000 offspring.
    rng = Xoshiro(41)
    N, no = 8, 40_000
    st = SimState([UInt32[] for _ = 1:2N], UInt32[], no, rng)
    st.rel = zeros(Int, 6, no)
    st.rnew = zeros(Int, 6, no)
    st.rel[1:2, 1:3] = [1 1 1; 2 3 4]
    push!(st.forced, 0)
    w = [1.0, 1.0, 1.0, 0.5, 1.0, 2.0, 3.0, 0.0]
    for (f, p) in ((Neutral(), [0, 0, 0, 1, 1, 1, 1, 1] / 5), (IndexFitness(w), [0, 0, 0, 0.5, 1, 2, 3, 0] / 6.5))
        f isa Neutral || (st.w[1:N] .= w)
        _mate!(AvoidRelatives(Monoecious(selfing=0)), f, st, no, rng)
        pm = st.pm[1:no, :]
        q = pm[pm[:, 1] .== 1, 2]
        c = [count(==(k), q) for k = 1:N]
        @test all(c[p .== 0] .== 0)
        e = length(q) .* p
        @test sum((c[p .> 0] .- e[p .> 0]) .^ 2 ./ e[p .> 0]) < quantile(Chisq(count(p .> 0) - 1), 0.999)
        @test st.forced[end] == 0
    end

    # With max_tries = 1 nothing is redrawn: the realization of the inner
    # scheme, with its matings of relatives counted.
    chr = [100_000]
    for inner in (RandomMating(nmale=5), Monoecious(selfing=0.2))
        a = fisher_wright(30, 30, chr, 1.0; mating=inner, result=true, rng=Xoshiro(5))
        b = fisher_wright(30, 30, chr, 1.0; mating=AvoidRelatives(inner; max_tries=1), result=true, rng=Xoshiro(5))
        @test a.active_haplotypes == b.active_haplotypes
        @test sum(b.related_matings) > 0 && sum(a.related_matings) == 0
    end

    # Against the pedigree: in every generation, the matings of relatives
    # at the level avoided are exactly those counted, and none with ample
    # retries where unrelated mates abound.
    for (m, some) in ((AvoidRelatives(RandomMating(nmale=5); level=:cousin), true),
        (AvoidRelatives(RandomMating(); level=:halfsib), false),
        (AvoidRelatives(Monoecious(); level=:cousin), false),
        (AvoidRelatives(Monoecious(selfing=0.0); level=:fullsib, max_tries=1), true))
        res = fisher_wright(30, 25, chr, 1.0; mating=m, result=true, pedigree=typemax(Int), rng=rng)
        ped = res.pedigree
        par = Dict(r.id => (r.sire, r.dam) for r in eachrow(ped))
        counts = [count(r -> r.sire != r.dam && related_by_pedigree(par, r.sire, r.dam, m.level),
            eachrow(ped[ped.generation .== g, :])) for g = 1:25]
        @test counts == res.related_matings
        @test (sum(counts) > 0) == some
    end

    # Selfing is decided as without avoidance: random selfing stays at 1/N.
    res = fisher_wright(20, 400, chr, 1.0; mating=AvoidRelatives(Monoecious(); level=:cousin), result=true, rng=rng)
    @test abs(mean(res.selfing) - 1 / 20) / sqrt(1 / 20 * 19 / 20 / (20 * 400)) < 3.5

    # A continuation starts the pedigree anew and appends to the record.
    more = fisher_wright(20, 5, chr, 1.0; init=res, mating=AvoidRelatives(), result=true, rng=rng)
    @test length(more.related_matings) == 405

    # Within the reserved capacity, a generation allocates task overhead only.
    rmap = uniform_recombination_map([1_000, 1_000]; M=1e3)
    st = SimState([UInt32.(findall(rand(rng, Bool, 2_000))) for _ = 1:80], UInt32[], 40, rng)
    run!(st, m, nt) = _simulate!(st, ConstantSize(40), m, Neutral(), nt, rmap, Poisson(1.0), UInt32(1):UInt32(2_000);
        mutate=false, track=false, fixation_interval=1, verbose=false, progress_interval=100, rng=rng)
    for m in (AvoidRelatives(; level=:cousin), AvoidRelatives(Monoecious(); level=:halfsib))
        run!(st, m, 50)
        @test @allocated(run!(st, m, 50)) / 50 < 2^14
    end
end

@testset "Custom mating" begin
    chr = [100_000]
    # Callbacks that draw as the built-in schemes do give their realizations.
    mono = CustomMating() do ctx
        for i in 1:ctx.noffspring
            ctx.pm[i, 1] = rand(ctx.rng, 1:ctx.n)
            ctx.pm[i, 2] = rand(ctx.rng, 1:ctx.n)
        end
    end
    sires, dams = Int[], Int[]
    twosex = CustomMating(; sexes=true) do ctx
        while true
            rand!(ctx.rng, ctx.sex)
            0 < count(ctx.sex) < ctx.n && break
        end
        empty!(sires), empty!(dams)
        for k in 1:ctx.n
            push!(ctx.sex[k] ? sires : dams, k)
        end
        for i in 1:ctx.noffspring
            ctx.pm[i, 1] = rand(ctx.rng, sires)
            ctx.pm[i, 2] = rand(ctx.rng, dams)
        end
    end
    for (builtin, custom) in ((Monoecious(), mono), (RandomMating(), twosex))
        a = fisher_wright(30, 30, chr, 1.0; mating=builtin, result=true, pedigree=2, rng=Xoshiro(3))
        b = fisher_wright(30, 30, chr, 1.0; mating=custom, result=true, pedigree=2, rng=Xoshiro(3))
        @test a.active_haplotypes == b.active_haplotypes
        @test (a.males, a.selfing, a.paternal_pairs) == (b.males, b.selfing, b.paternal_pairs)
        @test isequal(a.pedigree, b.pedigree)
    end

    # The context: generation numbers continued through init, IDs, parents
    # and grandparents as in the pedigree, haplotypes, fitness and phenotype.
    seen = []
    record = CustomMating() do ctx
        push!(seen, (ctx.generation, ctx.n, ctx.noffspring, collect(ctx.ids), copy(ctx.relatives),
            copy(ctx.haplotypes), copy(ctx.fitness), copy(ctx.phenotypic), copy(ctx.genotypic)))
        for i in 1:ctx.noffspring
            ctx.pm[i, 1] = rand(ctx.rng, 1:ctx.n)
            ctx.pm[i, 2] = rand(ctx.rng, 1:ctx.n)
        end
    end
    rng = Xoshiro(51)
    first = fisher_wright(20, 3, chr, 1.0; mating=record, popsize=[20, 25, 22], result=true,
        pedigree=typemax(Int), rng=rng)
    @test [s[1] for s in seen] == 1:3 && [s[2] for s in seen] == [20, 20, 25] && [s[3] for s in seen] == [20, 25, 22]
    ped = first.pedigree
    for (g, s) in enumerate(seen)
        @test s[4] == ped.id[ped.generation .== g-1]
        @test all(s[7] .== 1) && all(iszero, s[8])
    end
    # Parents (rows 1-2) of the cohort of generation 2, as indices into
    # generation 1, and their grandparents (rows 3-6) into generation 0.
    idx(id, g) = id - minimum(ped.id[ped.generation .== g]) + 1
    gen2 = ped[ped.generation .== 2, :]
    @test seen[3][5][1, :] == idx.(gen2.sire, 1) && seen[3][5][2, :] == idx.(gen2.dam, 1)
    par1 = Dict(r.id => (r.sire, r.dam) for r in eachrow(ped[ped.generation .== 1, :]))
    @test seen[3][5][3, :] == [idx(par1[s][1], 0) for s in gen2.sire]
    @test seen[3][5][6, :] == [idx(par1[d][2], 0) for d in gen2.dam]
    @test all(iszero, seen[1][5]) && all(iszero, seen[2][5][3:6, :])
    empty!(seen)
    more = fisher_wright(22, 2, chr, 1.0; init=first, mating=record, result=true, rng=rng)
    @test [s[1] for s in seen] == [4, 5] && seen[1][4] == first.first_id .+ (0:21)
    @test seen[1][6] == first.active_haplotypes

    empty!(seen)
    qtl = TraitLoci([10_000, 60_000], [1.0, -0.5])
    tf = TraitFitness(qtl; σe=1.0, fitness=Stabilizing(0.0, 4.0))
    sel = fisher_wright(20, 3, chr, 1.0; mating=record, fitness=tf, result=true, rng=rng)
    @test [mean(s[8]) for s in seen] ≈ sel.mean_phenotypic
    # Relative fitness: scaled by its largest value in the cohort.
    @test all(s -> s[7] ≈ exp.(-s[8] .^ 2 ./ 8) ./ maximum(exp.(-s[8] .^ 2 ./ 8)), seen)

    # Checks of the table.
    bad = CustomMating(ctx -> (ctx.pm[1:ctx.noffspring, :] .= 1; ctx.pm[2, 2] = ctx.n + 1))
    @test_throws ArgumentError fisher_wright(10, 2, chr, 1.0; mating=bad)
    wrongsex = CustomMating(ctx -> (ctx.sex .= true; ctx.pm[1:ctx.noffspring, :] .= 1); sexes=true)
    @test_throws ArgumentError fisher_wright(10, 2, chr, 1.0; mating=wrongsex)
    # Without sexes, selfing is allowed and recorded.
    selfer = CustomMating(ctx -> (for i in 1:ctx.noffspring; ctx.pm[i, :] .= rand(ctx.rng, 1:ctx.n); end))
    r = fisher_wright(10, 4, chr, 1.0; mating=selfer, result=true, rng=rng)
    @test r.selfing == ones(4) && r.males == zeros(Int, 4)
    @test effective_size(r) ≈ 1 / mean((r.paternal_pairs .+ r.maternal_pairs) ./ (4 * 10 * 9) .* 2)

    # A generation of an allocation-free callback allocates task overhead only.
    rmap = uniform_recombination_map([1_000, 1_000]; M=1e3)
    st = SimState([UInt32.(findall(rand(rng, Bool, 2_000))) for _ = 1:80], UInt32[], 40, rng)
    run!(st, m, nt) = _simulate!(st, ConstantSize(40), m, Neutral(), nt, rmap, Poisson(1.0), UInt32(1):UInt32(2_000);
        mutate=false, track=false, fixation_interval=1, verbose=false, progress_interval=100, rng=rng)
    run!(st, mono, 50)
    @test @allocated(run!(st, mono, 50)) / 50 < 2^14
end

@testset "Pedigree" begin
    chr, kw = [100_000, 100_000], (M=5e4, mut_base=1e5, result=true)
    sizes = [12, 8, 15, 15, 10]
    r = fisher_wright(10, 5, chr, 5.0; kw..., popsize=sizes, pedigree=typemax(Int), rng=Xoshiro(1))
    ped = r.pedigree
    @test names(ped) == ["id", "sire", "dam", "sex", "generation"]
    # IDs in birth order, one row per individual of every generation.
    @test ped.id == 1:sum(r.census) && ped.generation == reduce(vcat, [fill(g, n) for (g, n) in zip(0:5, r.census)])
    @test r.first_id == sum(r.census) - r.census[end] + 1
    # Founders have no parents, the last generation no sex yet; parents are
    # of the previous generation, sires male and dams female.
    @test all(==(0), ped.sire[ped.generation .== 0]) && all(==(0), ped.dam[ped.generation .== 0])
    @test all(ismissing, ped.sex[ped.generation .== 5]) && !any(ismissing, ped.sex[ped.generation .< 5])
    row = Dict(id => i for (i, id) in enumerate(ped.id))
    off = ped[ped.generation .> 0, :]
    @test all(ped.generation[row[s]] == g - 1 && ped.sex[row[s]] === true for (s, g) in zip(off.sire, off.generation))
    @test all(ped.generation[row[d]] == g - 1 && ped.sex[row[d]] === false for (d, g) in zip(off.dam, off.generation))
    # The other records agree: males among the parents, and pairs of
    # offspring sharing a sire or a dam.
    @test [count(x -> x === true, ped.sex[ped.generation .== g-1]) for g = 1:5] == r.males
    pairs(v) = (c = values(Dict(k => count(==(k), v) for k in unique(v))); sum(n -> n * (n - 1), c))
    @test [pairs(off.sire[off.generation .== g]) for g = 1:5] == r.paternal_pairs
    @test [pairs(off.dam[off.generation .== g]) for g = 1:5] == r.maternal_pairs

    # The pedigree describes the transmission: every haplotype of the start
    # carries a private marker, there is no mutation, and every individual
    # carries markers of its founder ancestors only.
    N0 = 12
    marker(k) = UInt32(1_000 * k)                    # haplotype k of the start
    start = [[marker(k)] for k = 1:2N0]
    r = fisher_wright(N0, 4, chr, 5.0; kw..., init=start, mutate=false, pedigree=typemax(Int), rng=Xoshiro(2))
    ped = r.pedigree
    row = Dict(id => i for (i, id) in enumerate(ped.id))
    founders(id) = ped.sire[row[id]] == 0 ? [id] : vcat(founders(ped.sire[row[id]]), founders(ped.dam[row[id]]))
    for i = 1:length(r.active_haplotypes)÷2
        id = r.first_id + i - 1
        allowed = Set(marker(k) for f in founders(id) for k in (2f - 1, 2f))
        @test issubset(union(r.active_haplotypes[2i-1], r.active_haplotypes[2i]), allowed)
    end
    # ... and exactly in one generation: the paternal haplotype comes from
    # the sire, the maternal one from the dam.
    r = fisher_wright(N0, 1, chr, 5.0; kw..., init=start, mutate=false, pedigree=1, rng=Xoshiro(3))
    ped = r.pedigree[r.pedigree.generation .== 1, :]
    @test ped.id == r.first_id .+ (0:N0-1)
    for i = 1:N0
        s, d = ped.sire[i], ped.dam[i]
        @test issubset(r.active_haplotypes[2i-1], [marker(2s - 1), marker(2s)])
        @test issubset(r.active_haplotypes[2i], [marker(2d - 1), marker(2d)])
    end

    # A window of the last k generations: k + 1 generations, the oldest as
    # founders of the record.
    w = fisher_wright(10, 6, chr, 5.0; kw..., pedigree=2, rng=Xoshiro(4))
    full = fisher_wright(10, 6, chr, 5.0; kw..., pedigree=typemax(Int), rng=Xoshiro(4))
    @test w.active_haplotypes == full.active_haplotypes      # recording does not change the run
    @test sort(unique(w.pedigree.generation)) == 4:6
    @test all(==(0), w.pedigree.sire[w.pedigree.generation .== 4])
    tail = full.pedigree[full.pedigree.generation .>= 5, :]
    @test isequal(w.pedigree[w.pedigree.generation .>= 5, :], tail)
    @test w.pedigree.sex[w.pedigree.generation .== 4] == full.pedigree.sex[full.pedigree.generation .== 4]

    # Without the keyword: nothing recorded, IDs counted, and the run is the
    # same.
    none = fisher_wright(10, 6, chr, 5.0; kw..., rng=Xoshiro(4))
    @test nrow(none.pedigree) == 0 && none.first_id == full.first_id
    @test none.active_haplotypes == full.active_haplotypes

    # A continued run continues generations and IDs and appends to the
    # pedigree of init, within the window; init is left alone.
    a = fisher_wright(10, 3, chr, 5.0; kw..., pedigree=typemax(Int), rng=Xoshiro(5))
    keep = copy(a.pedigree)
    b = fisher_wright(10, 3, chr, 5.0; kw..., init=a, pedigree=typemax(Int), rng=Xoshiro(6))
    @test b.pedigree.id == 1:70 && b.first_id == 61 && maximum(b.pedigree.generation) == 6
    @test isequal(b.pedigree[1:30, :], a.pedigree[1:30, :])     # generations 0 to 2
    @test !any(ismissing, b.pedigree.sex[b.pedigree.generation .== 3])   # sexes of the start filled in
    @test isequal(a.pedigree, keep)
    c = fisher_wright(10, 3, chr, 5.0; kw..., init=a, pedigree=4, rng=Xoshiro(6))
    @test sort(unique(c.pedigree.generation)) == 2:6 && all(==(0), c.pedigree.sire[c.pedigree.generation .== 2])
    @test isequal(c.pedigree[c.pedigree.generation .>= 3, :], b.pedigree[b.pedigree.generation .>= 3, :])
    # From plain haplotypes the numbering starts at 1.
    plain = fisher_wright(10, 2, chr, 5.0; kw..., init=a.active_haplotypes, pedigree=1, rng=Xoshiro(7))
    @test plain.first_id == 21 && minimum(plain.pedigree.generation) == 1

    # Family sizes: every dam has one sire.
    f = fisher_wright(20, 4, chr, 5.0; kw..., mating=Hierarchical(2, 3), pedigree=typemax(Int), rng=Xoshiro(8))
    off = f.pedigree[f.pedigree.generation .> 0, :]
    @test all(length(unique(off.sire[off.dam .== d])) == 1 for d in unique(off.dam))
    @test all(g -> length(unique(off.sire[off.generation .== g])) == 2, 1:4)

    # Errors, and the Species method.
    @test_throws ArgumentError fisher_wright(10, 2, chr, 5.0; pedigree=1)
    @test_throws ArgumentError fisher_wright(10, 2, chr, 5.0; kw..., pedigree=-1)
    sp = GenericSpecies("Ped", Int32(10), UInt32.(chr), UInt32(50_000))
    @test isequal(fisher_wright(sp, 6, 5.0; mut_base=1e5, result=true, pedigree=2, rng=Xoshiro(4)).pedigree, w.pedigree)
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

@testset "Heterozygosity decay under family-size mating" begin
    # Drift alone from a polymorphic start, under five FamilySize schemes
    # (two with breeder subsets),
    # against the exact recursion with a = (P_f + P_m)/4 from the closed
    # forms of the splits (`family_a`): 20 parents, 15 generations, 400
    # replicates of 1,000 sites on 20 chromosomes of 1 Morgan. Each scheme
    # is also separated from RandomMating with the same sexes.
    N, chr, reps, t = 20, fill(50, 20), 400, 15
    L = sum(chr)
    rng = Xoshiro(20261116)
    cases = (
        (FamilySize(dams=:equal, sires=:equal, nmale=5), RandomMating(nmale=5)),
        (FamilySize(dams=0.5, sires=:random, nmale=10), RandomMating(nmale=10)),
        (FamilySize(dams=:random, sires=:random), RandomMating()),
        (Hierarchical(3, 2), RandomMating()),
        (FamilySize(dams=1.0, sires=:random, nsires=4, ndams=8, nmale=8), RandomMating(nmale=8)),
    )
    binom = [Float64(binomial(big(N), k)) for k = 1:N-1]
    function expected(a)
        hB = hW = 0.5
        for _ = 1:t
            hB, hW = (1 - a) * hB + a / 2 * hW, hB
        end
        return hW / 2N + (1 - 1 / N) * hB
    end
    function coalescence(m, N)
        nm = m isa FamilySize ? m.sex.nmale : m.nmale
        f(k) = m isa FamilySize ? family_a(m, N, k, N) : 1 / 4k + 1 / 4(N - k)
        nm > 0 && return f(nm)
        # Binomial sex conditioned on enough of each sex for the breeders.
        lo, hi = m isa FamilySize ? (max(m.nsires, 1), N - max(m.ndams, 1)) : (1, N - 1)
        return sum(binom[k] * f(k) for k = lo:hi) / sum(binom[lo:hi])
    end
    for (m, alt) in cases
        H = zeros(reps)
        for r = 1:reps
            start = [UInt32.(findall(rand(rng, Bool, L))) for _ = 1:2N]
            haps, _ = fisher_wright(N, t, chr, 1.0; M=50, init=start, mutate=false, mating=m, rng=rng)
            H[r] = site_heterozygosity(haps, L)
        end
        se = sqrt(var(H) / reps)
        @test abs(mean(H) - expected(coalescence(m, N))) / se < 3.5
        @test abs(mean(H) - expected(coalescence(alt, N))) / se > 3.5
    end
end

@testset "Heterozygosity decay under monoecious mating" begin
    # Drift alone from a polymorphic start, hermaphrodites with selfing rate
    # σ (random selfing: σ = 1/N), against the exact recursion of the
    # heterozygosity within (h_W) and between (h_B) individuals:
    #   h_B' = (1 - 1/N) h_B + h_W/(2N),  h_W' = (1 - σ) h_B + σ h_W/2,
    # since the two genes of different offspring come from independent,
    # uniform parents, and the two genes of one offspring from one parent
    # with probability σ. Both Σ2pq and the observed heterozygosity are
    # tested; 20 parents, 15 generations, 400 replicates of 1,000 sites on 20
    # chromosomes of 1 Morgan. Each case is also separated from no selfing.
    N, chr, reps, t = 20, fill(50, 20), 400, 15
    L = sum(chr)
    rng = Xoshiro(20261117)
    function expected(σ)
        hB = hW = 0.5
        for _ = 1:t
            hB, hW = (1 - 1 / N) * hB + hW / 2N, (1 - σ) * hB + σ * hW / 2
        end
        return hW / 2N + (1 - 1 / N) * hB, hW
    end
    observed(haps) = mean(length(symdiff(haps[2i-1], haps[2i])) for i = 1:N) / L
    for s in (nothing, 0.0, 0.5, 1.0)
        H, O = zeros(reps), zeros(reps)
        for r = 1:reps
            start = [UInt32.(findall(rand(rng, Bool, L))) for _ = 1:2N]
            haps, _ = fisher_wright(N, t, chr, 1.0; M=50, init=start, mutate=false,
                mating=Monoecious(selfing=s), rng=rng)
            H[r] = site_heterozygosity(haps, L)
            O[r] = observed(haps)
        end
        σ = s === nothing ? 1 / N : s
        eH, eO = expected(σ)
        @test abs(mean(H) - eH) / sqrt(var(H) / reps) < 3.5
        @test abs(mean(O) - eO) / sqrt(var(O) / reps) < 3.5
        σ >= 0.5 && @test abs(mean(O) - expected(0.0)[2]) / sqrt(var(O) / reps) > 3.5
    end
end

# Fixation of one new copy of a selected allele in the engine: `reps`
# replicates of N individuals, each run to loss or fixation.
function engine_fixation(N, s, h, reps, rng)
    rmap = uniform_recombination_map([100]; M=1e8)
    f = SiteFitness(FitnessLoci([50], s; h=h))
    st = SimState([UInt32[] for _ = 1:2N], UInt32[], N, rng)
    fixed = 0
    for _ = 1:reps
        foreach(empty!, st.sel)
        push!(st.sel[rand(rng, 1:2N)], UInt32(50))
        while true
            _simulate!(st, ConstantSize(N), RandomMating(), f, 1, rmap, Poisson(1.0), UInt32(1):UInt32(100);
                mutate=false, track=false, fixation_interval=1, verbose=false, progress_interval=100, rng=rng)
            c = count(!isempty, view(st.sel, 1:2N))
            c == 0 && break
            c == 2N && (fixed += 1; break)
        end
        resize!(st.census, 1)
        empty!(st.males)
    end
    return fixed / reps
end

@testset "Fixation probability under selection" begin
    # One copy of a selected allele in N = 30, followed to loss or fixation,
    # 4,000 replicates per case. The reference is an independent, naive
    # simulation of the same two-sex model with 2,000,000 replicates
    # (docs/src/validation/scripts/fixation-reference.jl); selection within
    # each sex and the binomial sex ratio put it 6 to 13% below the ideal
    # Wright-Fisher value at this N. z-scores of the difference of two
    # proportions. The three cases with s = 0.2 differ only in dominance and
    # are at least 12 standard errors apart, so an error in the fitness of
    # the heterozygote is caught; s = 0 checks that equal fitness is drift.
    reps = 4_000
    rng = Xoshiro(20261114)
    for (s, h, ref, se_ref) in FIXATION_REFERENCE
        u = engine_fixation(30, s, h, reps, rng)
        @test abs(u - ref) / sqrt(ref * (1 - ref) / reps + se_ref^2) < 3.5
    end
end

@testset "Selection under recurrent mutation" begin
    # 10 deleterious sites (s = -0.1, h = 0.5, μ = 0.002 per site), one per
    # chromosome so that they recombine freely, N = 30, 300 generations from
    # a population without variation, 1,000 runs. At the end, the mean
    # derived-allele frequency per site and the proportion of sites without
    # the derived allele are compared with an independent naive simulation
    # of one site of the same model, 400,000 replicates
    # (docs/src/validation/validate-selection-balance.txt, case 1). The
    # variance of log fitness among parents is about 0.002, so selection at
    # the other sites changes the effective size by less than 1%. A wrong
    # dominance (h = 0.1: q = 0.0567) is 5.3 standard errors away.
    N, K, T, runs = 30, 10, 300, 1_000
    loci = FitnessLoci(collect(50:100:100K), -0.1; h=0.5, μ=2e-3)
    f = SiteFitness(loci)
    rng = Xoshiro(20261115)
    q, zero = zeros(runs), zeros(runs)
    for r = 1:runs
        res = fisher_wright(N, T, fill(100, K), 0.02; mut_base=1e12, fitness=f, result=true, rng=rng)
        c = zeros(Int, K)
        for hap in res.selected_haplotypes, x in hap
            c[searchsortedfirst(loci.pos, x)] += 1
        end
        for x in res.selected_substitutions
            c[searchsortedfirst(loci.pos, x)] = 2N
        end
        q[r] = mean(c) / 2N
        zero[r] = count(==(0), c) / K
    end
    z(x, ref, se_ref) = (mean(x) - ref) / sqrt(var(x) / runs + se_ref^2)
    @test abs(z(q, 0.06502, 0.00024)) < 3.5
    @test abs(z(zero, 0.50406, 0.00079)) < 3.5
    # The test can tell h = 0.5 from h = 0.1 (reference q = 0.05668).
    @test z(q, 0.05668, 0.00022) > 3.5
end

@testset "Neutral genetic variance of a trait" begin
    # The genetic variance of a neutral additive trait under mutation and
    # drift against its exact expectation E[α²](2H_B - H_W), α = a/2, with
    # H_B and H_W from the two-sex recursion with 2U new alleles per pair of
    # genes (docs/src/validation/validate-trait.txt). N = 50, 20,000 loci,
    # a ~ Normal(0, 0.2²) drawn anew for each of 100 replicates, μ = 1e-5,
    # 200 generations from no variation. A heterozygote coded as a instead
    # of a/2 would about quadruple the variance (most derived alleles are
    # rare, so mostly heterozygous), and a mutation rate off by a factor 2
    # double it.
    N, T, K, μ, σa, reps = 50, 200, 20_000, 1e-5, 0.2, 100
    chr = fill(1_000_000, 10)
    rng = Xoshiro(20261117)
    V = zeros(reps, T)
    for r = 1:reps
        loci = trait_loci(chr, K; a=Normal(0, σa), μ=μ, rng=rng)
        res = fisher_wright(N, T, chr, 0.02; mut_base=1e12, result=true, rng=rng,
            fitness=TraitFitness(loci; fitness=P -> 1.0))
        V[r, :] .= res.var_genotypic
    end
    binom = [Float64(binomial(big(N), k)) for k = 1:N-1]
    a = sum(binom[k] * (1 / 4k + 1 / 4(N - k)) for k = 1:N-1) / sum(binom)
    hB = hW = 0.0
    E = zeros(T)
    for g = 1:T
        E[g] = σa^2 / 4 * (2hB - hW)
        hB, hW = (1 - a) * hB + a / 2 * hW + 2K * μ, hB + 2K * μ
    end
    for g in (10, 50, 200)
        @test abs(mean(V[:, g]) - E[g]) / (std(V[:, g]) / sqrt(reps)) < 3.5
    end
end

@testset "Heterozygosity decay under weighted parents" begin
    # Fitness that is not heritable (it goes with the position in the
    # cohort) only raises the variance of family size. Two genes in
    # different offspring then derive from the same sire with probability
    # Σpᵢ² over the males, pᵢ = wᵢ/Σw, instead of 1/N_m, so the recursion
    # holds with a = E[Σ_sires pᵢ² + Σ_dams pᵢ²]/4, the expectation being
    # over which individuals are male. N = 12 with exactly 6 males and
    # weights 1, …, 12: a is enumerated over the 924 subsets and gives
    # Ne ≈ 1/a = 9.41, against 12 for equal weights. 400 replicates of 1,000
    # sites on 20 chromosomes of 1 Morgan, 20 generations, no mutation.
    N, nm, chr, reps, t = 12, 6, fill(50, 20), 400, 20
    L = sum(chr)
    w = collect(1.0:N)
    a = let acc = 0.0, k = 0
        for mask = 0:2^N-1
            count_ones(mask) == nm || continue
            male = [(mask >> (i - 1)) & 1 == 1 for i = 1:N]
            pm, pf = w[male] ./ sum(w[male]), w[.!male] ./ sum(w[.!male])
            acc += (sum(abs2, pm) + sum(abs2, pf)) / 4
            k += 1
        end
        acc / k
    end
    @test 1 / a ≈ 9.4147 atol = 1e-3
    function expected(a)
        hB = hW = 0.5
        for _ = 1:t
            hB, hW = (1 - a) * hB + a / 2 * hW, hB
        end
        return hW / 2N + (1 - 1 / N) * hB
    end
    rng = Xoshiro(20261113)
    cbp = UInt32.(cumsum(chr))
    rmap = uniform_recombination_map(chr; M=50)
    H = zeros(reps)
    for r = 1:reps
        start = [UInt32.(findall(rand(rng, Bool, L))) for _ = 1:2N]
        st = SimState(_initial_state(start, 2N, cbp, false)..., N, rng)
        _simulate!(st, ConstantSize(N), RandomMating(nmale=nm), IndexFitness(w), t, rmap,
            Poisson(1.0), UInt32(1):UInt32(L); mutate=false, track=false, fixation_interval=1,
            verbose=false, progress_interval=100, rng=rng)
        H[r] = site_heterozygosity(st.prt[1:2N], L)
    end
    se = sqrt(var(H) / reps)
    @test abs(mean(H) - expected(a)) / se < 3.5
    # The test separates weighted from uniform parents, a = 1/(4·6) + 1/(4·6).
    @test (expected(1 / 12) - mean(H)) / se > 3.5
end

include("properties.jl")
include("golden.jl")
