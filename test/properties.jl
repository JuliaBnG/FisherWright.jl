# Property-based tests: each docstring contract is stated for all inputs and
# checked against a deliberately naive reference by Supposition.jl, which
# shrinks any failure to a minimal counterexample.
#
# Properties marked `broken=true` / `@test_broken` are known violations of the
# documented contract. They turn into "Unexpected Pass" errors once fixed, at
# which point the marker should be removed.

using Supposition, Supposition.Data
using DataFrames: nrow
using FisherWright: _fixed_mutations, _with_fixed, _intersect_sorted!

const PROP_RNG_SEED = 20261029
prop_rng() = Xoshiro(PROP_RNG_SEED)

# Supposition 0.3 counts a failing @check in the enclosing testset only on
# Julia >= 1.11; on 1.10 it prints the counterexample and the testset still
# passes. Checks therefore run with record=false and their outcome is asserted.
function prop_passed(report)
    res = something(report.result)
    res isa Supposition.Pass && return true
    @error "Property $(report.description) failed" result = res
    return false
end

# --- Generators --------------------------------------------------------------

sorted_vec(maxpos, maxlen) =
    map(sort!, Data.Vectors(Data.Integers(UInt32(1), UInt32(maxpos)); max_size=maxlen))
sorted_unique_vec(maxpos, maxlen) = map(v -> unique!(sort!(v)), sorted_vec(maxpos, maxlen))
haplotypes(maxpos, maxlen, maxhap) =
    Data.Vectors(sorted_unique_vec(maxpos, maxlen); min_size=1, max_size=maxhap)

# Recombination map as (chromosome lengths, per-chromosome interval lengths and rates).
const chr_lengths = Data.Vectors(Data.Integers(1, 12); min_size=1, max_size=4)

function map_from(lengths, split_seed)
    rng = Xoshiro(split_seed)
    cbp = UInt32.(cumsum(lengths))
    ends, rates = Vector{UInt32}[], Vector{Float64}[]
    start = UInt32(1)
    for e in cbp
        # random strictly increasing interval ends inside [start, e], ending at e
        inner = sort!(unique!(rand(rng, start:e, rand(rng, 0:2))))
        filter!(<(e), inner)
        push!(ends, vcat(inner, e))
        push!(rates, rand(rng, length(ends[end])) .* 4)
        start = e + UInt32(1)
    end
    return RecombinationMap(cbp, ends, rates)
end

# --- Naive references --------------------------------------------------------

ref_dense(muts, pos) = BitMatrix([p in h for p in pos, h in muts])

# --- Sorted-set kernels ------------------------------------------------------

@testset "Properties: sorted-set kernels" begin
    # Docstring allows non-decreasing inputs (duplicates); existing tests only
    # use unique inputs.
    @test prop_passed(@check db=false record=false rng=prop_rng() function merge_is_sorted_union(
        a=sorted_vec(80, 30), b=sorted_vec(80, 30))
        merge_sorted(a, b) == sort!(union(a, b))
    end)
    @test prop_passed(@check db=false record=false rng=prop_rng() function merge_reuses_dirty_buffer(
        a=sorted_vec(80, 30), b=sorted_vec(80, 30), junk=sorted_vec(80, 50))
        dest = copy(junk)
        merge_sorted!(dest, a, b) === dest && dest == sort!(union(a, b))
    end)
    # Both branches (binary probe when 8m < n, linear merge otherwise).
    @test prop_passed(@check db=false record=false rng=prop_rng() function intersect_is_set_intersection(
        a=sorted_unique_vec(300, 10), b=sorted_unique_vec(300, 200))
        d = UInt32[]
        _intersect_sorted!(d, a, b) == sort!(intersect(a, b)) &&
            _intersect_sorted!(d, b, a) == sort!(intersect(a, b))
    end)
end

# --- Fixation bookkeeping ----------------------------------------------------

@testset "Properties: fixation bookkeeping" begin
    # Removing fixed sites and merging them back must be lossless.
    @test prop_passed(@check db=false record=false rng=prop_rng() function fixation_roundtrip(
        haps=haplotypes(60, 15, 8), common=sorted_unique_vec(60, 5))
        haps = [merge_sorted(h, common) for h in haps]
        fixed = reduce(intersect, haps)
        cleaned, subs = _fixation_step(deepcopy(haps), UInt32[])
        inplace = deepcopy(haps)
        subs! = _fixation_step!(inplace, UInt32[])
        subs == sort(fixed) &&
            subs! == subs && inplace == cleaned &&
            all(h -> isempty(intersect(h, subs)), cleaned) &&
            isempty(_fixed_mutations(cleaned)) &&
            _with_fixed(FisherWrightResult(cleaned, UInt32[60], subs)) == haps
    end)
end

# --- Dense export ------------------------------------------------------------

@testset "Properties: dense export" begin
    cbp = UInt32[100, 200, 300]
    # Sizes reach the threaded atomic-bit path (nlc >= 64, nhp >= 16).
    @test prop_passed(@check db=false record=false rng=prop_rng() max_examples=300 function muts2bitarray_matches_reference(
        muts=haplotypes(300, 30, 40))
        pos = sort!(unique!(reduce(vcat, muts; init=UInt32[])))
        xy, lmp = muts2bitarray(muts, cbp; include_fixed=true)
        ref = ref_dense(muts, pos)
        cnt = vec(sum(ref; dims=2))
        poly = (cnt .> 0) .& (cnt .< length(muts))
        xyp, lmpp = muts2bitarray(muts, cbp)
        xy == ref && lmp.pos == pos &&
            lmp.chr == Int8.(searchsortedfirst.(Ref(cbp), pos)) &&
            lmp.frq == Float32.(cnt) ./ length(muts) &&
            all(lmp.ref .!= lmp.alt) &&
            xyp == ref[poly, :] && lmpp.pos == pos[poly]
    end)
    @test prop_passed(@check db=false record=false rng=prop_rng() max_examples=300 function chip_matches_reference(
        muts=haplotypes(300, 30, 40), chip=sorted_unique_vec(300, 80))
        extract_chip_bitarray(muts, chip) == ref_dense(muts, chip)
    end)

    # chr::Int8 is the 1-based chromosome index from cbp: more than 127
    # chromosomes, or positions past the genome end, have no valid index.
    @test_throws ArgumentError muts2bitarray([UInt32[1295], UInt32[]], UInt32.(10:10:1300))
    @test_throws ArgumentError muts2bitarray([UInt32[5, 25], UInt32[]], UInt32[10])
    @test_throws ArgumentError muts2bitarray([UInt32[5], UInt32[]], UInt32[])
    _, lmp127 = muts2bitarray([UInt32[1265], UInt32[]], UInt32.(10:10:1270))
    @test lmp127.chr == Int8[127]
    @test size(muts2bitarray([UInt32[], UInt32[]], UInt32[])[1]) == (0, 2)
end

# --- Recombination -----------------------------------------------------------

# A breakpoint `c` switches parent between positions c-1 and c. Within a
# chromosome spanning [s, e] the meaningful crossover gaps are s+1:e; the
# chromosome-assortment breakpoint for chromosome i < n is cbp[i] + 1.
function valid_breakpoint(c, cbp)
    i = searchsortedfirst(cbp, c)
    i > length(cbp) && return false
    s = i == 1 ? UInt32(1) : cbp[i-1] + UInt32(1)
    return s < c <= cbp[i] || (i > 1 && c == s)  # c == s: assortment point
end

@testset "Properties: crossover sampling" begin
    @test prop_passed(@check db=false record=false rng=prop_rng() function cobp_sorted_within_genome(
        lengths=chr_lengths, seed=Data.Integers{UInt64}())
        rmap = map_from(lengths, seed)
        co = cobp!(UInt32[], rmap; rng=Xoshiro(seed))
        issorted(co) && all(c -> 1 <= c <= rmap.cbp[end], co)
    end)

    # Every crossover must fall in an in-chromosome gap or be an assortment
    # point (regression: the span was s:e-1 instead of s+1:e).
    @test prop_passed(@check db=false record=false rng=prop_rng() function cobp_breakpoints_are_valid_gaps(
        lengths=chr_lengths, seed=Data.Integers{UInt64}())
        rmap = map_from(lengths, seed)
        co = cobp!(UInt32[], rmap; rng=Xoshiro(seed))
        all(c -> valid_breakpoint(c, rmap.cbp), co)
    end)

    # Every gap is reachable, including the one before the last bp of each
    # interval: a 2 bp chromosome recombines; a 1 bp chromosome has no gaps.
    two = RecombinationMap(UInt32[2], [UInt32[2]], [[100.0]])
    rng = Xoshiro(1)
    @test all(1:200) do _
        co = cobp!(UInt32[], two; rng=rng)
        !isempty(co) && all(==(UInt32(2)), co)
    end
    one = RecombinationMap(UInt32[1, 3], [UInt32[1], UInt32[3]], [[100.0], [100.0]])
    @test all(_ -> all(c -> c in UInt32[2, 3], cobp!(UInt32[], one; rng=rng)), 1:200)
    multi = RecombinationMap(UInt32[6], [UInt32[2, 4, 6]], [[50.0, 50.0, 50.0]])
    seen = Set{UInt32}()
    for _ = 1:200
        union!(seen, cobp!(UInt32[], multi; rng=rng))
    end
    @test seen == Set(UInt32(2):UInt32(6))

    # Constructor must reject non-increasing chromosome ends, including the
    # case where prev_end + 1 would wrap at typemax(UInt32).
    @test_throws ArgumentError RecombinationMap(UInt32[typemax(UInt32), 5],
        [UInt32[typemax(UInt32)], UInt32[5]], [[0.1], [0.1]])
    @test_throws ArgumentError RecombinationMap(UInt32[0], [UInt32[0]], [[0.1]])
    @test RecombinationMap(UInt32[typemax(UInt32)],
        [UInt32[typemax(UInt32)]], [[0.1]]) isa RecombinationMap
end

# --- Mating ------------------------------------------------------------------

@testset "Properties: random mating" begin
    @test prop_passed(@check db=false record=false rng=prop_rng() function mates_have_correct_sex(
        n1=Data.Integers(2, 30), n2=Data.Integers(1, 30), seed=Data.Integers{UInt64}())
        pm, sex = random_mate(n1, n2; rng=Xoshiro(seed), max_tries=1_000)
        all(i -> sex[i], pm[:, 1]) && all(i -> !sex[i], pm[:, 2])
    end)

    # fisher_wright accepts ne = 2, where one sex only has probability 1/2 per
    # draw; a 10-try cap used to throw about once per 1024 generations.
    Random.seed!(3)
    res2 = fisher_wright(2, 20_000, [100], 0.02; result=true)
    @test length(res2.active_haplotypes) == 4
end

# --- Whole simulation --------------------------------------------------------

@testset "Properties: simulation output invariants" begin
    @test prop_passed(@check db=false record=false rng=prop_rng() max_examples=40 function simulation_invariants(
        ne=Data.Integers(2, 12), nt=Data.Integers(1, 40),
        lengths=Data.Vectors(Data.Integers(1, 500); min_size=1, max_size=3),
        seed=Data.Integers{UInt64}())
        Random.seed!(seed)
        res = fisher_wright(ne, nt, lengths, 15.0; mut_base=1e3, M=200, result=true)
        G = res.chromosome_ends[end]
        H = res.active_haplotypes
        length(H) == 2ne && G == sum(lengths) &&
            all(h -> issorted(h) && allunique(h) && all(p -> 1 <= p <= G, h), H) &&
            issorted(res.substitutions) && allunique(res.substitutions) &&
            isempty(_fixed_mutations(H))
    end)

    # A recurrent mutation at a substituted site is a no-op in both return
    # modes (regression: with result = true it re-entered as segregating).
    # Saturated toy genome, so such hits are frequent.
    res = fisher_wright(10, 300, [200, 200], 15.0; mut_base=1e3, result=true, rng=Xoshiro(1))
    @test length(res.substitutions) > 100
    @test all(h -> isempty(intersect(h, res.substitutions)), res.active_haplotypes)
    a, _ = fisher_wright(10, 300, [200, 200], 15.0; mut_base=1e3, rng=Xoshiro(1))
    @test _with_fixed(res) == a

    # Same rng state gives the same realization: across return modes, across
    # repeated calls (the first call in a session used to differ), and with a
    # fixation scan every generation or only at the end.
    args = (10, 100, [1_000_000, 1_000_000], 10.0)
    Random.seed!(7)
    a, _ = fisher_wright(args...)
    Random.seed!(7)
    r = fisher_wright(args...; result=true)
    @test !isempty(r.substitutions)  # the comparisons below are not vacuous
    @test _with_fixed(r) == a
    @test _with_fixed(fisher_wright(args...; result=true, fixation_interval=7, rng=Xoshiro(3))) ==
          fisher_wright(args...; rng=Xoshiro(3))[1]
    Random.seed!(7)
    @test fisher_wright(args...)[1] == a

    # ... and does not depend on the thread count.
    other = Threads.nthreads() == 1 ? 3 : 1
    code = "using FisherWright, Random; " *
           "print(hash(fisher_wright($(args[1]), $(args[2]), $(args[3]), $(args[4]); rng=Xoshiro(11))[1]))"
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -t $other -e $code`
    @test parse(UInt, read(cmd, String)) == hash(fisher_wright(args...; rng=Xoshiro(11))[1])
end
