using Test
using Random
using FisherWright: fisher_wright, muts2bitarray, extract_chip_bitarray, merge_sorted, merge_sorted!, recombine,
    cobp, cobp!, FisherWrightResult, to_haplotype, RecombinationMap,
    uniform_recombination_map, random_mate, random_mate!, _fixed_mutations,
    _drop_mutations, _drop_mutations!, _fixation_step, _fixation_step!, _extract_fixed,
    _intersect_sorted!, _remove_sorted!
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
        counts_chr1[k] = count(x -> x < 50_000_000, cbuf)
        counts_chr2[k] = count(x -> 50_000_000 < x < 150_000_000, cbuf)
        assortment[k] = count(x -> x == 50_000_000, cbuf)
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

# END population genetics validation
