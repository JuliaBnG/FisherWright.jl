# Golden-output regression: fixed seeds must reproduce the realizations of
# v0.4.0 exactly, in both return modes. Refactors that are meant to leave the
# default path unchanged (docs/dev/roadmap.md, design principle 1) are checked
# against these digests.
#
# The digest is FNV-1a (64 bit) over the UInt32 words of the output, so it does
# not depend on `Base.hash`, which may change between Julia versions.
#
# The realization itself is not portable across Julia versions in general.
# `randn` and the tail of `randexp` changed between Julia 1.10 and 1.13, and
# Distributions' Poisson sampler uses `randn` for means >= 6 (PoissonADSampler).
# Below 6 it sums `randexp` draws (PoissonCountSampler); a draw in the changed
# tail (> 7.7) ends the sum either way, so those counts agree across versions.
# The portable cases therefore keep every Poisson mean (mutations per
# haplotype, crossovers per chromosome) below 6 and must match on every Julia
# version. One large-mean case is pinned per verified Julia minor version and
# skipped on others.
#
# Checked against v0.4.0 on Julia 1.10.12 and 1.13.3 (Distributions 0.25.130
# and 0.25.132), at 1 and 4 threads.
#
# A deliberate change to the realization must update the digests and say so in
# the changelog. On a mismatch the test prints the new digests.

using Random: Xoshiro

function fnv1a(h::UInt64, w::UInt32)
    for k = 0:3
        h = (h ⊻ UInt64((w >> (8k)) & 0xff)) * 0x00000100000001b3
    end
    return h
end

function golden_digest(haps::Vector{Vector{UInt32}}, tail::Vector{UInt32})
    h = 0xcbf29ce484222325
    for v in (haps..., tail)
        h = fnv1a(h, UInt32(length(v)))
        for w in v
            h = fnv1a(h, w)
        end
    end
    return h
end

golden_digest((haps, cbp)::Tuple) = golden_digest(haps, vcat(cbp, typemax(UInt32)))
golden_digest(r::FisherWrightResult) =
    golden_digest(r.active_haplotypes, vcat(r.chromosome_ends, typemax(UInt32), r.substitutions))

# Each case covers a different part of the generation loop: several
# chromosomes with frequent crossovers and assortment; fixation scans every
# generation and every 7th; a saturated genome where recurrent hits on
# segregating and substituted sites are frequent; ne = 2, where the
# both-sexes rejection loop retries often. Poisson means: 5.0 mutations and
# 2.0 / 3.0 crossovers; 4.8 mutations; 0.5 mutations and 2.0 crossovers.
# The two `result = true` multi-chromosome cases share a digest by design:
# a final scan at g = nt makes the fixation interval invisible in the output.
const GOLDEN_CASES = [
    (:multi_chr_plain, (20, 60, [200_000, 300_000], 2.0), (M=1e5, mut_base=2e5), 0x2f1ad10bc6acd9f2),
    (:multi_chr_result, (20, 60, [200_000, 300_000], 2.0), (M=1e5, mut_base=2e5, result=true), 0xc792203be6071cb1),
    (:multi_chr_interval7, (20, 60, [200_000, 300_000], 2.0), (M=1e5, mut_base=2e5, result=true, fixation_interval=7), 0xc792203be6071cb1),
    (:saturated, (10, 300, [200, 200], 12.0), (mut_base=1e3, result=true), 0xd1d4ffac3d897d1e),
    (:ne2, (2, 500, [1_000], 5.0), (M=500, mut_base=1e4, result=true), 0xd8528809d519b03c),
]

# 50 mutations per haplotype: the PoissonADSampler path. For this seed the
# digest happens to agree on 1.10 and 1.13 (the changed `randn` tail is rare),
# but that is not guaranteed, so it is pinned per minor version.
const GOLDEN_LARGE_MEAN = (
    args=(20, 60, [200_000, 300_000], 10.0),
    kw=(M=1e5, mut_base=1e5, result=true),
    digests=Dict(v"1.10" => 0xa0f94ba96603b0fc, v"1.13" => 0xa0f94ba96603b0fc),
)

golden_run(args, kw, seed) = fisher_wright(args...; kw..., rng=Xoshiro(seed))

const GOLDEN_SEED = 20261102

@testset "Golden output (v0.4.0 realizations)" begin
    got = [golden_digest(golden_run(args, kw, GOLDEN_SEED)) for (_, args, kw, _) in GOLDEN_CASES]
    for ((name, _, _, want), d) in zip(GOLDEN_CASES, got)
        @test (name, d) == (name, want)
    end
    if any(((c, d),) -> c[4] != d, zip(GOLDEN_CASES, got))
        @info "Golden digests on this run" Tuple(c[1] => repr(d) for (c, d) in zip(GOLDEN_CASES, got))
    end

    minor = VersionNumber(VERSION.major, VERSION.minor)
    want = get(GOLDEN_LARGE_MEAN.digests, minor, nothing)
    if want === nothing
        @test_skip "no large-mean golden digest for Julia $minor"
    else
        d = golden_digest(golden_run(GOLDEN_LARGE_MEAN.args, GOLDEN_LARGE_MEAN.kw, GOLDEN_SEED))
        @test d == want
        d == want || @info "Golden large-mean digest on Julia $minor" repr(d)
    end
end
