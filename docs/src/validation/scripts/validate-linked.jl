# Linked selection in FisherWright against a matched SLiM model
# (validate-linked.slim; docs/dev/roadmap.md, validation 7).
#
# Model: N = 200 diploids, two sexes with exactly N/2 males (SLiM's WF
# sexual model; RandomMating(nmale = N/2)), soft selection within each sex,
# 4 chromosomes of 1 Mb and 0.1 Morgan (M = 10^7 bp), neutral mutation
# 10^-6 per bp per haplotype per generation, NT = 1,600 generations (8N) from
# a monomorphic start. Selected sites: K = 2,000, in 5 "genes" of 10 kb per
# chromosome (100 sites per gene; genes centred 200 kb apart), mutation rate
# 2.5 × 10^-5 per site per haplotype, so U = 2Kμ = 0.1 per diploid genome.
# Scenarios:
#   neutral: s = 0 at every site (the machinery without selection);
#   bgs:     s ~ -Gamma(0.5, 0.04) per site (mean -0.02), h = 0.25;
#   sweep:   as bgs, but 40 random sites beneficial with s = 0.05, h = 0.5;
#   ben:     the 40 beneficial sites of sweep only, in place;
#   benfree: the same 40 sites, each on a chromosome of its own (1 bp, after
#            the 4 Mb), unlinked from the rest: selection without linkage.
# Neutral statistics are taken over the 4 Mb only, with the distance bins of
# the 2,000 sites in every scenario.
# The coefficient is fixed per site, as in FisherWright; SLiM reuses the
# segregating allele object of a site (mutation() callback), so both run the
# same biallelic model; neutral sites are biallelic in both as well (in SLiM,
# allele objects of a position are pooled and none is converted to a
# substitution). Statistics of the whole population at the end:
# Σ2pq per bp of the neutral mutations (π) and the number of segregating
# neutral sites (S) by distance to the nearest selected site, neutral
# substitutions, segregating and fixed selected sites, Σ of derived
# frequencies, mean fitness and mean log fitness (fixed alleles included),
# and Ne from neutral π. For the π rows the theory column is the neutral
# expectation for biallelic sites (neutral_pi below, a = 1/N); its z is
# shown for the neutral scenario.
#
# Usage:
#   julia --project validate-linked.jl setup DIR
#       writes DIR/sites-<scenario>.tsv and DIR/bins.tsv for SLiM;
#   julia --project -t 16 validate-linked.jl compare DIR R SEED [SCENARIOS]
#       runs R FisherWright replicates per scenario and compares them with
#       DIR/<scenario>.tsv from SLiM.

using FisherWright, Distributions, Random, Statistics, Printf

const N, NT = 200, 1_600
const CL, NCHR = 1_000_000, 4
const CHR = fill(CL, NCHR)
const MBP = 1e7                    # bp per Morgan
const MUN, MUS = 1e-6, 2.5e-5      # neutral per bp, selected per site
const NGENE, GENE, PERGENE = 5, 10_000, 100
const BINS = [0, 1_000, 5_000, 20_000, 50_000]   # lower bounds of distance bins
const BINLAB = ["<1k", "1-5k", "5-20k", "20-50k", ">=50k"]
const SCEN = ["neutral", "bgs", "sweep", "ben", "benfree"]

# Selected positions (1-based), the same for every scenario.
function positions()
    rng = Xoshiro(20261113)
    pos = Int[]
    for c = 0:NCHR-1, g = 0:NGENE-1
        lo = c * CL + 100_000 + 200_000g - GENE ÷ 2
        append!(pos, lo .+ sort!(randperm(rng, GENE)[1:PERGENE]))
    end
    return sort!(pos)
end

function loci(scen)
    pos = positions()
    K = length(pos)
    if scen in ("ben", "benfree")
        _, ben = loci("sweep")
        p = scen == "ben" ? pos[ben] : NCHR * CL .+ (1:count(ben))
        return FitnessLoci(p, 0.05; h=0.5, μ=MUS), trues(count(ben))
    end
    scen == "neutral" && return FitnessLoci(pos, 0.0; h=0.5, μ=MUS), falses(K)
    rng = Xoshiro(7)
    s = max.(-rand(rng, Gamma(0.5, 0.04), K), -1.0)
    h = fill(0.25, K)
    ben = falses(K)
    if scen == "sweep"
        ben[randperm(rng, K)[1:40]] .= true
        s[ben] .= 0.05
        h[ben] .= 0.5
    end
    return FitnessLoci(pos, s; h=h, μ=MUS), ben
end

# Distance bin of every bp (1-based), to the nearest selected site of the same
# chromosome.
function binof()
    pos = positions()
    b = Vector{Int8}(undef, NCHR * CL)
    for c = 0:NCHR-1
        sp = filter(p -> c * CL < p <= (c + 1) * CL, pos)
        k = 1
        for x = c*CL+1:(c+1)*CL
            while k < length(sp) && abs(sp[k+1] - x) <= abs(sp[k] - x)
                k += 1
            end
            b[x] = searchsortedlast(BINS, abs(sp[k] - x))
        end
    end
    return b
end

function setup(dir)
    mkpath(dir)
    for scen in SCEN
        l, ben = loci(scen)
        open(joinpath(dir, "sites-$scen.tsv"), "w") do io
            for k in eachindex(l.pos)
                println(io, l.pos[k] - 1, '\t', l.s[k], '\t', l.h[k], '\t', ben[k] ? 3 : 2)
            end
        end
    end
    b = binof()
    open(joinpath(dir, "bins.tsv"), "w") do io      # runs of one bin, 0-based starts
        for x in eachindex(b)
            (x == 1 || b[x] != b[x-1]) && println(io, x - 1, '\t', b[x])
        end
    end
end

const NB = length(BINS)
const STATS = [["pi/bp " * l for l in BINLAB]; ["S " * l for l in BINLAB];
    "pi/bp all"; "Ne (pi)"; "neutral subs"; "sel segregating"; "sel fixed -";
    "sel fixed +"; "sum sel freq"; "mean w"; "mean log w"]

# Expected Σ2pq per bp at generation NT under neutrality (a = 1/N, two sexes)
# for biallelic sites: ancestral to derived at rate μ per gene, no back
# mutation, so a hit on a gene that carries the derived allele, or on a fixed
# site, changes nothing. For a pair of genes of different individuals (B) or
# of one individual (W), d is the probability that exactly one is derived and
# b that both are; q is the probability that one gene is derived. The
# infinite-sites recursion (every hit a new difference) is 0.2% higher here.
function neutral_pi()
    a, μ = 1 / N, MUN
    q = dB = bB = dW = bW = 0.0
    mut(d, b) = (d * (1 - μ) + 2μ * (1 - μ) * (1 - d - b), b + μ * d + μ^2 * (1 - d - b))
    for _ = 1:NT
        # B pairs: other individuals (1 - a), one gene of one parent (a/2),
        # the two genes of one parent (a/2); W pairs: a sire and a dam.
        nd, nb = (1 - a) * dB + a / 2 * dW, (1 - a) * bB + a / 2 * q + a / 2 * bW
        (dW, bW) = mut(dB, bB)
        (dB, bB) = mut(nd, nb)
        q += μ * (1 - q)
    end
    return dW / 2N + (1 - 1 / N) * dB
end

# Chromosome lengths: the 4 Mb, and one bp per unlinked selected site.
chrs(scen) = scen == "benfree" ? [CHR; ones(Int, length(first(loci(scen)).pos))] : CHR

function fw_stats(res, l::FitnessLoci, ben, bin, nbp)
    L = NCHR * CL
    haps = res.active_haplotypes
    nh = length(haps)
    cnt = Dict{UInt32,Int}()
    for h in haps, p in h
        p <= L && (cnt[p] = get(cnt, p, 0) + 1)
    end
    H, S = zeros(NB), zeros(NB)
    for (p, c) in cnt
        c < nh || continue
        q = c / nh
        b = bin[p]
        H[b] += 2q * (1 - q)
        S[b] += 1
    end
    π = H ./ nbp
    πall = sum(H) / sum(nbp)
    # Selected sites.
    sel, subs = res.selected_haplotypes, res.selected_substitutions
    idx = Dict(UInt32(p) => k for (k, p) in pairs(l.pos))
    scnt = Dict{UInt32,Int}()
    for h in sel, p in h
        scnt[p] = get(scnt, p, 0) + 1
    end
    seg = count(c -> c < nh, values(scnt))
    fixed = vcat(UInt32.(subs), [p for (p, c) in scnt if c == nh])
    fneg = count(p -> !ben[idx[p]], fixed)
    fpos = length(fixed) - fneg
    sumf = sum(values(scnt); init=0) / nh + length(subs)
    lhom = log1p.(l.s)
    lhet = log1p.(l.h .* l.s)
    base = sum(p -> lhom[idx[p]], subs; init=0.0)
    lw = zeros(nh ÷ 2)
    for i = 1:nh÷2
        a, b = sel[2i-1], sel[2i]
        x = base
        for p in a
            x += insorted(p, b) ? lhom[idx[p]] / 2 : lhet[idx[p]]
        end
        for p in b
            x += insorted(p, a) ? lhom[idx[p]] / 2 : lhet[idx[p]]
        end
        lw[i] = x
    end
    return [π; S; πall; N * πall / neutral_pi(); count(<=(L), res.substitutions); seg; fneg; fpos; sumf;
        mean(exp.(lw)); mean(lw)]
end

function run_fw(scen, R, seed)
    l, ben = loci(scen)
    fit = SiteFitness(l)
    bin = binof()
    nbp = [count(==(b), bin) for b = 1:NB]
    out = zeros(R, length(STATS))
    Threads.@threads for r = 1:R
        rng = Xoshiro(hash((seed, r)))
        res = fisher_wright(N, NT, chrs(scen), 1.0; M=MBP, mut_base=1 / MUN, fitness=fit,
            mating=RandomMating(nmale=N ÷ 2), result=true, rng=rng)
        out[r, :] = fw_stats(res, l, ben, bin, nbp)
    end
    return out, nbp
end

function read_slim(path)
    # SLiM reports every statistic but Ne, which follows from π here.
    v = [parse(Float64, x) for l in eachline(path) if startswith(l, "LINKED\t")
         for x in split(l, '\t')[3:end]]
    v = permutedims(reshape(v, length(STATS) + 2, :))
    k = 2NB + 1                                  # pi/bp all
    stats = hcat(v[:, 1:k], N .* v[:, k] ./ neutral_pi(), v[:, k+1:end-3])
    return stats, v[:, end-2:end]
end

function compare(dir, R, seed, only)
    println("FisherWright ", pkgversion(FisherWright), ", Julia ", VERSION, ", threads ", Threads.nthreads())
    println("N = $N ($(N ÷ 2) males), $NT generations, $NCHR chromosomes of $(CL ÷ 1000) kb and ",
        CL / MBP, " Morgan, μ = $MUN per bp, K = $(length(positions())) selected sites, ",
        "μ_sel = $MUS, R = $R replicates per simulator, seed $seed")
    π0 = neutral_pi()
    @printf("neutral expectation of π per bp at generation %d: %.4e (4Nμ = %.4e)\n\n", NT, π0, 4N * MUN)
    @printf("%-8s %-16s %22s %22s %7s %11s %7s\n", "scenario", "stat", "FisherWright", "SLiM", "z", "theory", "z(FW)")
    for scen in SCEN
        only === nothing || scen in only || continue
        fw, nbp = run_fw(scen, R, seed)
        sl, hits = read_slim(joinpath(dir, scen * ".tsv"))
        println("# ", scen, "  (SLiM replicates: ", size(sl, 1), "; bp per distance bin: ", join(nbp, ", "), ")")
        for (k, name) in enumerate(STATS)
            a, b = fw[:, k], sl[:, k]
            ma, mb = mean(a), mean(b)
            sa, sb = std(a) / sqrt(length(a)), std(b) / sqrt(length(b))
            if sa == 0 && sb == 0
                @printf("%-8s %-16s %22s %22s\n", scen, name, @sprintf("%.5g", ma), @sprintf("%.5g", mb))
                continue
            end
            z = (ma - mb) / sqrt(sa^2 + sb^2)
            ispi = startswith(name, "pi/bp")
            t = ispi ? @sprintf("%.4e", π0) : name == "Ne (pi)" ? @sprintf("%.1f", N) : ""
            zt = (ispi || name == "Ne (pi)") && scen == "neutral" ?
                 @sprintf("%.2f", (ma - (ispi ? π0 : N)) / sa) : ""
            @printf("%-8s %-16s %22s %22s %7.2f %11s %7s\n", scen, name,
                @sprintf("%.5g ± %.3g", ma, sa), @sprintf("%.5g ± %.3g", mb, sb), z, t, zt)
        end
        ev, reused, noop = sum(hits[:, 1]), sum(hits[:, 2]), sum(hits[:, 3])
        @printf("%-8s SLiM selected mutation events: %d; at a segregating site %.2f%%, on a carried or fixed site %.2f%%\n",
            scen, ev, 100reused / ev, 100noop / ev)
    end
end

if ARGS[1] == "setup"
    setup(ARGS[2])
else
    compare(ARGS[2], parse(Int, ARGS[3]), parse(Int, ARGS[4]), length(ARGS) > 4 ? split(ARGS[5]) : nothing)
end
