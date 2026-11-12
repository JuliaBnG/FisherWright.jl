# Pedigree record of a run (docs/dev/roadmap.md, "Results and bookkeeping").
# Individuals are numbered in birth order across the run: the starting
# population has IDs first_id:first_id+N0-1, and every generation's offspring
# take the next N_g IDs. The record keeps the matings of the last `k`
# generations; the parents of the oldest recorded generation are listed as
# founders.

mutable struct PedRecord
    k::Int                                # generations of matings to keep; 0 for none
    gen::Int                              # generation of the current cohort
    first::Int                            # ID of individual 1 of the current cohort
    next::Int                             # next free ID
    bgen::Int                             # the founders of the record: generation,
    bfirst::Int                           #   first ID, size and sex
    bn::Int
    bsex::Vector{Union{Missing,Bool}}
    cgen::Vector{Int}                     # recorded cohorts: generation of birth,
    cfirst::Vector{Int}                   #   first ID, sire and dam of every member,
    csire::Vector{Vector{Int}}            #   and sex (missing until it reproduces)
    cdam::Vector{Vector{Int}}
    csex::Vector{Vector{Union{Missing,Bool}}}
end

PedRecord(n::Int) = PedRecord(0, 0, 1, n + 1, 0, 1, n, fill(missing, n), Int[], Int[], Vector{Int}[], Vector{Int}[],
    Vector{Union{Missing,Bool}}[])

# Start the record for a cohort of `n` individuals of generation `gen` whose
# first ID is `first`, keeping `k` generations of matings.
function _pedigree_start!(p::PedRecord, k::Int, gen::Int, first::Int, n::Int)
    p.k, p.gen, p.first, p.next = k, gen, first, first + n
    p.bgen, p.bfirst, p.bn, p.bsex = gen, first, n, fill(missing, n)
    return p
end

# One generation: the current cohort (st.n parents, sexes in st.sex) has
# produced `no` offspring with the matings in st.pm.
function _record_pedigree!(p::PedRecord, st, no::Int)
    pf = p.first
    if p.k > 0
        sex = Vector{Union{Missing,Bool}}(view(st.sex, 1:st.n))
        isempty(p.cgen) ? (p.bsex = sex) : (p.csex[end] = sex)
        push!(p.cgen, p.gen + 1)
        push!(p.cfirst, p.next)
        push!(p.csire, [pf + st.pm[i, 1] - 1 for i = 1:no])
        push!(p.cdam, [pf + st.pm[i, 2] - 1 for i = 1:no])
        push!(p.csex, fill(missing, no))
        if length(p.cgen) > p.k
            # The oldest recorded cohort becomes the founders of the record.
            p.bgen, p.bfirst = popfirst!(p.cgen), popfirst!(p.cfirst)
            popfirst!(p.csire)
            popfirst!(p.cdam)
            p.bsex = popfirst!(p.csex)
            p.bn = length(p.bsex)
        end
    end
    p.gen += 1
    p.first = p.next
    p.next += no
    return nothing
end

function _pedigree_frame(p::PedRecord)
    p.k == 0 && return _empty_pedigree()
    id = collect(p.bfirst:p.bfirst+p.bn-1)
    sire, dam = zeros(Int, p.bn), zeros(Int, p.bn)
    sex = copy(p.bsex)
    gen = fill(p.bgen, p.bn)
    for j in eachindex(p.cgen)
        n = length(p.csire[j])
        append!(id, p.cfirst[j]:p.cfirst[j]+n-1)
        append!(sire, p.csire[j])
        append!(dam, p.cdam[j])
        append!(sex, p.csex[j])
        append!(gen, fill(p.cgen[j], n))
    end
    return DataFrame(id=id, sire=sire, dam=dam, sex=sex, generation=gen)
end

# The pedigree of a run continued from `init`, appended to that of `init`
# and limited to the last k generations of matings.
function _pedigree_result(p::PedRecord, init, gen0::Int, first0::Int, n0::Int)
    new = _pedigree_frame(p)
    (p.k == 0 || !(init isa FisherWrightResult) || nrow(init.pedigree) == 0) && return new
    # The starting cohort is the last cohort of init: its rows in init gain
    # the sexes recorded now, and are not repeated.
    old = copy(init.pedigree)
    start = (first0 .<= old.id .< first0 + n0)
    if p.bgen == gen0
        newstart = new[new.generation .== gen0, :]
        old[start, :sex] = newstart.sex[old.id[start] .- first0 .+ 1]
        new = new[new.generation .> gen0, :]
    end
    all = vcat(old, new)
    keep = all.generation .>= p.gen - p.k
    all = all[keep, :]
    oldest = all.generation .== minimum(all.generation)
    if minimum(all.generation) == p.gen - p.k
        all[oldest, :sire] .= 0
        all[oldest, :dam] .= 0
    end
    return all
end
