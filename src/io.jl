# A versioned binary file format for FisherWrightResult (docs/dev/roadmap.md,
# 0.10.0). Serialization is tied to the Julia and package versions; this
# format is defined byte by byte, so a file stays readable by later releases,
# each of which keeps a reader for every earlier format version.
#
# Layout (all integers and floats little-endian):
#
#   magic        8 bytes  "FWRESULT"
#   version      UInt32   format version (1)
#   payload      the fields below, in order
#   checksum     UInt32   CRC-32C of the payload
#
# Payload of format version 1, with `vec` = Int64 length then the elements,
# and `ragged` = Int64 count n, Int64 offsets[n + 1] (offsets[1] = 0), then the
# concatenated UInt32 elements:
#
#   active_haplotypes        ragged UInt32
#   chromosome_ends          vec UInt32
#   substitutions            vec UInt32
#   selected_haplotypes      ragged UInt32
#   selected_substitutions   vec UInt32
#   ngen, T                  Int64, Int64
#   census                   ngen + 1 Int64
#   males                    ngen Int64
#   mean_fitness, mean_log_fitness, var_log_fitness   ngen Float64 each
#   mean_genotypic           T × ngen Float64 (column-major)
#   var_genotypic            T × T × ngen Float64
#   mean_phenotypic          T × ngen Float64
#   var_phenotypic           T × T × ngen Float64
#   paternal_pairs, maternal_pairs, selfing           ngen Float64 each
#   related_matings          ngen Int64
#   first_id                 Int64
#   pedigree rows m          Int64
#   id, sire, dam            m Int64 each
#   sex                      m Int8 (1 male, 0 female, -1 missing)
#   generation               m Int64

const _FWR_MAGIC = b"FWRESULT"
const _FWR_VERSION = UInt32(1)

# Bytes are written and read through these, which keep the running checksum.
mutable struct _CRCIO{S<:IO}
    io::S
    crc::UInt32
end

function _put(c::_CRCIO, x::AbstractArray{T}) where {T<:Union{UInt32,Int64,Float64,Int8}}
    isempty(x) && return nothing
    bytes = reinterpret(UInt8, x isa Array ? vec(x) : collect(x))
    c.crc = crc32c(bytes, c.crc)
    write(c.io, bytes)
    return nothing
end
_put(c::_CRCIO, x::T) where {T<:Union{UInt32,Int64,Float64,Int8}} = _put(c, T[x])
_put_vec(c::_CRCIO, v::AbstractVector) = (_put(c, Int64(length(v))); _put(c, v))

function _put_ragged(c::_CRCIO, v::AbstractVector{Vector{UInt32}})
    _put(c, Int64(length(v)))
    off = Vector{Int64}(undef, length(v) + 1)
    off[1] = 0
    for (i, h) in pairs(v)
        off[i+1] = off[i] + length(h)
    end
    _put(c, off)
    for h in v
        _put(c, h)
    end
    return nothing
end

function _get(c::_CRCIO, ::Type{T}, n::Integer) where {T}
    n >= 0 || throw(ArgumentError("corrupt FisherWright result file: negative length"))
    x = Vector{T}(undef, n)
    n == 0 && return x
    bytes = reinterpret(UInt8, x)
    nb = readbytes!(c.io, bytes, length(bytes))
    nb == length(bytes) || throw(EOFError())
    c.crc = crc32c(bytes, c.crc)
    return x
end
_get(c::_CRCIO, ::Type{T}) where {T} = _get(c, T, 1)[1]
_get_vec(c::_CRCIO, ::Type{T}) where {T} = _get(c, T, _get(c, Int64))

function _get_ragged(c::_CRCIO)
    n = _get(c, Int64)
    off = _get(c, Int64, n + 1)
    (off[1] == 0 && issorted(off)) || throw(ArgumentError("corrupt FisherWright result file: offsets"))
    data = _get(c, UInt32, off[end])
    return [data[off[i]+1:off[i+1]] for i = 1:n]
end

"""
    write_result(path::AbstractString, result::FisherWrightResult)
    write_result(io::IO, result::FisherWrightResult)

Write `result` in FisherWright's binary result format, which [`read_result`](@ref) reads back. The format is
versioned and defined byte by byte (little-endian, with a CRC-32C checksum), not tied to the Julia or package
version, so a file written by this release stays readable by later ones. It holds every field of the result:
haplotypes, substitutions, the selected stream, the record of every generation and the pedigree.

A result is written as is; save a long burn-in once and start replicates from it with `init`.

# Examples
```julia
using FisherWright

res = fisher_wright(50, 100, [100_000], 1.0; result = true)
write_result("burnin.fwr", res)
isequal(read_result("burnin.fwr"), res)   # true
```
"""
write_result(path::AbstractString, r::FisherWrightResult) = open(io -> write_result(io, r), path, "w")

function write_result(io::IO, r::FisherWrightResult)
    ENDIAN_BOM == 0x04030201 || error("the result format is little-endian; this platform is not")
    write(io, _FWR_MAGIC)
    write(io, htol(_FWR_VERSION))
    c = _CRCIO(io, UInt32(0))
    _put_ragged(c, r.active_haplotypes)
    _put_vec(c, r.chromosome_ends)
    _put_vec(c, r.substitutions)
    _put_ragged(c, r.selected_haplotypes)
    _put_vec(c, r.selected_substitutions)
    rec = r.records
    ng, T = ngenerations(rec), ntraits(rec)
    _put(c, Int64(ng))
    _put(c, Int64(T))
    _put(c, Vector{Int64}(rec.census))
    _put(c, Vector{Int64}(rec.males))
    _put(c, rec.mean_fitness)
    _put(c, rec.mean_log_fitness)
    _put(c, rec.var_log_fitness)
    _put(c, rec.mean_genotypic)
    _put(c, rec.var_genotypic)
    _put(c, rec.mean_phenotypic)
    _put(c, rec.var_phenotypic)
    _put(c, rec.paternal_pairs)
    _put(c, rec.maternal_pairs)
    _put(c, rec.selfing)
    _put(c, Vector{Int64}(rec.related_matings))
    _put(c, Int64(r.first_id))
    ped = r.pedigree
    _put(c, Int64(nrow(ped)))
    for col in (:id, :sire, :dam)
        _put(c, Vector{Int64}(ped[!, col]))
    end
    _put(c, Int8[ismissing(s) ? Int8(-1) : Int8(s) for s in ped.sex])
    _put(c, Vector{Int64}(ped.generation))
    write(io, htol(c.crc))
    return nothing
end

"""
    read_result(path::AbstractString) -> FisherWrightResult
    read_result(io::IO) -> FisherWrightResult

Read a result written by [`write_result`](@ref), of this or an earlier format version. Throws an `ArgumentError`
if the file is not a FisherWright result, is of a later format version than this release reads, or fails its
checksum.
"""
read_result(path::AbstractString) = open(read_result, path)

function read_result(io::IO)
    ENDIAN_BOM == 0x04030201 || error("the result format is little-endian; this platform is not")
    magic = read(io, length(_FWR_MAGIC))
    magic == _FWR_MAGIC || throw(ArgumentError("not a FisherWright result file"))
    version = ltoh(read(io, UInt32))
    1 <= version <= _FWR_VERSION || throw(ArgumentError(
        "result format version $version; this release of FisherWright reads versions 1 to $(_FWR_VERSION)"))
    c = _CRCIO(io, UInt32(0))
    r = _read_payload(c, Val(Int(version)))
    ltoh(read(io, UInt32)) == c.crc || throw(ArgumentError("FisherWright result file fails its checksum"))
    return r
end

function _read_payload(c::_CRCIO, ::Val{1})
    haps = _get_ragged(c)
    ends = _get_vec(c, UInt32)
    subs = _get_vec(c, UInt32)
    sel = _get_ragged(c)
    ssubs = _get_vec(c, UInt32)
    ng, T = _get(c, Int64), _get(c, Int64)
    (ng >= 0 && T >= 1) || throw(ArgumentError("corrupt FisherWright result file: record size"))
    census = _get(c, Int64, ng + 1)
    males = _get(c, Int64, ng)
    fm, flm, flv = _get(c, Float64, ng), _get(c, Float64, ng), _get(c, Float64, ng)
    gm = reshape(_get(c, Float64, T * ng), T, ng)
    gv = reshape(_get(c, Float64, T * T * ng), T, T, ng)
    pm = reshape(_get(c, Float64, T * ng), T, ng)
    pv = reshape(_get(c, Float64, T * T * ng), T, T, ng)
    pp, mp, sf = _get(c, Float64, ng), _get(c, Float64, ng), _get(c, Float64, ng)
    related = _get(c, Int64, ng)
    first_id = _get(c, Int64)
    m = _get(c, Int64)
    id, sire, dam = _get(c, Int64, m), _get(c, Int64, m), _get(c, Int64, m)
    sex = Union{Missing,Bool}[s < 0 ? missing : s == 1 for s in _get(c, Int8, m)]
    gen = _get(c, Int64, m)
    records = GenerationRecords(census, males, fm, flm, flv, gm, gv, pm, pv, pp, mp, sf, related)
    ped = DataFrame(id=id, sire=sire, dam=dam, sex=sex, generation=gen)
    return FisherWrightResult(haps, ends, subs, sel, ssubs, records, first_id, ped)
end
