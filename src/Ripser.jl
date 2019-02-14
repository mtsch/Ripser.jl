module Ripser
export ripser

using SparseArrays: AbstractSparseMatrix, findnz

const depsfile = joinpath(dirname(@__DIR__), "deps", "deps.jl")
if isfile(depsfile)
    include("../deps/deps.jl")
else
    error("Ripser in not properly installed. Please run Pkg.build(\"Ripser\").")
end
__init__() = check_deps()

# The value_t type from ripser source code.
const Cvalue_t = Cfloat

# RawResult contains pointers to all output arrays required by ripser.
struct RawResult{T}
    dim_max        ::Int
    n_intervals    ::Ref{Ptr{Cint}}
    births_deaths  ::Ref{Ptr{Cvalue_t}}
    cocycle_length ::Ref{Ptr{Cint}}
    cocycles       ::Ref{Ptr{Cint}}
end

RawResult{T}(dim_max) where T =
    RawResult{T}(dim_max,
                 Ref{Ptr{Cint}}(),
                 Ref{Ptr{Cvalue_t}}(),
                 Ref{Ptr{Cint}}(),
                 Ref{Ptr{Cint}}())

# Converts a Matrix to Vector{Cvalue_t} since ripser expects a flat array as input.
function flatten_distmat(dists)
    dists_flat = Cvalue_t[]
    for i in 1:size(dists, 1)-1, j in i+1:size(dists, 1)
        push!(dists_flat, Cvalue_t(dists[j, i]))
    end
    dists_flat
end

function isprime(n)
    if iseven(n) || n < 2
        n == 2
    else
        p = 3
        q = n / p
        while p ≤ q
            iszero(n % p) && return false
            p += 2
            q = n / p
        end
        true
    end
end

function check_args(dists, modulus, dim_max, threshold)
    size(dists, 1) == size(dists, 2) || throw(ArgumentError("distance matrix must be square"))
    isprime(modulus) || throw(ArgumentError("modulus must be a prime number"))
    dim_max ≥ 0      || throw(ArgumentError("dim_max must be non-negative"))
    threshold > 0    || throw(ArgumentError("threshold must be positive"))
end

function split_cocycle(cocycle, dim)
    c = Tuple{Vector{Int}, Int}[]
    for i in 1:dim+2:length(cocycle)-dim-1
        push!(c, (cocycle[i:i+dim] .+ 1, cocycle[i+dim+1]))
    end
    c
end

# Unpack RawResult{T} to barcode and cocycles (if return_cocycles is true)
function unpack_results(raw::RawResult{T}, return_cocycles) where T
    dim_max = raw.dim_max

    n_intervals = unsafe_wrap(Vector{Cint}, raw.n_intervals[], raw.dim_max + 1,
                              own = true)
    intervals = unsafe_wrap(Matrix{Cvalue_t}, raw.births_deaths[], (2, sum(n_intervals)),
                            own = true)
    cocycle_length = unsafe_wrap(Vector{Cint}, raw.cocycle_length[], sum(n_intervals),
                                 own = true)
    if sum(cocycle_length) > 0
        cocycles_flat = unsafe_wrap(Vector{Cint}, raw.cocycles[], sum(cocycle_length),
                                    own = true)
    else
        cocycles_flat = Cint[]
    end

    if !return_cocycles
        barcodes = Matrix{T}[]
        start = 0
        for int in n_intervals
            push!(barcodes, T.(intervals[:, start+1:start+int]))
            start += int
        end

        map(barcodes) do bc
            collect(vec(reinterpret(Tuple{T, T}, bc)))
        end
    else
        barcodes = Matrix{T}[]
        cocycles = Vector{Vector{Tuple{Vector{Int}, Int}}}[]

        start_bc = 0
        start_cc = 0
        for int in n_intervals
            push!(barcodes, T.(intervals[:, start_bc+1:start_bc+int]))
            push!(cocycles, Int[])
            for i in 1:int
                len = cocycle_length[start_bc + i]
                push!(cocycles[end],
                      split_cocycle(cocycles_flat[start_cc+1:start_cc+len],
                                    length(cocycles)-1))
                start_cc += len
            end
            start_bc += int
        end

        map(barcodes) do bc
            collect(vec(reinterpret(Tuple{T, T}, bc)))
        end, cocycles
    end
end

"""
    ripser(dists; modulus = 2, dim_max = 1, threshold = Inf, cocycles = false)
"""
function ripser(dists     ::AbstractMatrix{T};
                modulus   ::Integer = 2,
                dim_max   ::Integer = 1,
                threshold ::Real = Inf,
                cocycles  ::Bool = false) where T<:AbstractFloat
    check_args(dists, modulus, dim_max, threshold)

    res = RawResult{T}(dim_max)
    dists_flat = flatten_distmat(dists)

    ripser_fptr = Libdl.dlsym(Libdl.dlopen(libripser), :c_rips_dm)
    n_edges = ccall(ripser_fptr,
                    Cint,
                    (Ptr{Ptr{Cint}}, Ptr{Ptr{Cvalue_t}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}},
                     Ptr{Cvalue_t}, Cint,
                     Cint, Cint, Cvalue_t, Cint),
                    res.n_intervals, res.births_deaths, res.cocycle_length, res.cocycles,
                    dists_flat, length(dists_flat),
                    modulus, dim_max, threshold, cocycles)

    unpack_results(res, cocycles)
end

function ripser(dists     ::AbstractSparseMatrix{T};
                modulus   ::Integer = 2,
                dim_max   ::Integer = 1,
                threshold ::Real = Inf,
                cocycles  ::Bool = false) where T<:AbstractFloat
    check_args(dists, modulus, dim_max, threshold)

    J, I, V = findnz(dists)
    I .-= 1
    J .-= 1
    res = RawResult{T}(dim_max)

    ripser_fptr = Libdl.dlsym(Libdl.dlopen(libripser), :c_rips_dm_sparse)
    n_edges = ccall(ripser_fptr,
                    Cint,
                    (Ptr{Ptr{Cint}}, Ptr{Ptr{Cvalue_t}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}},
                     Ptr{Cint}, Ptr{Cint}, Ptr{Cvalue_t}, Cint, Cint,
                     Cint, Cint, Cvalue_t, Cint),
                    res.n_intervals, res.births_deaths, res.cocycle_length, res.cocycles,
                    Cint.(I), Cint.(J), Cvalue_t.(V), length(I), size(dists, 1),
                    modulus, dim_max, threshold, cocycles)

    unpack_results(res, cocycles)
end
end
