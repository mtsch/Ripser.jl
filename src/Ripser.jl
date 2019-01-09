module Ripser

using Suppressor
using PersistenceBarcodes

export ripser

using SparseArrays
using Libdl

#=
const libripser = joinpath(pathof(Ripser), "deps",
                           @static Sys.iswindows() ? "libripser.dll" : "libripser.so")
=#
const libripser = joinpath("..", "deps",
                           @static Sys.iswindows() ? "libripser.dll" : "libripser.so")

# The value_t type from ripser source code.
const Cvalue_t = Cfloat

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

# Ripser accepts a one-dimensional array.
# This function converts a Matrix to Vector{Cvalue_t}.
function flatten_distmat(dist)
    dist_flat = Cvalue_t[]
    for i in 1:size(dist, 1)-1, j in i+1:size(dist, 1)
        push!(dist_flat, Cvalue_t(dist[j, i]))
    end
    dist_flat
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

function check_args(dist, modulus, dim_max, threshold)
    size(dist, 1) == size(dist, 2) || throw(ArgumentError("distance matrix must be square"))
    isprime(modulus) || throw(ArgumentError("modulus must be a prime number"))
    dim_max ≥ 0      || throw(ArgumentError("dim_max must be non-negative"))
    threshold > 0    || throw(ArgumentError("threshold must be positive"))
end


function unpackresults(raw::RawResult{T}) where T
    dim_max = raw.dim_max

    n_intervals = unsafe_wrap(Vector{Cint}, raw.n_intervals[], raw.dim_max + 1, own = false)
    @show n_intervals
    intervals = unsafe_wrap(Matrix{Cvalue_t}, raw.births_deaths[], (2, sum(n_intervals)), own = false)
    cocycle_length = unsafe_wrap(Vector{Cint}, raw.cocycle_length[], sum(n_intervals), own = false)
    if sum(cocycle_length) > 0
        cocycles = unsafe_wrap(Vector{Cint}, raw.cocycles[], sum(cocycle_length), own = false)
    else
        cocycles = Cint[]
    end

    barcodes = Matrix{T}[]
    start = 0
    for int in n_intervals
        push!(barcodes, T.(intervals[:, start+1:start+int]))
        start += int
    end

    barcodes
end

function ripser(dist::AbstractMatrix{T};
                modulus = 2,
                dim_max = 1,
                threshold = Inf,
                cocycles = false) where T
    check_args(dist, modulus, dim_max, threshold)

    res = RawResult{T}(dim_max)
    dist_flat = flatten_distmat(dist)

    ripser_fptr = Libdl.dlsym(Libdl.dlopen(libripser), :c_rips_dm)
    n_edges = ccall(ripser_fptr,
                    Cint,
                    (Ptr{Ptr{Cint}}, Ptr{Ptr{Cvalue_t}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}},
                     Ptr{Cvalue_t}, Cint,
                     Cint, Cint, Cvalue_t, Cint),
                    res.n_intervals, res.births_deaths, res.cocycle_length, res.cocycles,
                    dist_flat, length(dist_flat),
                    modulus, dim_max, threshold, cocycles)

    unpackresults(res)
end

function ripser(dist::AbstractSparseMatrix{T},
                modulus = 2,
                dim_max = 1,
                threshold = Inf,
                cocycles = false) where T
    check_args(dist, modulus, dim_max, threshold)

    I, J, V = findnz(dist)
    res = RawResult{T}(dim_max)

    ripser_fptr = Libdl.dlsym(Libdl.dlopen(libripser), :c_rips_dm_sparse)
    n_edges = ccall(ripser_fptr,
                    Cint,
                    (Ptr{Ptr{Cint}}, Ptr{Ptr{Cvalue_t}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}},
                     Ptr{Cint}, Ptr{Cint}, Ptr{Cvalue_t}, Cint, Cint,
                     Cint, Cint, Cvalue_t, Cint),
                    res.n_intervals, res.births_deaths, res.cocycle_length, res.cocycles,
                    I, J, V, length(I), size(dist, 1),
                    modulus, dim_max, threshold, cocycles)

    unpackresults(res)
end

using Distances
using Plots

function ncirc(n)
    φ = 2π * rand(n)
    pts = [cos.(φ)'; sin.(φ)']
    pairwise(Euclidean(), pts)
end

function plotres(r)
    plt = plot([0, 1.5], [0, 1.5])
    for (i,a) in enumerate(r)
        a = a[:, isfinite.(a[2, :])]
        scatter!(plt, a[1,:], a[2,:], label = i)
    end
    plt
end

#S = ncirc(50)

end
