#======
 Types
======#
struct PersistenceDiagram{T<:AbstractFloat}
    data::Vector{Vector{Tuple{T, T}}}

    PersistenceDiagram{T}(data) where T =
        new(convert(Vector{Vector{Tuple{T, T}}}, data))
end

PersistenceDiagram(data) = PersistenceDiagram{Float64}(data)

Base.getindex(pd::PersistenceDiagram, i) = pd.data[i]
Base.show(io::IO, pd::PersistenceDiagram) =
    print(io, "$(dim(pd))d PersistenceDiagram")
Base.length(pd::PersistenceDiagram) = length(pd.data)

dim(pd::PersistenceDiagram) = length(pd.data) - 1

function Base.print(io::IO, pd::PersistenceDiagram)
    for d in 0:dim(pd)
        println(io, "persistence intervals in dim $d:")
        for pair in pd.data[d+1]
            b = string(pair[1])
            d = isfinite(pair[2]) ? string(pair[2]) : " "
            println(io, "[$b,$d)")
        end
    end
end

function Base.:(==)(pd1::PersistenceDiagram, pd2::PersistenceDiagram)
    length(pd1) == length(pd2) || return false
    for i in 1:length(pd1)
        pd1[i] == pd2[i] || return false
    end
    true
end

function Base.parse(::Type{PersistenceDiagram{T}}, str) where T
    lines = split(str, '\n')
    i = 1
    while !ismatch(r"persistence", lines[i]) i += 1 end

    out = Vector{Tuple{T, T}}[]
    while i <= length(lines) && lines[i] != ""
        if ismatch(r"persistence", lines[i])
            push!(out, [])
        else
            int = parse.(T, matchall(r"[0-9.]+", String(lines[i])))
            length(int) == 1 && push!(int, convert(T, Inf))
            push!(out[end], (int[1], int[2]))
        end
        i += 1
    end
    PersistenceDiagram(out)
end

Base.parse(::Type{PersistenceDiagram}, str) =
    parse(PersistenceDiagram{Float64}, str)

#==============
 Plots recipes
==============#
"""
    getmax(persistencediagram, dim)

Get the max death time in a dimension, ignoring Inf.
"""
function getmax(pd::PersistenceDiagram, d::Int)
    maximum(filter(isfinite, map(x -> x[2], pd.data[d+1])))
end

@recipe function f(pd::PersistenceDiagram; dims=nothing)
    if dims == nothing
        dims = 0:dim(pd)
    end
    if length(dims) == 1
        maxval = getmax(pd, dims)
    else
        maxval = maximum(map(d -> getmax(pd, d), dims))
    end
    infty = round(maxval, RoundUp) +
        maxval > 1 ? length(digits(round(Int, maxval))) : 0
    padding = maxval * 0.1

    xlims := (-padding, maxval)
    ylims := (-padding, infty + padding)

    xlabel := "birth"
    ylabel := "death"

    # x = y line
    @series begin
        seriestype := :path
        label := ""
        color := :black
        [-padding, infty], [-padding, infty]
    end
    # infinity
    @series begin
        seriestype := :path
        label := "infinity"
        color := :grey
        [-padding, infty], [infty, infty]
    end
    # Births and deaths.
    for dim in dims
        @series begin
            seriestype := :scatter
            label := "dim = $dim"

            xs = map(x -> x[1], pd.data[dim+1])
            ys = map(y -> y[2], pd.data[dim+1])
            map!(y -> isfinite(y) ? y : infty, ys, ys)
            xs, ys
        end
    end
end
