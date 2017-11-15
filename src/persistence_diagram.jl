#======
 Types
======#
struct PersistenceDiagram{T}
    data::Vector{Vector{Tuple{T, T}}}
end

Base.getindex(pd::PersistenceDiagram, i) = pd.data[i]
Base.show(io::IO, pd::PersistenceDiagram) =
    print("$(dim(pd))d PersistenceDiagram")
Base.length(pd::PersistenceDiagram) = length(pd.data)

dim(pd::PersistenceDiagram) = length(pd.data) - 1

function Base.print(pd::PersistenceDiagram)
    for d in 0:dim(pd)
        println("persistence intervals in dim $d:")
        for pair in pd.data[d+1]
            b = string(pair[1])
            d = isfinite(pair[2]) ? string(pair[2]) : " "
            println("[$b,$d)")
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

#==============
 Plots recipes
==============#
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
