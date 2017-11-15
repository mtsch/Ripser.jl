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

#==============
 Plots recipes
==============#
function getmax(pd::PersistenceDiagram, d::Int)
    maximum(filter(isfinite, map(x -> x[2], pd.data[d+1])))
end

@recipe function f(pd::PersistenceDiagram, dims=1)
    if length(dims) == 1
        infty = getmax(pd, dims)
    else
        infty = maximum(map(d -> getmax(pd, d), dims))
    end
    infty = round(infty, RoundUp) + length(digits(round(Int, infty)))
    margin = infty * 0.1

    lo = -margin
    hi = infty + margin

    xlims := (lo, hi)
    ylims := (lo, hi)

    xlabel := "birth"
    ylabel := "death"

    # x = y line
    @series begin
        seriestype := :path
        label := ""
        color := :black
        [lo, infty], [lo, infty]
    end
    # infinity
    @series begin
        seriestype := :path
        label := "infinity"
        color := :grey
        [lo, infty], [infty, infty]
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
