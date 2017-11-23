"""
    getmax(persistencediagram, dim)

Get the max death time in a dimension, ignoring Inf.
"""
function getmax(pd::PersistenceDiagram, d::Int)
    maximum(filter(isfinite, map(x -> x[2], pd.data[d+1])))
end

function check_dims(pd, dims)
    if dims == nothing
        dims = 0:dim(pd)
    elseif maximum(dims) > dim(pd)
        error("dims can't be larger than diagram dimension!")
    elseif any(dims .< 0)
        error("dims can't be negative!")
    end
    dims
end

function getinf_maxval(pd::PersistenceDiagram, dims)
    if length(dims) == 1
        maxval = getmax(pd, dims)
    else
        maxval = maximum(map(d -> getmax(pd, d), dims))
    end
    infty = round(maxval, RoundUp) +
        (maxval > 1 ? length(digits(round(Int, maxval))) : 0)

    infty, maxval
end

@recipe function f(pd::PersistenceDiagram; dims=nothing)
    _dims = check_dims(pd, dims)

    infty, maxval = getinf_maxval(pd, _dims)
    padding = maxval * 0.1

    xlims --> (-padding, maxval)
    ylims --> (-padding, infty + padding)

    xlabel := "birth"
    ylabel := "death"

    # Births and deaths.
    for dim in _dims
        @series begin
            seriestype := :scatter
            label := "dim = $dim"

            xs = map(x -> x[1], pd.data[dim+1])
            ys = map(y -> y[2], pd.data[dim+1])
            map!(y -> isfinite(y) ? y : infty, ys, ys)
            xs, ys
        end
    end
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
end

@userplot Barcode

@recipe function f(b::Barcode; dims=nothing)
    if length(b.args) ≠ 1 || !(typeof(first(b.args)) <: PersistenceDiagram)
        error("barcode is expecting a single PersistenceDiagram argument. " *
              "Got: $(typeof(b.args))")
    end
    pd = first(b.args)
    T  = eltype(pd)

    _dims = check_dims(pd, dims)

    infty, maxval = getinf_maxval(pd, _dims)
    xlim --> [0, maxval]

    h = 1
    for dim in _dims
        @series begin
            seriestype := :path
            label := "dim = $dim"
            linewidth --> 3

            xs = T[]; ys = T[]
            for (birth, death) in pd.data[dim+1]
                death = isfinite(death) ? death : infty
                append!(xs, [birth, death, NaN])
                append!(ys, [h, h, NaN])
                h += 1
            end
            xs, ys
        end
    end
    xlim --> [1, h]

    @series begin
        seriestype := :path
        label := "infinity"
        color := :grey
        [infty, infty], [1, h]
    end
end
