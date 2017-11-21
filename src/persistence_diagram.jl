struct PersistenceDiagram{T<:AbstractFloat}
    data::Vector{Vector{Tuple{T, T}}}

    PersistenceDiagram{T}(data) where T =
        new(convert(Vector{Vector{Tuple{T, T}}}, data))
end

PersistenceDiagram(data) = PersistenceDiagram{Float64}(data)

Base.getindex(pd::PersistenceDiagram, i) = pd.data[i]
Base.show(io::IO, pd::PersistenceDiagram{T}) where T =
    print(io, "$(dim(pd))d PersistenceDiagram{$T}")
Base.length(pd::PersistenceDiagram) = length(pd.data)
Base.eltype(pd::PersistenceDiagram{T}) where T = T

dim(pd::PersistenceDiagram) = length(pd.data) - 1

function Base.print(io::IO, pd::PersistenceDiagram{T}) where T
    println(io, "PersistenceDiagram{$T}:")
    for d in 0:dim(pd)
        println(io, "  persistence intervals in dim $d:")
        for pair in pd.data[d+1]
            b = string(pair[1])
            d = isfinite(pair[2]) ? string(pair[2]) : " "
            println(io, "   [$b,$d)")
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
    while i ≤ length(lines) && lines[i] ≠ ""
        if ismatch(r"persistence", lines[i])
            push!(out, [])
        else
            int = parse.(T, matchall(r"[0-9.]+", String(lines[i])))
            length(int) == 1 && push!(int, convert(T, Inf))
            push!(out[end], (int[1], int[2]))
        end
        i += 1
    end
    PersistenceDiagram{T}(out)
end

Base.parse(::Type{PersistenceDiagram}, str) =
    parse(PersistenceDiagram{Float64}, str)
