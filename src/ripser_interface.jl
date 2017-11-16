"""
    ripser(mat; dim_max = 1, thresh = Inf, modulus = 2)

Run Ripser on given matrix and return a PersistenceDiagram.
"""
function ripser(mat::AbstractMatrix{<:Real};
                dim_max = 1, thresh = Inf, modulus = 2)

    dist = Float32[]
    for i in 2:size(mat, 1)
        for j in 1:(i-1)
            push!(dist, mat[i, j])
        end
    end
    ripser_interface(dist, dim_max, thresh, modulus)
end

function ripser(mat::UpperTriangular{<:Real};
                dim_max = 1, thresh = Inf, modulus = 2)

    ripser(mat'; dim_max = dim_max, thresh = thresh, modulus = modulus)
end

function ripser_interface(dist::Vector{Float32}, dim_max = 1,
                          thresh = Inf, modulus = 2)

    shlib_path  = joinpath(Pkg.dir("Ripser"), "deps", "libripser.so")
    ripser_fptr = Libdl.dlsym(Libdl.dlopen(shlib_path), :ripser)

    origSTDOUT = STDOUT
    (r, w) = redirect_stdout()
    ccall(ripser_fptr,
          Void,
          (Cint, Ptr{Float32}, Cint, Float32, Int16),
          length(dist), dist,
          dim_max, thresh, modulus)

    close(w)
    res_str = String(map(Char, readavailable(r)))
    close(r)
    redirect_stdout(origSTDOUT)

    PersistenceDiagram(parse_output(res_str))
end

"""
    parse_output(str)

Parse the string output of Ripser and return a vector of vectors of tuples.
"""
function parse_output(str)
    lines = split(str, '\n')
    i = 1
    while !ismatch(r"persistence", lines[i]) i += 1 end

    out = Vector{Tuple{Float64, Float64}}[]
    while i <= length(lines) && lines[i] != ""
        if ismatch(r"persistence", lines[i])
            push!(out, [])
        else
            int = parse.(Float64, matchall(r"[0-9.]+", String(lines[i])))
            length(int) == 1 && push!(int, Inf)
            push!(out[end], (int[1], int[2]))
        end
        i += 1
    end
    out
end

"""
    read_lowertridist(filename)

Read a lower triangular matrix in comma separated format from filename.
Returns `LowerTriangular{Float32}`.
The file should NOT contain the 1s on the diagonal.
"""
function read_lowertridist(filename)
    text = split.(split(readstring(filename), '\n'), ',')
    if text[1] == [""]
        text = text[2:end]
    end
    n = length(text) + 1
    res  = LowerTriangular(zeros(n, n))

    for (i, line) in enumerate(text)
        for (j, entry) in enumerate(line)
            entry == "" && continue
            res[i+1, j] = parse(Float32, entry)
        end
    end
    res += I

    res
end
