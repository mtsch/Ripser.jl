# Function pointer for ripser.
const shlib_path  = joinpath(Pkg.dir("Ripser"), "deps", "libripser.so")
const ripser_fptr = Libdl.dlsym(Libdl.dlopen(shlib_path), :ripser)

"""
    ripser(mat; dim_max = 1, thresh = Inf, modulus = 2)

Run ripser on given matrix and return a PersistenceDiagram.
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

function read_lowertridist(path)
    mat = (map(x -> x == "" ? 0 : x, readcsv(path)))
    res = zeros(size(mat, 1) + 1, size(mat, 1) + 1)
    res[2:end, 1:end-1] .= mat
    for i in 1:size(res, 1)
        res[i, i] = 1
    end
    LowerTriangular(res)
end
