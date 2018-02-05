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
    ripser_interface(dist, dim_max, thresh, modulus, eltype(mat))
end

function ripser(mat::UpperTriangular{<:Real};
                dim_max = 1, thresh = Inf, modulus = 2)

    ripser(mat'; dim_max = dim_max, thresh = thresh, modulus = modulus)
end

function ripser_interface(dist::Vector{Float32}, dim_max::Integer,
                          thresh::Real, modulus::Integer, T::Type)
    isprime(modulus) || error("modulus must be a prime number!")
    dim_max ≥ 0      || error("dim_max must be non-negative!")
    thresh > 0       || error("thresh must be positive!")

    libripser = is_windows() ? "libripser.dll" : "libripser.so"
    shlib_path = joinpath(Pkg.dir("Ripser"), "deps", libripser)
    ripser_fptr = Libdl.dlsym(Libdl.dlopen(shlib_path), :ripser)

    origSTDOUT = STDOUT
    (r, w) = redirect_stdout()
    res_str = ""
    try
        ccall(ripser_fptr,
              Void,
              (Cint, Ptr{Float32}, Cint, Float32, Int16),
              length(dist), dist,
              dim_max, thresh, modulus)
    finally
        close(w)
        res_str = String(map(Char, readavailable(r)))
        close(r)
        redirect_stdout(origSTDOUT)
    end

    parse(PersistenceDiagram{T}, res_str)
end

"""
    read_lowertridist(filename)

Read a lower triangular matrix in comma separated format from filename.
Returns `LowerTriangular{Float32}`.
The file should NOT contain the 0s on the diagonal.
"""
function read_lowertridist(filename)
    text = split.(split(readstring(filename), '\n'), ',')
    while text[1] == [""] shift!(text) end
    while text[end] == [""] pop!(text) end

    n = length(text) + 1
    res  = LowerTriangular(zeros(n, n))

    for (i, line) in enumerate(text)
        for (j, entry) in enumerate(line)
            entry == "" && continue
            res[i+1, j] = parse(Float32, entry)
        end
    end

    res
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
