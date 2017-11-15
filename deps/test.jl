function ripser(dist::Vector{Vector{Float32}}; dim_max = 1,
                thresh = Inf, modulus = 2)

    so_path = joinpath(Pkg.dir("Ripser"), "deps", "ripser-wrapper.so")
    println(so_path)

    origSTDOUT = STDOUT
    (r, w) = redirect_stdout()

    ccall((:ripser, "./ripser-wrapper.so"),
          Void,
          (Cint, Ptr{Ptr{Float32}}, Cint, Float32, Int16),
          length(dist), dist,
          dim_max, Float32(thresh), Int16(modulus))

    close(w)
    res_str = String(map(Char, readavailable(r)))
    close(r)
    redirect_stdout(origSTDOUT)

    parse_output(res_str)
end

function parse_output(str)
    lines = split(str, '\n')
    i = 3
    out = Vector{Tuple{Float64, Float64}}[]
    while lines[i] != ""
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

function readdistmat(path="ripser/examples/projective_plane.lower_distance_matrix")
    conts = readcsv(path)
    out = Vector{Float32}[]

    for i in 1:size(conts, 1)
        push!(out, [])
        for j in 1:i
            push!(out[end], conts[i, j])
        end
    end

    out
end

ripser(readdistmat())
