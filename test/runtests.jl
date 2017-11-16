using Ripser
using Base.Test
using Plots; unicodeplots()

const example_dir = joinpath(Pkg.dir("Ripser"), "examples")

"""
Run the standalone version of ripser and return PersistenceDiagram.
"""
function runripser(filename; dim_max=1)
    origSTDOUT = STDOUT
    (r, w) = redirect_stdout()
    run(`../deps/ripser.out $filename --dim $dim_max`)
    close(w)
    res_str = String(map(Char, readavailable(r)))
    close(r)
    redirect_stdout(origSTDOUT)

    parse(PersistenceDiagram, res_str)
end

@testset "Ripser output parsing" begin
    str1 = """
           some
           header
           text
           persistence intervals in dim 0:
            [0,0)
           """
    @test parse(PersistenceDiagram, str1) ==
        PersistenceDiagram([[(0, 0)]])

    str2 = """
           persistence intervals in dim 0:
            [0,1)
            [1,2)
            [3, )
           persistence intervals in dim 1:
            [4, )
            [0,4)
            [1,1)
           persistence intervals in dim 2:
            [0, )
            [0,1)
            [1, )
           """
    @test parse(PersistenceDiagram, str2) ==
        PersistenceDiagram([[(0, 1),   (1, 2), (3, Inf)],
                            [(4, Inf), (0, 4), (1, 1)],
                            [(0, Inf), (0, 1), (1, Inf)]])
end

@testset "Printing" begin
    # Print and parse should not change the diagram.
    example_file = joinpath(example_dir, "100", "torus_100.ldm")
    pdiag = PersistenceDiagram([[(1.0, Inf)]])
    pdiag2 = ripser(read_lowertridist(example_file))
    @test sprint(print, pdiag) ==
        "PersistenceDiagram{Float64}:\n  persistence intervals in dim 0:\n   [1.0, )\n"
    @test sprint(show, pdiag2) ==
        "1d PersistenceDiagram{Float64}"

    @test parse(PersistenceDiagram, "$pdiag2") == pdiag2
end

@testset "Compare with standalone" begin
    for d in readdir(example_dir)
        if d == "1000"
            dim_max = 1
        elseif d == "100"
            dim_max = 2
        else
            dim_max = 4
        end
        for f in readdir(joinpath(example_dir, d))
            file = joinpath(example_dir, d, f)
            @test runripser(file, dim_max = dim_max) ==
                ripser(read_lowertridist(file), dim_max = dim_max)
        end
    end
end

@testset "Matrix types" begin
    # Test if using different kinds of matrices returns the same result.
    for f in readdir(joinpath(example_dir, "20"))
        file = joinpath(example_dir, "20", f)
        mat = read_lowertridist(file)
        @test ripser(mat) ==          # LowerTriangular
              ripser(mat') ==         # UpperTriangular
              ripser(sparse(mat)) ==  # SparseMatrixCSC
              ripser(full(mat)) ==    # Matrix
              ripser(Symmetric(mat')) # Symmetric
    end
end

@testset "Plotting does not crash" begin
    for f in readdir(joinpath(example_dir, "100"))
        file = joinpath(example_dir, "100", f)
        mat = read_lowertridist(file)
        @test plot(ripser(mat)) != nothing
    end
end
