using Ripser
using Base.Test
using Plots; unicodeplots()

const example_dir = joinpath(Pkg.dir("Ripser"), "examples")

"""
Run the standalone version of ripser and return PersistenceDiagram.
"""
function runripser(filename; dim_max=1, modulus=2, thresh=Inf)
    origSTDOUT = STDOUT
    (r, w) = redirect_stdout()
    if !isfinite(thresh)
        run(`../deps/ripser.out $filename --dim $dim_max --modulus $modulus`)
    else
        run(`../deps/ripser.out $filename --dim $dim_max --modulus $modulus --threshold $thresh`)
    end
    close(w)
    res_str = String(map(Char, readavailable(r)))
    close(r)
    redirect_stdout(origSTDOUT)

    parse(PersistenceDiagram, res_str)
end

@testset "basics" begin
    @testset "ripser output parsing" begin
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

    @testset "printing" begin
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

    @testset "invalid argument errors" begin
        example_file = joinpath(example_dir, "100", "torus_100.ldm")
        data = read_lowertridist(example_file)

        @test all(Ripser.isprime.([2, 3, 5, 7, 11, 13, 17, 19, 23, 29]))
        @test all(.!Ripser.isprime.([1, 6, 8, 9, 12, 14, 15, 18, 20, 21, 27, 33]))

        @test_throws ErrorException ripser(data, dim_max = -1)
        @test_throws ErrorException ripser(data, modulus = 6)
        @test_throws ErrorException ripser(data, modulus = 15)
        @test_throws ErrorException ripser(data, thresh  = -0.1)
        @test_throws ErrorException ripser(data, thresh  = 0)
        @test_throws ErrorException ripser(data, thresh  = -5)
    end

    @testset "float types" begin
        for T in [Float64, Float32, Float16]
            diagram = ripser(rand(T, 10, 10))
            @test eltype(diagram) == T
            @test typeof(diagram) == PersistenceDiagram{T}
        end
    end

    @testset "matrix types" begin
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

    @testset "plotting" begin
        for f in readdir(joinpath(example_dir, "100"))
            file = joinpath(example_dir, "100", f)
            mat = read_lowertridist(file)
            diagram = ripser(mat)
            @test plot(diagram) ≠ nothing
            @test plot(diagram, dims = 0) ≠ nothing
            @test barcode(diagram) ≠ nothing
            @test barcode(diagram, dims = 0) ≠ nothing
            @test_throws ErrorException plot(diagram, dims = 1:100)
            @test_throws ErrorException plot(diagram, dims = 100)
            @test_throws ErrorException plot(diagram, dims = -1)
            @test_throws ErrorException barcode(diagram, dims = 1:100)
            @test_throws ErrorException barcode(diagram, dims = 100)
            @test_throws ErrorException barcode(diagram, dims = -1)
        end
    end
end

@testset "compare with standalone" begin
    for dir in readdir(example_dir), m in [2,3,7], t in [0.1, 2.5, 5, Inf]
        if dir == "1000"
            dim_max = 1
        elseif dir == "100"
            dim_max = [1, 2]
        else
            dim_max = [1, 2, 3, 4]
        end
        for f in readdir(joinpath(example_dir, dir)), d in dim_max
            file = joinpath(example_dir, dir, f)
            data = read_lowertridist(file)
            @test runripser(file, dim_max = d, modulus = m, thresh  = t) ==
                  ripser(data, dim_max = d, modulus = m, thresh  = t)
        end
    end
end
