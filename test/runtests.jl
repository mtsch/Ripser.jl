using Ripser
using Base.Test
using Plots; unicodeplots()

const example_dir = joinpath(Pkg.dir("Ripser"), "deps", "ripser", "examples")

"""
Run the standalone version of ripser and return PersistenceDiagram.
"""
function runripser(filename)
    origSTDOUT = STDOUT
    (r, w) = redirect_stdout()
    run(`../deps/ripser.out $filename`)
    close(w)
    res_str = String(map(Char, readavailable(r)))
    close(r)
    redirect_stdout(origSTDOUT)

    PersistenceDiagram(Ripser.parse_output(res_str))
end

@testset "Ripser output parsing" begin
    str1 = """
           some
           header
           text
           persistence intervals in dim 0:
            [0,0)
           """
    @test Ripser.parse_output(str1) == [[(0, 0)]]

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
    @test Ripser.parse_output(str2) == [[(0, 1),   (1, 2), (3, Inf)],
                                        [(4, Inf), (0, 4), (1, 1)],
                                        [(0, Inf), (0, 1), (1, Inf)]]
end

@testset "Printing" begin
    # Print and parse should not change the diagram.
    example_file = joinpath(example_dir,
                            "projective_plane.lower_distance_matrix")

    pdiag = PersistenceDiagram([[(1.0, Inf)]])
    @test "$pdiag" == "persistence intervals in dim 0:\n[1.0, )\n"

    pdiag2 = ripser(read_lowertridist(example_file))
    str = "$pdiag2"

    @test PersistenceDiagram(Ripser.parse_output(str)) == pdiag2

end

@testset "Compare with standalone" begin
    #TODO: more examples
    # Example 2 is broken - ripser does not parse it correctly.
    for f in readdir(example_dir)[[1,3]]
        file = joinpath(example_dir, f)
        @test runripser(file) == ripser(read_lowertridist(file))
    end
end

@testset "Matrix types" begin
    # Test if using different kinds of matrices returns the same result.
    for f in readdir(example_dir)
        file = joinpath(example_dir, f)
        mat = read_lowertridist(file)
        @test ripser(mat) ==          # LowerTriangular
              ripser(mat') ==         # UpperTriangular
              ripser(sparse(mat)) ==  # SparseMatrixCSC
              ripser(full(mat)) ==    # Matrix
              ripser(Symmetric(mat')) # Symmetric
    end
end

@testset "Plotting does not crash" begin
    for f in readdir(example_dir)
        file = joinpath(example_dir, f)
        mat = read_lowertridist(file)
        @test plot(ripser(mat)) != nothing
    end
end
