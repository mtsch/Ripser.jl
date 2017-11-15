using Ripser
using Base.Test
# TODO: need more test files.

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

@testset "Compare output with Ripser standalone." begin
    # Compare the outputs of standalone ripser and this wrapper.

    example_dir = "../deps/ripser/examples"

    # Example 2 is broken - ripser thinks there are 15 points.
    for f in readdir(example_dir)[[1,3]]
        file = joinpath(example_dir, f)
        @test runripser(file) == ripser(read_lowertridist(file))
    end

end

@testset "Different matrices." begin
    # Test if using different kinds of matrices returns the same result.

    example_dir = "../deps/ripser/examples"

    @testset "Transposition" begin
        for f in readdir(example_dir)
            file = joinpath(example_dir, f)
            mat = read_lowertridist(file)
            @test ripser(mat) == ripser(mat')
        end
    end

    @testset "Sparse" begin
        for f in readdir(example_dir)
            file = joinpath(example_dir, f)
            mat = read_lowertridist(file)
            @test ripser(mat) == ripser(sparse(mat))
        end
    end

    @testset "Full" begin
        for f in readdir(example_dir)
            file = joinpath(example_dir, f)
            mat = read_lowertridist(file)
            @test ripser(mat) == ripser(full(mat))
        end
    end

    @testset "Symmetric" begin
        for f in readdir(example_dir)
            file = joinpath(example_dir, f)
            mat = read_lowertridist(file)
            @test ripser(mat) == ripser(Symmetric(mat'))
        end
    end
end
