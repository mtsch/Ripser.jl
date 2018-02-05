using Ripser
using Base.Test

const example_dir = joinpath(Pkg.dir("Ripser"), "examples")

"""
Run the standalone version of ripser and return PersistenceBarcode.
"""
function runripser(filename; dim_max=1, modulus=2, thresh=Inf)
    origSTDOUT = STDOUT
    (r, w) = redirect_stdout()
    if !isfinite(thresh)
        run(`../deps/ripser.out $filename --dim $dim_max --modulus $modulus`)
    else
        run(`../deps/ripser.out $filename
             --dim $dim_max --modulus $modulus --threshold $thresh`)
    end
    close(w)
    res_str = String(map(Char, readavailable(r)))
    close(r)
    redirect_stdout(origSTDOUT)

    Ripser.parsebarcode(Float64, res_str)
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
        @test Ripser.parsebarcode(Float64, str1) ==
            PersistenceBarcode([PersistencePair(0.0, 0.0)])

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
        @test Ripser.parsebarcode(Float64, str2) ==
            PersistenceBarcode(PersistencePair.([0.0, 1.0, 3.0], [1, 2, Inf]),
                               PersistencePair.([4.0, 0.0, 1.0], [Inf, 4, 1]),
                               PersistencePair.([0.0, 0.0, 1.0], [Inf, 1, Inf]))
    end

    @testset "invalid argument errors" begin
        example_file = joinpath(example_dir, "100", "torus_100.ldm")
        data = read_lowertridist(example_file)

        @test all(Ripser.isprime.([2, 3, 5, 7, 11, 13, 17, 19, 23, 29]))
        @test all(.!Ripser.isprime.([1,  6,  8,  9,  12, 14,
                                     15, 18, 20, 21, 27, 33]))

        @test_throws ErrorException ripser(data, dim_max = -1)
        @test_throws ErrorException ripser(data, modulus = 6)
        @test_throws ErrorException ripser(data, modulus = 15)
        @test_throws ErrorException ripser(data, thresh  = -0.1)
        @test_throws ErrorException ripser(data, thresh  = 0)
        @test_throws ErrorException ripser(data, thresh  = -5)
    end

    @testset "float types" begin
        for T in [Float64, Float32, Float16]
            barcode = ripser(rand(T, 10, 10))
            @test barcode isa PersistenceBarcode{T, Void}
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
