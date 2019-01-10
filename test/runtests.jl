using Ripser
using Test
using SparseArrays
using Distances

@testset "Ripser" begin
    @testset "flatten_distmat" begin
        mat = [1 2 3; 2 1 4; 3 4 1]
        @test Ripser.flatten_distmat(mat) == [2, 3, 4]
        @test eltype(Ripser.flatten_distmat(mat)) == Ripser.Cvalue_t
    end

    @testset "isprime" begin
        @test !Ripser.isprime(1)
        @test Ripser.isprime(2)
        @test Ripser.isprime(3)
        @test !Ripser.isprime(4)
        @test Ripser.isprime(5)
        @test !Ripser.isprime(6)
        @test Ripser.isprime(7)
        @test !Ripser.isprime(8)
        @test !Ripser.isprime(9)
        @test !Ripser.isprime(10)
        @test Ripser.isprime(7867)
        @test !Ripser.isprime(7869)
    end

    @testset "check_args" begin
        mat_square = zeros(3, 3)
        mat_rect = zeros(2, 3)
        @test_throws ArgumentError Ripser.check_args(mat_rect, 2, 1, Inf)
        @test_throws ArgumentError Ripser.check_args(mat_square, -1, 1, Inf)
        @test_throws ArgumentError Ripser.check_args(mat_square, 4, 1, Inf)
        @test_throws ArgumentError Ripser.check_args(mat_square, 2, -1, Inf)
        @test_throws ArgumentError Ripser.check_args(mat_square, 2, 1, 0)
        @test_throws ArgumentError Ripser.check_args(mat_square, 2, 1, -1)
    end

    @testset "ripser dense" begin
        square = Float64[0 1 2 1;
                         1 0 1 2;
                         2 1 0 1;
                         1 2 1 0]
        r = ripser(square)
        @test r[1] == [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, Inf)]
        @test r[2] == [(1.0, 2.0)]

        two_squares = [square fill(Inf, 4, 4); fill(Inf, 4, 4) square]
        r = ripser(two_squares)
        @test r[1] == vcat(fill((0.0, 1.0), 6), fill((0.0, Inf), 2))
        @test r[2] == [(1.0, 2.0), (1.0, 2.0)]

        n = 100
        dim = 3
        circ = Matrix{Float64}(undef, (2, n))
        for i in 1:n
            t = 2rand()
            circ[1, i] = sinpi(t)
            circ[2, i] = cospi(t)
        end
        r, c = ripser(pairwise(Euclidean(), circ), dim_max = dim, cocycles = true)
        @test length(r[1]) == n
        @test all(iszero, first.(r[1]))
        @test r[1][end][2] == Inf
        @test length(r[2]) == 1

        for d in 1:dim+1
            @test length(r[d]) == length(c[d])
        end
    end
end
