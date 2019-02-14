using Ripser
using Test
using SparseArrays
using Distances

const Float = Int == Int64 ? Float64 : Float32

function circle(n)
    circ = Matrix{Float}(undef, (2, n))
    for i in 1:n
        φ = 2π * rand()
        circ[1, i] = sin(φ)
        circ[2, i] = cos(φ)
    end
    circ
end

function torus(n)
    tor = Matrix{Float}(undef, (3, n))
    for i in 1:n
        ϑ, φ = 2π .* rand(2)
        tor[1, i] = (2 + cos(ϑ))cos(φ)
        tor[2, i] = (2 + cos(ϑ))sin(φ)
        tor[3, i] = sin(ϑ)
    end
    tor
end

function sparsify(M, τ)
    copy(M)
    M[M .> τ] .= 0
    sparse(M)
end

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
        # Simple cases, where we know what kind of barcode to expect.
        square = Float[0 1 2 1;
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

        # Create a circle with radius one, sampled randomly.
        n = 100
        dim = 3
        circ = circle(n)
        r, c = ripser(pairwise(Euclidean(), circ), dim_max = dim, cocycles = true)

        # Zero-dimensional classes start at 0 and exactly one of them is infinite.
        zeroth = r[1]
        @test length(zeroth) == n
        @test all(iszero, first.(zeroth))
        @test all(isfinite, last.(zeroth[1:end-1]))
        @test last(zeroth[end]) == Inf

        # We expect one 1-dimensional class.
        @test length(r[2]) == 1

        # There should be one cocycle for each bar in barcode.
        @test map(length, r) == map(length, c)
        # All simplices in cocycles have value 1.
        for i in 2:dim+1
            for coc in c[i]
                @test all(v -> length(v[1]) == i, coc)
                @test all(v -> v[2] == 1, coc)
            end
        end

        # Create two circles. The barcode should have two long lines.
        dim = 2
        circs = hcat(circle(n), 5circle(n))
        r, c = ripser(pairwise(Euclidean(), circs), modulus = 3, dim_max = dim, cocycles = true)
        lengths = map(x -> x[2] - x[1], r[2])
        @test length(r[1]) == 2n
        @test length(filter(l -> l > 1, lengths)) == 2

        # There should be one cocycle for each bar in barcode.
        @test map(length, r) == map(length, c)
        # All cocycles have value 1 or 2.
        for i in 2:dim+1
            for coc in c[i]
                @test all(v -> length(v[1]) == i, coc)
                @test all(v -> -1 ≤ v[2] ≤ 1, coc)
            end
        end

        # More points, torus.
        n = 1000
        τ = 3.0
        dim = 1
        tor = torus(n)
        M = pairwise(Euclidean(), tor)
        r, c = ripser(M, dim_max = dim, modulus = 29)
    end

    @testset "ripser sparse" begin
        # Circle, few points, higher dimension.
        n = 100
        τ = 1.0
        dim = 3
        circ = circle(n)
        M = pairwise(Euclidean(), circ)
        r, c = ripser(sparsify(M, τ), dim_max = dim, modulus = 5, cocycles = true)

        @test length(r) == dim + 1
        @test length(c) == dim + 1

        @test map(length, r) == map(length, c)

        for i in 2:dim+1
            for coc in c[i]
                @test all(v -> length(v[1]) == i, coc)
                @test all(v -> -2 ≤ v[2] ≤ 2, coc)
            end
        end

        # Torus, more points.
        n = 300
        τ = 3.0
        dim = 1
        tor = torus(n)
        M = pairwise(Euclidean(), tor)
        r = ripser(sparsify(M, τ), dim_max = dim, modulus = 29)

        @test length(r) == dim + 1
        @test length(r[1]) == n
    end
end
