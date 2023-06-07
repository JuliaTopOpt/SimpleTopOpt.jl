using MAT
using SimpleTopOpt
using Test
using LinearAlgebra
using Statistics

"""
Integration test suite for Top88
"""

@testset "60x40 Integration Suite" begin
    # First comparison on 60x40
    vars = matread("mat_cases/top88_60_40_04_3_2_1.mat")
    x1 = vars["x"]
    x1h = top88(60, 40, 0.4, 3.0, 2.0, true, false)

    @test size(x1h) == size(x1)
    ss = size(x1h)
    num_elements = ss[1] * ss[2]

    @test (mean(x1h) - mean(x1)) ≈ 0 atol=1e-10
    @test (norm(x1h - x1))/num_elements ≈ 0 atol=1e-10
    @test (norm(x1h - x1)) == 0
end

@testset "30x30 Integration Suite" begin
    # Second comparison on 30x30
    vars = matread("mat_cases/top88_30_30_04_3_2_1.mat")
    x2 = vars["ans"]
    x2h = top88(30, 30, 0.4, 3.0, 2.0, true, false)

    @test size(x2h) == size(x2)
    ss = size(x2h)
    num_elements = ss[1] * ss[2]

    @test (mean(x2h) - mean(x2)) ≈ 0 atol=1e-10
    @test (norm(x2h-x2))/num_elements ≈ 0 atol=1e-10
    @test (norm(x2h-x2)) == 0
end