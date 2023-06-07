using Test
using SimpleTopOpt.Top88
using MAT
using LinearAlgebra

"""
Top88 Unit Test Suite
"""

@testset "Single iteration integration tests" begin

    vars = matread("mat_cases/top88_unit_tests/SINGLE_30_30.mat")
    x = top88(30, 30, 0.4, 3.0, 2.0, true, false, 1)
    @test norm(vars["ans"] - x) == 0

    vars = matread("mat_cases/top88_unit_tests/SINGLE_60_40.mat")
    x = top88(60, 40, 0.4, 3.0, 2.0, true, false, 1)
    @test norm(vars["ans"] - x) == 0

end

@testset "Cases for prepare_filter" begin

    varsH = matread("mat_cases/top88_unit_tests/FILTER_20_20_H.mat")
    varsHs = matread("mat_cases/top88_unit_tests/FILTER_20_20_Hs.mat")
    Hh = varsH["H"]; Hsh = varsHs["Hs"]
    H, Hs = prepare_filter(20, 20, 2.0)

    @test norm(Hh - H)   == 0
    @test norm(Hsh - Hs) == 0

end