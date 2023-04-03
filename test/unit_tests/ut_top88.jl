using Test
using SimpleTopOpt.Top88
using MAT
using LinearAlgebra

@testset "Single iteration integration tests" begin
    vars = matread("mat_cases/top88_unit_tests/SINGLE_30_30.mat")
    x,_,_ = top88(30, 30, 0.4, 3.0, 2.0, true, false, 1)
    @test norm(vars["ans"] - x) == 0

    vars = matread("mat_cases/top88_unit_tests/SINGLE_60_40.mat")
    x,_,_ = top88(60, 40, 0.4, 3.0, 2.0, true, false, 1)
    @test norm(vars["ans"] - x) == 0
end