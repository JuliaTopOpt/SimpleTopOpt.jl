using Test
using SimpleTopOpt.TopH
using MAT
using LinearAlgebra

@testset "Cases for OC" begin

    vars_1 = matread("mat_cases/toph_unit_tests/OC_20_20_04.mat")
    xnew = vars_1["ans"]
    @test norm(OC(20, 20, ones(20, 20), 0.4, zeros(20,20)) - xnew) == 0

    vars_2 = matread("mat_cases/toph_unit_tests/OC_20_20_08.mat")
    xnew = vars_2["ans"]
    @test norm(OC(20, 20, 2*ones(20, 20), 0.8, zeros(20,20)) - xnew) == 0

    vars_3 = matread("mat_cases/toph_unit_tests/OC_40_40_01.mat")
    xnew = vars_3["ans"]
    @test norm(OC(40, 40, 2*ones(40, 40), 0.1, zeros(40,40)) - xnew) == 0

end

@testset "Cases for check" begin

    vars_1 = matread("mat_cases/toph_unit_tests/CHECK_20_20_ones.mat")
    dnc = vars_1["ans"]
    @test norm(check(20, 20, 0.5, ones(20, 20), ones(20, 20)) - dnc) == 0

    vars_2 = matread("mat_cases/toph_unit_tests/CHECK_40_40_ones.mat")
    dnc = vars_2["ans"]
    @test norm(check(40, 40, 0.5, ones(40, 40), ones(40, 40)) - dnc) == 0


end

@testset "Cases for FE" begin

    vars_20 = matread("mat_cases/toph_unit_tests/FE_20_20.mat")
    U = vars_20["ans"]
    @test norm(FE(20, 20, ones(20, 20), 0.5) - U) == 0

    vars_40 = matread("mat_cases/toph_unit_tests/FE_40_40.mat")
    U = vars_40["ans"]
    @test norm(FE(40, 40, ones(40, 40), 0.5) - U) == 0

    vars_40_1 = matread("mat_cases/toph_unit_tests/FE_40_40_1.mat")
    U = vars_40_1["ans"]
    @test norm(FE(40, 40, ones(40, 40), 1.0) - U) == 0

end
