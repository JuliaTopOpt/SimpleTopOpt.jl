using Test
using SimpleTopOpt
using MAT
using LinearAlgebra
using SimpleTopOpt.TopFlow: call_RES

@testset "Test RES" begin

    # NOTE -- these are the inputs as of the default run
    dxv = 0.0333 * ones(1,900)
    dyv = 0.0333 * ones(1,900)
    muv = ones(1, 900)
    rhov = ones(1, 900)
    alpha = 250 * ones(1, 900)
    vars = matread("mat_cases/topflow_unit_tests/MATLAB_FFI/vals/S_of_edofMatT_topFlow_RES_1.mat")
    S = vars["SSS"]

    t_sR = SimpleTopOpt.TopFlow.call_RES(dxv, dyv, muv, rhov, alpha, S)

    vars = matread("mat_cases/topflow_unit_tests/MATLAB_FFI/cases/sR_RES_case1.mat")

    @test norm(vars["sR"] - t_sR) ≈ 0
end

@testset "Test JAC" begin
    dxv = 0.0333 * ones(1,900)
    dyv = 0.0333 * ones(1,900)
    muv = ones(1, 900)
    rhov = ones(1, 900)
    alpha = 250 * ones(1, 900)
    vars = matread("mat_cases/topflow_unit_tests/MATLAB_FFI/vals/S_of_edofMatT_topFlow_RES_1.mat")
    S = vars["SSS"] 

    t_sJ = SimpleTopOpt.TopFlow.JAC(dxv, dyv, muv, rhov, alpha, S)

    vars = matread("mat_cases/topflow_unit_tests/MATLAB_FFI/cases/sJ_JAC_case1.mat")

    @test norm(vars["sJ"] - t_sJ) ≈ 0

end

@testset "Test PHI" begin
    


    t_phi = PHI(dxv, dyv, muv, alpha, S)

end

