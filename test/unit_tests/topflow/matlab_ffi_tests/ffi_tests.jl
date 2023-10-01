using Test
using SimpleTopOpt
using MAT
using LinearAlgebra
using SimpleTopOpt.TopFlow

@testset "for RES" begin

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

    @test size(t_sR) == (12, 900)
    @test norm(vars["sR"] - t_sR) / (900 * 12) ≈ 0 atol=1e-7
end

@testset "for JAC" begin
    dxv = 0.0333 * ones(1,900)
    dyv = 0.0333 * ones(1,900)
    muv = ones(1, 900)
    rhov = ones(1, 900)
    alpha = 250 * ones(1, 900)
    vars = matread("mat_cases/topflow_unit_tests/MATLAB_FFI/vals/S_of_edofMatT_topFlow_RES_1.mat")
    S = vars["SSS"] 

    t_sJ = SimpleTopOpt.TopFlow.call_JAC(dxv, dyv, muv, rhov, alpha, S)

    vars = matread("mat_cases/topflow_unit_tests/MATLAB_FFI/cases/sJ_JAC_case1.mat")

    @test size(t_sJ) == (144, 900)
    @test norm(vars["sJ"] - t_sJ) / (900 * 144) ≈ 0 atol=1e-7
end

@testset "for PHI" begin
    dxv = 0.0333 * ones(1, 900)
    dyv = 0.0333 * ones(1, 900)
    muv = ones(1, 900)
    alpha = 250 * ones(1, 900)
    vars = matread("mat_cases/topflow_unit_tests/MATLAB_FFI/vals/S_of_edofMatT_topFlow_RES_1.mat")
    S = vars["SSS"] 

    t_obj = sum(SimpleTopOpt.TopFlow.call_PHI(dxv, dyv, muv, alpha, S))
    r_obj = 10.8551

    # NOTE -- currently fails
    @test (t_obj - r_obj) ≈ 0
end

@testset "for dPHIdg" begin
    dxv = 0.0333 * ones(1, 900)
    dyv = 0.0333 * ones(1, 900)
    muv = ones(1, 900)
    rhov = ones(1, 900)
    alpha = 250 * ones(1, 900)
    dalpha = -1.1137e3 * ones(1, 900)
    # NOTE -- S(edofMat') has not changed since the last time this was called
    vars = matread("mat_cases/topflow_unit_tests/MATLAB_FFI/vals/S_of_edofMatT_topFlow_RES_1.mat")
    S = vars["SSS"] 

    t_dphidg = SimpleTopOpt.TopFlow.call_dPHIdg(dxv, dyv, muv, alpha, dalpha, S)

    vars = matread("mat_cases/topflow_unit_tests/MATLAB_FFI/cases/dphidg_dPHIdg_case1.mat")
    r_dphidg = vars["dphidg"]

    @test size(t_dphidg) == size(r_dphidg)
    @test norm(t_dphidg - r_dphidg) / 900 ≈ 0
end


@testset "for dPHIds" begin
    dxv = 0.0333 * ones(1, 900)
    dyv = 0.0333 * ones(1, 900)
    muv = ones(1, 900)
    alpha = 250 * ones(1, 900)
    # NOTE -- S(edofMat') has not changed since the last time this was called
    vars = matread("mat_cases/topflow_unit_tests/MATLAB_FFI/vals/S_of_edofMatT_topFlow_RES_1.mat")
    S = vars["SSS"] 

    t_dphids = SimpleTopOpt.TopFlow.call_dPHIds(dxv, dyv, muv, alpha, S)

    vars = matread("mat_cases/topflow_unit_tests/MATLAB_FFI/cases/ans_dPHIds_case1.mat")
    r_dphids = vars["ans"]

    @test size(t_dphids) == (8, 900)
    @test norm(t_dphids - r_dphids) / (900 * 8) ≈ 0 atol=1e-7
end

@testset "for dRESdg" begin
    dxv = 0.0333 * ones(1, 900)
    dyv = 0.0333 * ones(1, 900)
    muv = ones(1, 900)
    rhov = ones(1, 900)
    alpha = 250 * ones(1, 900)
    dalpha = -1.1137e3 * ones(1, 900)
    # NOTE -- S(edofMat') has not changed since the last time this was called
    vars = matread("mat_cases/topflow_unit_tests/MATLAB_FFI/vals/S_of_edofMatT_topFlow_RES_1.mat")
    S = vars["SSS"] 

    t_dresdg = SimpleTopOpt.TopFlow.call_dRESdg(dxv, dyv, muv, rhov, alpha, dalpha, S)

    vars = matread("mat_cases/topflow_unit_tests/MATLAB_FFI/cases/sR_dRESdg_case1.mat")
    r_dresdg = vars["sR"]

    @test size(r_dresdg) == size(t_dresdg)
    @test norm(r_dresdg - t_dresdg) / (12 * 900) ≈ 0
end