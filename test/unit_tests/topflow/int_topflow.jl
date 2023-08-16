using Test
using MAT
using LinearAlgebra

using SimpleTopOpt
using SimpleTopOpt.TopFlow


@testset "Double Pipe Problem; default parameters" begin
    vars = matread("mat_cases/topflow_unit_tests/DoublePipeIntegrationCases/paper_params.mat")
    r_xPhys = vars["xPhys"]

    Lx = 1.0; Ly = 1.0; nely = 30
    volfrac = 1/3 ; Uin = 1e0; rho = 1e0
    mu = 1e0; conit = 50

    tfdc = TopflowDomain(Lx, Ly, nely)
    fea = SimpleTopOpt.TopflowFEA(tfdc)
    optimizer = OCParameters(200, 0.2)
    dpbc = SimpleTopOpt.DoublePipeBC(tfdc, fea, Uin)
    dpc = DoublePipeContainer(tfdc, volfrac, optimizer)

    t_xPhys = SimpleTopOpt.TopFlow.topflow(dpc)

    @test size(r_xPhys) == size(t_xPhys)
    @test norm(r_xPhys - t_xPhys) / (30 * 30) ≈ 0
end


@testset "Pipe Bend Problem; default parameters" begin
    # TODO
end