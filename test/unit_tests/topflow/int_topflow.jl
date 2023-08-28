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

    domain = TopflowDomain(Lx, Ly, nely)
    fea = SimpleTopOpt.TopflowFEA(domain)
    optimizer = OCParameters(200, 0.2)
    dpbc = SimpleTopOpt.DoublePipeBC(domain, fea, Uin)
    dpc = DoublePipeContainer(domain, volfrac, optimizer)

    sol = SimpleTopOpt.TopFlow.topflow(dpc)
    t_xPhys = sol.xPhys

    @test size(r_xPhys) == size(t_xPhys)
    
    residuals = r_xPhys - t_xPhys

    @test norm(residuals) / (30 * 30) ≈ 0 atol=1e-10
    @test sum(abs.(residuals)) / (30 * 30) ≈ 0 atol=1e-10
end


@testset "Pipe Bend Problem; default parameters" begin
    vars = matread("mat_cases/topflow_unit_tests/PipeBendIntegrationCases/paper_params.mat")   
    r_xPhys = vars["xPhys"]

    Lx = 1.0; Ly = 1.0; nely = 30
    volfrac = 1/3 ; Uin = 1e0; rho = 1e0
    mu = 1e0; conit = 50

    domain = TopflowDomain(Lx, Ly, nely)
    feaw = SimpleTopOpt.TopflowFEA(domain)
    optimizer = OCParameters(200, 0.2)
    pbbc = SimpleTopOpt.PipeBendBC(domain, fea, Uin)
    pbc = PipeBendContainer(domain, volfrac, optimizer)

    sol = SimpleTopOpt.TopFlow.topflow(pbc)
    t_xPhys = sol.xPhys

    @test size(r_xPhys) == size(t_xPhys)

    residuals = r_xPhys - t_xPhys

    @test norm(residuals) / (30*30) ≈ 0 atol=1e-10
    @test sum(abs.(residuals)) / (30*30) ≈ 0 atol=1e-10
end