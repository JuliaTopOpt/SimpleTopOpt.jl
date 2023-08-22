using Test
using SimpleTopOpt.TopFlow
using SimpleTopOpt
using LinearAlgebra
using MAT

@testset "Unit tests for Topflow construction with default parameters" begin

    Lx = 1.0
    Ly = 1.0
    nely = 30

    volfrac = 1 / 3

    Uin = 1e0
    rho = 1e0
    mu = 1e0

    conit = 50

    # Create domain
    tfdc = TopflowDomain(Lx, Ly, nely)

    @test tfdc.nelx == 30
    @test_throws AssertionError TopflowDomain(-1.0, 1.0, 30)
    @test_throws AssertionError TopflowDomain(1.0, -1.0, 30)
    @test_throws AssertionError TopflowDomain(1.0, 1.0, 0)

    @testset "Easy parameter checks" begin
        # Brinkman penalization
        bkman_param = BrinkmanPenalizationParameters(mu)
        @test bkman_param.alphamin == 2.5e-04
        @test bkman_param.alphamax == 25000

        # Continuation strategy
        cont = SimpleTopOpt.TopflowContinuation(volfrac, bkman_param, conit)

        @test cont.ainit ≈ 250.0
        @test cont.qinit ≈ 197.0002

        @test size(cont.qavec) == (1, 4)
        @test cont.qanum == 4

        # Skipping TopflowOptNSParams for now...

    end


    ### Finite Element Container tests

    fea = SimpleTopOpt.TopflowFEA(tfdc)

    @testset "Finite Element Container tests" begin

        @test tfdc.dx ≈ 0.0333 atol = 1e-4
        @test tfdc.dy ≈ 0.0333 atol = 1e-4

        @test fea.neltot == 900
        @test fea.doftot == 2883

        @test size(fea.edofMat) == (900, 12)
        @test fea.edofMat[1, 1] isa Int64
        vars = matread("mat_cases/topflow_unit_tests/TopflowFEA/edofMat_standard.mat")

        @test size(fea.iJ) == (129600, 1)
        @test fea.iJ[1, 1] isa Int64
        vars = matread("mat_cases/topflow_unit_tests/TopflowFEA/iJ_standard.mat")
        @test norm(vars["iJ"] - fea.iJ) ≈ 0

        @test size(fea.jJ) == (129600, 1)
        @test fea.jJ[1, 1] isa Int64
        vars = matread("mat_cases/topflow_unit_tests/TopflowFEA/jJ_standard.mat")
        @test norm(vars["jJ"] - fea.jJ) ≈ 0

        @test size(fea.iR) == (10800, 1)
        @test fea.iR[1, 1] isa Int64
        vars = matread("mat_cases/topflow_unit_tests/TopflowFEA/iR_standard.mat")
        @test norm(vars["iR"] - fea.iR) ≈ 0

        @test size(fea.jR) == (10800, 1)
        @test fea.jR[1, 1] isa Int64
        vars = matread("mat_cases/topflow_unit_tests/TopflowFEA/jR_standard.mat")
        @test norm(vars["jR"] - fea.jR) ≈ 0

        @test size(fea.jE) == (12, 900)
        @test fea.jE[1, 1] isa Int64
        vars = matread("mat_cases/topflow_unit_tests/TopflowFEA/jE_standard.mat")
        @test norm(vars["jE"] - fea.jE) ≈ 0
    end

    optimizer = OCParameters(200, 0.2)

    ### Problem 1

    dpbc = SimpleTopOpt.DoublePipeBC(tfdc, fea, Uin)

    @testset "Double Pipe BC construction" begin
        @test size(dpbc.fixedDofs) == (1, 264)
        vars = matread("mat_cases/topflow_unit_tests/DoublePipeBC/fixedDofs_standard.mat")
        @test norm(vars["fixedDofs"] - dpbc.fixedDofs) ≈ 0

        @test size(dpbc.DIR) == (2883, 1)
        vars = matread("mat_cases/topflow_unit_tests/DoublePipeBC/DIR_standard.mat")
        @test norm(vars["DIR"] - dpbc.DIR) ≈ 0
    end

    dpc = DoublePipeContainer(tfdc, volfrac, optimizer, Uin, rho, mu)
    @testset "Douple Pipe Container construction" begin
        @test dpc.Renum ≈ 0.166666666666666

        # TODO -- is there anything else to test with this?
    end

    ### Problem 2

    pbbc = SimpleTopOpt.PipeBendBC(tfdc, fea, Uin)

    @testset "Pipe Bend BC construction" begin
        @test size(pbbc.fixedDofs) == (1, 254)
        vars = matread("mat_cases/topflow_unit_tests/PipeBendBC/fixedDofs_standard.mat")
        @test norm(vars["fixedDofs"] - pbbc.fixedDofs) ≈ 0

        @test size(pbbc.DIR) == (2883, 1)
        vars = matread("mat_cases/topflow_unit_tests/PipeBendBC/DIR_standard.mat")
        @test norm(vars["DIR"] - pbbc.DIR) ≈ 0

        @test pbbc.inletLength == 6.0
    end

    pbc = PipeBendContainer(tfdc, cont, volfrac, optimizer, Uin, rho, mu)
    @testset "Pipe Bend Container construction" begin
        @test pbc.Renum ≈ 0.2

        @test_throws AssertionError PipeBendContainer(tfdc, 1.0, optimizer, Uin, rho, mu)
        @test_throws AssertionError PipeBendContainer(tfdc, -1.0, optimizer, Uin, rho, mu)
        @test_throws AssertionError PipeBendContainer(tfdc, 0.0, optimizer, Uin, rho, mu)

    end





end
