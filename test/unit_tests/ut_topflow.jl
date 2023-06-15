using Test
using SimpleTopOpt.TopFlow
using SimpleTopOpt



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
        bkman_param = BrinkmanPenalizationParamters(mu)
        @test bkman_param.alphamin == 2.5e-04
        @test bkman_param.alphamax == 25000

        # Continuation strategy
        cont = SimpleTopOpt.TopflowContinuation(mu, volfrac, bkman_param, conit)

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
        # TODO: exact matching!

        @test size(fea.iJ) == (129600, 1)
        @test fea.iJ[1, 1] isa Int64
        # TODO exact matching!
        @test size(fea.jJ) == (129600, 1)
        @test fea.jJ[1, 1] isa Int64
        # TODO exact matching!
        @test size(fea.iR) == (10800, 1)
        @test fea.iR[1, 1] isa Int64
        # TODO exact matching!
        @test size(fea.jR) == (10800, 1)
        @test fea.jR[1, 1] isa Int64
        # TODO exact matching!
        @test size(fea.jE) == (12, 900) # TODO -- LHS is (10800,1)
        # NB, size(repmat(1:neltot, 1, 12)) == (1, 10800) in MATLAB
        #  so, somehow, this ends up being repmat(1:neltot, 1, 12)'
        @test fea.jE[1, 1] isa Int64
        # TODO exact matching!
    end

    ocp = OCParameters(200, 0.2)

    ### Problem 1
    dpc = DoublePipeContainer(tfdc, volfrac, ocp)

    @testset "Double Pipe Container construction" begin

        @test dpc.Renum ≈ 0.1667 rtol = 1e-4
    end




    ### Problem 2





end
