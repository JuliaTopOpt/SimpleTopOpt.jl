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
    tf_domain = TopflowDomain(Lx, Ly, nely)

    @test tf_domain.nelx == 30
    @test_throws AssertionError TopflowDomain(-1.0, 1.0, 30)
    @test_throws AssertionError TopflowDomain(1.0, -1.0, 30)
    @test_throws AssertionError TopflowDomain(1.0, 1.0, 0)

    ### Easy parameter checks
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

    ### Finite Element Container tests
    # fea = SimpleTopOpt.TopflowFEA(tf_domain)

    @test tf_domain.dx ≈ 0.0333 atol = 1e-4
    @test tf_domain.dy ≈ 0.0333 atol = 1e-4

    nodx = tf_domain.nelx + 1
    nody = tf_domain.nely + 1
    nodtot = nodx * nody
    neltot = tf_domain.nelx * tf_domain.nely
    doftot = 3 * nodtot

    nodenrs = reshape(1:nodtot, nody, nodx)
    edofVecU = reshape(2 * nodenrs[1:end-1, 1:end-1] .+ 1, neltot, 1)
    edofMatU =
        repeat(edofVecU, 1, 8) + repeat([0 1 (2 * nely .+ [2 3 0 1]) -2 -1], neltot, 1)
    edofVecP = reshape(nodenrs[1:end-1, 1:end-1], neltot, 1)
    edofMatP = repeat(edofVecP, 1, 4) + repeat([1 (nely .+ [2 1]) 0], neltot, 1)

    edofMat = [edofMatU (2 * nodtot .+ edofMatP)]


    @test nodx == 31
    @test nody == 31
    @test nodtot == 961
    @test neltot == 900
    @test doftot == 2883

    # TODO -- export all the massive MATLAB things into here for comparison?
    @test size(nodenrs) == (31, 31)
    @test nodenrs[1, 1] isa Int64
    # TODO: exact matching!
    @test size(edofVecU) == (900, 1)
    @test edofVecU[1, 1] isa Int64
    # TODO: exact matching!
    @test size(edofMatU) == (900, 8)
    @test edofMatU[1, 1] isa Int64
    # TODO: exact matching!
    @test size(edofVecP) == (900, 1)
    @test edofVecP[1, 1] isa Int64
    # TODO: exact matching!
    @test size(edofMatP) == (900, 4)
    @test edofMatP[1, 1] isa Int64
    # TODO: exact matching!
    @test size(edofMat) == (900, 12)
    @test edofMat[1, 1] isa Int64
    # TODO: exact matching!










end
