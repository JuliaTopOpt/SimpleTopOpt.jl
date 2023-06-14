using Test
using SimpleTopOpt.TopFlow
using SimpleTopOpt



@testset "Unit tests for Topflow" begin

    # Create domain
    tf_domain = TopflowDomain()

    @test tf_domain.nelx == 30
    @test_throws AssertionError TopflowDomain(-1.0, 1.0, 30)
    @test_throws AssertionError TopflowDomain(1.0, -1.0, 30)
    @test_throws AssertionError TopflowDomain(1.0, 1.0, 0)


    # Easy parameter checks
    bkman_param = BrinkmanPenalizationParamters(1e0)
    @test bkman_param.alphamin == 2.5e-04
    @test bkman_param.alphamax == 25000




end
