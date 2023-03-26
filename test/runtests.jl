using SimpleTopOpt
using Test
using TopOpt

using BenchmarkTools


@testset "Bogus" begin
    @test 1+1 == 2
    @test 2+2 != 2
    @test 1.0 isa Real
end

@testset "ForcingTopOptPrecompileEtc" begin
    # geso.jl example
    E = 1.0 # Young’s modulus
    v = 0.3 # Poisson’s ratio
    f = 1.0; # downward force

    nels = (160, 40)
    problem = HalfMBB(Val{:Linear}, nels, (1.0, 1.0), E, v, f)

    solver = FEASolver(Direct, problem; xmin=0.01, penalty=TopOpt.PowerPenalty(3.0))

    comp = Compliance(solver)
    volfrac = Volume(solver)
    sensfilter = SensFilter(solver; rmin=4.0)
    geso = GESO(comp, volfrac, 0.5, sensfilter)

    x0 = ones(length(solver.vars))
    result = geso(x0)

    println("Made it to the end")
end





