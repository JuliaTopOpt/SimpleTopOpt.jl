using Test
using SimpleTopOpt.TopFlow.analyticElement:
    Symbols, doubleIntegrate
using SymbolicUtils
using LinearAlgebra
using Integrals
using SymbolicNumericIntegration

vars = Symbols()

ρ = vars.ρ
μ = vars.μ
α = vars.α
dα = vars.dα
ξ = vars.ξ
η = vars.η
dx = vars.dx
dy = vars.dy
u1 = vars.u1
u2 = vars.u2
u3 = vars.u3
u4 = vars.u4
u5 = vars.u5
u6 = vars.u6
u7 = vars.u7
u8 = vars.u8
p1 = vars.p1
p2 = vars.p2
p3 = vars.p3
p4 = vars.p4

function ftoc(expression, var, lower, upper)
    return simplify(
        SymbolicUtils.substitute(expression, Dict([var => upper])) - SymbolicUtils.substitute(expression, Dict([var => lower]))
    )

end


cases = [4, 2.5, π, 2.718, 2.00000000001]


@testset "Single variable" begin

    @testset "Linear" begin
        for i in cases
            f(xi, eta) = eta + xi

            eta = i
            prob_1 = IntegralProblem(f, -2, 2, eta)
            sol_num = solve(prob_1, QuadGKJL()).u

            exp1 = η + ξ
            exp2 = SymbolicUtils.substitute(exp1, Dict([η => i]))
            exp2 = integrate(exp2, ξ)
            @test exp2[2] == 0.0
            @test exp2[3] == 0.0
            expr = exp2[1]

            sol_sym = ftoc(expr, ξ, -2, 2)

            @test sol_num ≈ sol_sym

            # Serves as sanity check on the FToC 
        end
    end

    # NOTE -- a few of these fail due to numerical error eg exp2[3] is
    #   consistently non-zero...
    @testset "to the -1" begin
        for i in cases
            println("{To the -1} -- i is " * string(i))
            f(xi, eta) = 1 / (xi + eta)

            eta = i
            prob_1 = IntegralProblem(f, 0, 4, eta)
            sol_num = solve(prob_1, QuadGKJL()).u

            exp1 = 1 / (η + ξ)
            exp2 = SymbolicUtils.substitute(exp1, Dict([η => i]))
            exp2 = integrate(exp2, ξ)
            # NOTE: 2
            @test exp2[2] ≈ 0.0 atol=1e-11
            @test exp2[3] ≈ 0.0 atol=1e-11
            expr = exp2[1]

            sol_sym = ftoc(expr, ξ, 0, 4)

            @test sol_num ≈ sol_sym
        end
    end

    # NOTE: some of the 
    @testset "sqrts" begin
        for i in cases
            f(xi, eta) = √(xi + eta)
            exp1 = √(η + ξ)

            eta = i
            prob_1 = IntegralProblem(f, 0, 4, eta)
            sol_num = solve(prob_1, QuadGKJL()).u

            exp2 = SymbolicUtils.substitute(exp1, Dict([η => i]))
            exp2 = integrate(exp2, ξ)
            @test exp2[2] ≈ 0.0 atol=1e-12
            @test exp2[3] ≈ 0.0 atol=1e-12
            expr = exp2[1]

            sol_sym = ftoc(expr, ξ, 0, 4)

            @test sol_num ≈ sol_sym
        end
    end


end


@testset "Double integral" begin
    

end