using Test
using SimpleTopOpt
using SimpleTopOpt.TopFlow.analyticElement:
    Symbols,
    shapeFunctionDerivatives,
    shapeFunctionsAndMatrices,
    nodalCoordsTransforms,
    nodalDofs,
    stabilizationParameters,
    residualFormation,
    doubleIntegrate,
    formJe,
    computePhi,
    computePartialJeDF,
    computePartialPhiDF,
    computePartialPhiSF,
    jacobianConstruction
using SymbolicUtils
using LinearAlgebra

vars = Symbols()
xv, yv, Np, Nu = shapeFunctionsAndMatrices(vars)
x, y = nodalCoordsTransforms(vars, xv, yv, Np)
iJ, detJ, J = jacobianConstruction(vars, x, y)

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

@testset "Jacobian construction" begin
    @test size(J) == (2, 2)

    @test norm(
        SymbolicUtils.substitute(J, Dict([dx => 1, η => 1, dy => 1, ξ => 1])) -
        [1/2 0; 0 1/2],
    ) == 0
    @test norm(
        SymbolicUtils.substitute(J, Dict([dx => 1, η => 0, dy => 1, ξ => 0])) -
        [1/2 0; 0 1/2],
    ) == 0
    @test norm(
        SymbolicUtils.substitute(J, Dict([dx => 10, dy => 10, ξ => 0, η => 0])) -
        [5 0; 0 5],
    ) == 0

    sdJ = SymbolicUtils.substitute(detJ, Dict([dx => 1, dy => 1]))

    @test sdJ == 1 / 4
end

dNpdx, dNudx = shapeFunctionDerivatives(iJ, Np, Nu, vars)

@testset "Test dNpdx and dNudx" begin
    Np_1 = SymbolicUtils.substitute(dNpdx, Dict([ξ => 1, η => 1, dx => 1, dy => 1]))
    Np_2 = SymbolicUtils.substitute(dNpdx, Dict([ξ => 0, dy => 1, η => 0, dx => 1]))
    Nu_1 = SymbolicUtils.substitute(dNudx, Dict([ξ => 1, η => 1, dx => 1, dy => 1]))
    Nu_2 = SymbolicUtils.substitute(dNudx, Dict([ξ => 0, dy => 1, η => 0, dx => 1]))

    @test norm(Np_1 - [0 0 1 -1; 0 -1 1 0]) == 0
    @test norm(Np_2 - [-1/2 1/2 1/2 -1/2; -1/2 -1/2 1/2 1/2]) == 0

    @test norm(Nu_1[1, :, 1] - [0 0 0 0 1 0 -1 0]') == 0
    @test norm(Nu_1[2, :, 1] - [0 0 0 0 0 1 0 -1]') == 0
    @test norm(Nu_1[1, :, 2] - [0 0 -1 0 1 0 0 0]') == 0
    @test norm(Nu_1[2, :, 2] - [0 0 0 -1 0 1 0 0]') == 0

    @test norm(Nu_2[1, :, 1] - [-1 / 2 0 1 / 2 0 1 / 2 0 -1 / 2 0]') == 0
    @test norm(Nu_2[2, :, 1] - [0 -1 / 2 0 1 / 2 0 1 / 2 0 -1 / 2]') == 0
    @test norm(Nu_2[1, :, 2] - [-1 / 2 0 -1 / 2 0 1 / 2 0 1 / 2 0]') == 0
    @test norm(Nu_2[2, :, 2] - [0 -1 / 2 0 -1 / 2 0 1 / 2 0 1 / 2]') == 0
end

s, ux, px, dudx, dpdx = nodalDofs(vars, Nu, Np, dNudx, dNpdx)

@testset "Test ux, px, dudx, and dpdx" begin
    include("intermediates/ux_pd_dudx_dpdx.jl")
end

τ, ue = stabilizationParameters(vars, ux)

@testset "Stabilization parameters" begin
    include("intermediates/stabilization_params.jl")
end

Ru, Rp = residualFormation(vars, τ, Nu, ux, dNudx, dudx, dpdx, px, Np, dNpdx, detJ)

@testset "Testing simplified Ru and Rp" begin
    @testset "Ru" begin
        include("intermediates/preint_Ru.jl")
    end

    @testset "Rp" begin
        include("intermediates/preint_Rp.jl")
    end
end

Ru = doubleIntegrate(Ru, vars)
Rp = doubleIntegrate(Rp, vars)

# TODO -- test against Integrals.jl and see MATLAB as well; substitute all but integration variable
@testset "Integrated Ru and Rp" begin

    @testset "Ru" begin
        include("intermediates/integrated_Ru.jl")
    end

    @testset "Rp" begin
        include("intermediates/integrated_Rp.jl")
    end
end

Re = [Ru; Rp]
Je = formJe(Re, s)

@testset "Testing a few Je components" begin
    @test size(Re) == (12,)
    @test size(Je) == (12, 12)


    # TEST Re

    @test norm(
        SymbolicUtils.substitute(
            Re,
            Dict([
                u1 => 1,
                u2 => 1,
                u3 => 1,
                u4 => 1,
                u5 => 1,
                u6 => 1,
                u7 => 1,
                u8 => 1,
                p1 => 1,
                p2 => 1,
                p3 => 1,
                p4 => 1,
                μ => 1,
                η => 1,
                ξ => 1,
                α => 1,
                ρ => 1,
                dx => 1,
                dy => 1,
            ]),
        ) - [
            ((41)^(1 / 2) * (270 * (41)^(1 / 2) - 360)) / 14760,
            ((41)^(1 / 2) * (270 * (41)^(1 / 2) - 360)) / 14760,
            -1 / 4,
            3 / 4,
            -((41)^(1 / 2) * (90 * (41)^(1 / 2) - 360)) / 14760,
            -((41)^(1 / 2) * (90 * (41)^(1 / 2) - 360)) / 14760,
            3 / 4,
            -1 / 4,
            -((41)^(1 / 2)) / 41,
            0,
            -((41)^(1 / 2)) / 41,
            0,
        ],
    ) ≈ 0 atol = 1e-12


    # TEST Je

    sje = SymbolicUtils.substitute(
        Je,
        Dict([
            u1 => 1,
            u2 => 1,
            u3 => 1,
            u4 => 1,
            u5 => 1,
            u6 => 1,
            u7 => 1,
            u8 => 1,
            p1 => 1,
            p2 => 1,
            p3 => 1,
            p4 => 1,
            μ => 1,
            η => 1,
            ξ => 1,
            α => 1,
            ρ => 1,
            dx => 1,
            dy => 1,
        ]),
    )
    @test (
        sje[4, 1] - (-(41^(1 / 2) * ((3555 * 41^(1 / 2)) / 41 - 60)) / 14760 - 3 / 328)
    ) ≈ 0 atol = 1e-12
    @test (sje[12, 3] - ((41^(1 / 2) * (6 * 41^(1 / 2) + 6)) / 2952 - 41^(1 / 2) / 492)) ≈ 0 atol =
        1e-12
    @test (sje[5, 9] - (-(41^(1 / 2) * (30 * 41^(1 / 2) + 150)) / 14760)) ≈ 0 atol = 1e-12

    sje = SymbolicUtils.substitute(
        Je,
        Dict([
            u1 => 2,
            u2 => 1,
            u3 => 2,
            u4 => 1,
            u5 => 2,
            u6 => 1,
            u7 => 2,
            u8 => 1,
            p1 => 2,
            p2 => 1,
            p3 => 2,
            p4 => 1,
            μ => 2,
            η => 1,
            ξ => 2,
            α => 1,
            ρ => 2,
            dx => 1,
            dy => 2,
        ]),
    )

    @test (
        sje[11, 1] - (
            (215 * 1001^(1 / 2)) / 858858 -
            (5 * 1001^(1 / 2) * ((12 * 1001^(1 / 2)) / 5 + 132)) / 72072
        )
    ) ≈ 0 atol = 1e-12
    @test (
        sje[5, 7] - (
            (5 * 1001^(1 / 2) * (36 * 1001^(1 / 2) - 4080)) / 9018009 -
            (1001^(1 / 2) * ((257696 * 1001^(1 / 2)) / 1001 + 2580)) / 72072
        )
    ) ≈ 0 atol = 1e-12
    @test (
        sje[1, 8] - (
            (1001^(1 / 2) * ((39636 * 1001^(1 / 2)) / 1001 - 280)) / 72072 -
            (5 * 1001^(1 / 2) * (180 * 1001^(1 / 2) - 3120)) / 18036018
        )
    ) ≈ 0 atol = 1e-12
end

ϕ = computePhi(vars, ux, dudx)

@testset "Testing ϕ pre-integration" begin
    @test SymbolicUtils.substitute(
        ϕ,
        Dict([
            u1 => 1,
            u2 => 1,
            u3 => 1,
            u4 => 1,
            u5 => 1,
            u6 => 1,
            u7 => 1,
            u8 => 1,
            α => 1,
            ξ => 1,
            η => 1,
            μ => 1,
            dx => 1,
            dy => 1,
        ]),
    ) - 1 == 0

    @test SymbolicUtils.substitute(
        ϕ,
        Dict([
            u1 => 2,
            u2 => 2,
            u3 => 2,
            u4 => 2,
            u5 => 2,
            u6 => 1,
            u7 => 1,
            u8 => 1,
            α => 1,
            ξ => 1,
            η => 1,
            μ => 1,
            dx => 1,
            dy => 1,
        ]),
    ) - 9 / 2 == 0
end

ϕ = doubleIntegrate(detJ * ϕ, vars)

@testset "Testing integrated ϕ" begin
    @test SymbolicUtils.substitute(
        ϕ,
        Dict([
            u1 => 2,
            u2 => 2,
            u3 => 2,
            u4 => 2,
            u5 => 2,
            u6 => 1,
            u7 => 1,
            u8 => 1,
            α => 1,
            ξ => 1,
            η => 1,
            μ => 1,
            dx => 1,
            dy => 1,
        ]),
    ) - (38 / 9) ≈ 0 atol = 1e-12

    @test SymbolicUtils.substitute(
        ϕ,
        Dict([
            u1 => 1,
            u2 => 1,
            u3 => 1,
            u4 => 1,
            u5 => 1,
            u6 => 1,
            u7 => 1,
            u8 => 1,
            α => 1,
            ξ => 1,
            η => 1,
            μ => 1,
            dx => 1,
            dy => 1,
        ]),
    ) - 1 ≈ 0 atol = 1e-12

end

dphidg = computePartialPhiDF(vars, ϕ)
dphids = computePartialPhiSF(ϕ, s)
drdg = computePartialJeDF(vars, Re)

@testset "Testing dphidg" begin
    @test SymbolicUtils.substitute(
        dphidg,
        Dict([
            u1 => 2,
            u2 => 2,
            u3 => 2,
            u4 => 2,
            u5 => 2,
            u6 => 1,
            u7 => 1,
            u8 => 1,
            α => 1,
            ξ => 1,
            η => 1,
            μ => 1,
            dx => 1,
            dy => 1,
            dα => 1,
        ]),
    ) - 49 / 18 == 0

    @test SymbolicUtils.substitute(
        dphidg,
        Dict([
            u1 => 1,
            u2 => 1,
            u3 => 1,
            u4 => 1,
            u5 => 1,
            u6 => 1,
            u7 => 1,
            u8 => 1,
            α => 1,
            ξ => 1,
            η => 1,
            μ => 1,
            dx => 1,
            dy => 1,
            dα => 100,
        ]),
    ) - 100 == 0
end
