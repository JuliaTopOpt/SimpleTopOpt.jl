using Test
using SimpleTopOpt
using SimpleTopOpt.TopFlow.analyticElement:
    Symbols, shapeFunctionDerivatives, shapeFunctionsAndMatrices,
    nodalCoordsTransforms, nodalDofs, stabilizationParameters,
    residualFormation, doubleIntegrate, formJe, computePhi,
    computePartialJeDF, computePartialPhiDF, computePartialPhiSF,
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

    @test size(s) == (12,)

    @test norm(
        SymbolicUtils.substitute(
            ux,
            Dict([
                u1 => 1,
                u2 => 1,
                u3 => 1,
                u4 => 1,
                u5 => 1,
                u6 => 1,
                u7 => 1,
                u8 => 1,
                ξ => 1,
                η => 1,
            ]),
        ) - [1; 1],
    ) == 0

    @test norm(
        SymbolicUtils.substitute(
            ux,
            Dict([
                u1 => 1,
                u2 => 2,
                u3 => 3,
                u4 => 4,
                u5 => 5,
                u6 => 6,
                u7 => 7,
                u8 => 8,
                ξ => 9,
                η => 10,
            ]),
        ) - [-68; -67],
    ) == 0

    @test norm(
        SymbolicUtils.substitute(
            ux,
            Dict([
                u1 => 2,
                u2 => 4,
                u3 => 6,
                u4 => 8,
                u5 => 10,
                u6 => 12,
                u7 => 14,
                u8 => 16,
                ξ => 18,
                η => 20,
            ]),
        ) - [-640; -638],
    ) == 0

    @test SymbolicUtils.substitute(
        px,
        Dict([p1 => 1, p2 => 1, p3 => 1, p4 => 1, ξ => 1, η => 1]),
    ) - 1 == 0

    @test SymbolicUtils.substitute(
        px,
        Dict([p1 => 1, p2 => 2, p3 => 3, p4 => 4, ξ => 5, η => 6]),
    ) + 6.5 == 0

    @test norm(
        (SymbolicUtils.substitute(
            dudx,
            Dict([
                u1 => 1,
                u2 => 1,
                u3 => 1,
                u4 => 1,
                u5 => 1,
                u6 => 1,
                u7 => 1,
                u8 => 1,
                ξ => 1,
                η => 1,
                dx => 1,
                dy => 1,
            ]),
        )) - [0 0; 0 0],
    ) == 0

    @test norm(
        (SymbolicUtils.substitute(
            dudx,
            Dict([
                u1 => 1,
                u2 => 2,
                u3 => 3,
                u4 => 4,
                u5 => 5,
                u6 => 6,
                u7 => 7,
                u8 => 8,
                ξ => 9,
                η => 10,
                dx => 11,
                dy => 12,
            ]),
        )) - [-20/11 -7/6; -20/11 -7/6],
    ) ≈ 0 atol = 1e-12

    @test norm(
        (SymbolicUtils.substitute(
            dpdx,
            Dict([p1 => 1, p2 => 1, p3 => 1, p4 => 1, ξ => 1, η => 1, dx => 1, dy => 1]),
        )) - [0; 0],
    ) == 0

    @test norm(
        (SymbolicUtils.substitute(
            dpdx,
            Dict([p1 => 1, p2 => 2, p3 => 3, p4 => 4, ξ => 5, η => 6, dx => 7, dy => 8]),
        )) - [-6 / 7; -3 / 8],
    ) == 0
end

τ, ue = stabilizationParameters(vars, ux)

@testset "Stabilization parameters" begin

    @test SymbolicUtils.substitute(
        ue,
        Dict([u1 => 1, u2 => 1, u3 => 1, u4 => 1, u5 => 1, u6 => 1, u7 => 1, u8 => 1]),
    ) - √2 ≈ 0 atol = 1e-12

    @test SymbolicUtils.substitute(
        τ,
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
            ρ => 1,
            μ => 1,
            dx => 1,
            dy => 1,
        ]),
    ) - (1 / √(41)) ≈ 0 atol = 1e-12

    @test SymbolicUtils.substitute(
        τ,
        Dict([
            u1 => 1,
            u2 => 2,
            u3 => 3,
            u4 => 4,
            u5 => 5,
            u6 => 6,
            u7 => 7,
            u8 => 8,
            α => 9,
            ρ => 10,
            μ => 11,
            dx => 12,
            dy => 13,
        ]),
    ) - (3130 / √(13086113)) ≈ 0 atol = 1e-12

    @test SymbolicUtils.substitute(
        τ,
        Dict([
            u1 => 0,
            u2 => 1,
            u3 => 0,
            u4 => 1,
            u5 => 0,
            u6 => 1,
            u7 => 0,
            u8 => 1,
            α => 1,
            ρ => 1,
            μ => 1,
            dx => 1,
            dy => 1,
        ]),
    ) - (sqrt(39))/39 ≈ 0 atol = 1e-12

    @test SymbolicUtils.substitute(
        τ,
        Dict([
            u1 => 1,
            u2 => 0,
            u3 => 1,
            u4 => 0,
            u5 => 1,
            u6 => 0,
            u7 => 1,
            u8 => 0,
            α => 1,
            ρ => 1,
            μ => 1,
            dx => 1,
            dy => 1,
        ]),
    ) - (sqrt(39))/39 ≈ 0 atol = 1e-12
end

Ru, Rp = residualFormation(vars, τ, Nu, ux, dNudx, dudx, dpdx, px, Np, dNpdx, detJ)

@testset "Testing simplified Ru and Rp" begin
    @testset "Ru" begin
        @test norm(
            SymbolicUtils.substitute(
                Ru,
                Dict([
                    u1 => 2,
                    u2 => 2,
                    u3 => 2,
                    u4 => 2,
                    u5 => 2,
                    u6 => 2,
                    u7 => 2,
                    u8 => 2,
                    p1 => 2,
                    p2 => 2,
                    p3 => 2,
                    p4 => 2,
                    μ => 2,
                    η => 2,
                    ξ => 2,
                    α => 2,
                    ρ => 2,
                    dx => 2,
                    dy => 2,
                ]),
            ) - [
                (8 * sqrt(29)) / 29 + 1 / 2,
                (8 * sqrt(29)) / 29 + 1 / 2,
                -(16 * sqrt(29)) / 29 - 5 / 2,
                -(16 * sqrt(29)) / 29 - 3 / 2,
                (24 * sqrt(29)) / 29 + 15 / 2,
                (24 * sqrt(29)) / 29 + 15 / 2,
                -(16 * sqrt(29)) / 29 - 3 / 2,
                -(16 * sqrt(29)) / 29 - 5 / 2,
            ],
        ) ≈ 0 atol = 1e-12


        # NOTE: Fails with 1e-8 instead of 1e-12
        @test norm(
            SymbolicUtils.substitute(
                Ru,
                Dict([
                    u1 => 1,
                    u2 => 2,
                    u3 => 3,
                    u4 => 4,
                    u5 => 5,
                    u6 => 6,
                    u7 => 7,
                    u8 => 8,
                    p1 => 9,
                    p2 => 10,
                    p3 => 11,
                    p4 => 12,
                    μ => 13,
                    η => 14,
                    ξ => 15,
                    α => 16,
                    ρ => 17,
                    dx => 18,
                    dy => 19,
                ]),
            ) - [
                63614347069/2736 - (15084666695*7976825^(1/2)*8977188^(1/2))/2046798864,
                15946299869/684 - (45375872405*7976825^(1/2)*8977188^(1/2))/6140396592,
                (16169904317*7976825^(1/2)*8977188^(1/2))/2046798864 - 72701627557/2736,
                (48640353143*7976825^(1/2)*8977188^(1/2))/6140396592 - 9112228693/342,
                27962356445/912 - (524906207*7976825^(1/2)*8977188^(1/2))/62024208,
                3504674981/114 - (1578959453*7976825^(1/2)*8977188^(1/2))/186072624,
                (5412222403*7976825^(1/2)*8977188^(1/2))/682266288 - 24467224565/912,
                (16280393737*7976825^(1/2)*8977188^(1/2))/2046798864 - 6133142623/228,
            ],
        ) ≈ 0 atol = 1e-12

        @test norm(
            SymbolicUtils.substitute(
                Ru,
                Dict([
                    u1 => 19,
                    u2 => 18,
                    u3 => 17,
                    u4 => 16,
                    u5 => 15,
                    u6 => 14,
                    u7 => 13,
                    u8 => 12,
                    p1 => 11,
                    p2 => 10,
                    p3 => 11,
                    p4 => 12,
                    μ => 13,
                    η => 14,
                    ξ => 15,
                    α => 16,
                    ρ => 17,
                    dx => 18,
                    dy => 19,
                ]),
            ) - [
                (2210646225*7976825^(1/2)*29472388^(1/2))/407254816 + 138234168527/2736,
                (12143675925*7976825^(1/2)*29472388^(1/2))/2239901488 + 34515858955/684,
              - (2369731999*7976825^(1/2)*29472388^(1/2))/407254816 - 157981545527/2736,
              - (13017576987*7976825^(1/2)*29472388^(1/2))/2239901488 - 19723290635/342,
                  (2538508277*7976825^(1/2)*29472388^(1/2))/407254816 + 60761940895/912,
                 (13944710601*7976825^(1/2)*29472388^(1/2))/2239901488 + 7585865179/114,
                - (2379422503*7976825^(1/2)*29472388^(1/2))/407254816 - 53166819895/912,
              - (13070809539*7976825^(1/2)*29472388^(1/2))/2239901488 - 13275302657/228,
            ],
        ) ≈ 0 atol = 1e-12


        @test norm(
            SymbolicUtils.substitute(
                Ru,
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
                0,
                0,
                -(sqrt(41)) / 164,
                1 / 4 - sqrt(41) / 164,
                (sqrt(41)) / 82,
                (sqrt(41)) / 82,
                (1 / 4) - (sqrt(41) / 164),
                -sqrt(41) / 164,
            ],
        ) ≈ 0 atol = 1e-12

        @test norm(
            SymbolicUtils.substitute(
                Ru,
                Dict([
                    u1 => 0,
                    u2 => 0,
                    u3 => 1,
                    u4 => 1,
                    u5 => 1,
                    u6 => 1,
                    u7 => 0,
                    u8 => 0,
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
                0,
                0,
                - 38^(1/2)/76 - 1/4,
                1/4 - 38^(1/2)/76,
                38^(1/2)/38 + 1,
                38^(1/2)/38 + 1/2,
                - 38^(1/2)/76 - 1/4,
                - 38^(1/2)/76 - 1/4
            ],
        ) ≈ 0 atol = 1e-12

        @test norm(
            SymbolicUtils.substitute(
                Ru,
                Dict([
                    u1 => 2,
                    u2 => 2,
                    u3 => 1,
                    u4 => 1,
                    u5 => 1,
                    u6 => 1,
                    u7 => 2,
                    u8 => 2,
                    p1 => 1,
                    p2 => 1,
                    p3 => 1,
                    p4 => 0,
                    μ => 0,
                    η => 0,
                    ξ => 0,
                    α => 1,
                    ρ => 1,
                    dx => 1,
                    dy => 1,
                ]),
            ) - [
                3/32 - (3*10^(1/2))/160,
                (3*10^(1/2))/160 + 3/32,
                                  -3/32,
                                   3/32,
                (3*10^(1/2))/160 - 3/32,
              - (3*10^(1/2))/160 - 3/32,
                                   3/32,
                                  -3/32,
            ],
        ) ≈ 0 atol = 1e-12

        @test norm(
            SymbolicUtils.substitute(
                Ru,
                Dict([
                    u1 => 2,
                    u2 => 2,
                    u3 => 0,
                    u4 => 0,
                    u5 => 0,
                    u6 => 0,
                    u7 => 2,
                    u8 => 2,
                    p1 => 0,
                    p2 => 0,
                    p3 => 0,
                    p4 => 0,
                    μ => 0,
                    η => 0,
                    ξ => 0,
                    α => 1,
                    ρ => 1,
                    dx => 1,
                    dy => 1,
                ]),
            ) - [
                5^(1/2)/20 - 1/16,
                5^(1/2)/20 - 1/16,
                            -1/16,
                            -1/16,
              - 5^(1/2)/20 - 1/16,
              - 5^(1/2)/20 - 1/16,
                            -1/16,
                            -1/16,
            ],
        ) ≈ 0 atol = 1e-12

        @test norm(
            SymbolicUtils.substitute(
                Ru,
                Dict([
                    u1 => 1,
                    u2 => 1,
                    u3 => 1,
                    u4 => 1,
                    u5 => 1,
                    u6 => 1,
                    u7 => 1,
                    u8 => 1,
                    p1 => 0,
                    p2 => 0,
                    p3 => 0,
                    p4 => 0,
                    μ => 0,
                    η => 0,
                    ξ => 0,
                    α => 1,
                    ρ => 1,
                    dx => 1,
                    dy => 1,
                ]),
            ) - [
                1/16 - 5^(1/2)/20,
                1/16 - 5^(1/2)/20,
                            1/16,
                            1/16,
                5^(1/2)/20 + 1/16,
                5^(1/2)/20 + 1/16,
                            1/16,
                            1/16
            ],
        ) ≈ 0 atol = 1e-12

        @test norm(
            SymbolicUtils.substitute(
                Ru,
                Dict([
                    u1 => 0,
                    u2 => 0,
                    u3 => 0,
                    u4 => 0,
                    u5 => 0,
                    u6 => 0,
                    u7 => 0,
                    u8 => 0,
                    p1 => 1,
                    p2 => 1,
                    p3 => 1,
                    p4 => 1,
                    μ => 0,
                    η => 0,
                    ξ => 0,
                    α => 1,
                    ρ => 1,
                    dx => 1,
                    dy => 1,
                ]),
            ) - [
                1/8,
                1/8,
               -1/8,
                1/8,
               -1/8,
               -1/8,
                1/8,
               -1/8
            ],
        ) ≈ 0 atol = 1e-12


        @test norm(
            SymbolicUtils.substitute(
                Ru,
                Dict([
                    u1 => 1.1,
                    u2 => 1.2,
                    u3 => 1.3,
                    u4 => 1.4,
                    u5 => 1.5,
                    u6 => 1.6,
                    u7 => 1.7,
                    u8 => 1.8,
                    p1 => 0.7,
                    p2 => 0.8,
                    p3 => 0.9,
                    p4 => 1.0,
                    μ => 2,
                    η => 3,
                    ξ => 4,
                    α => 5,
                    ρ => 6,
                    dx => 7,
                    dy => 8,
                ]),
            ) - [
                (246412659*9707653^(1/2))/54362856800 + 188137/2800,
                (570735993*9707653^(1/2))/108725713600 + 434999/5600,
                - (344355861*9707653^(1/2))/54362856800 - 63379/560,
                - (797590047*9707653^(1/2))/108725713600 - 29081/224,
                (443853717*9707653^(1/2))/54362856800 + 127423/560,
                (1028044959*9707653^(1/2))/108725713600 + 58897/224,
                - (69182103*9707653^(1/2))/10872571360 - 378269/2800,
                - (160238181*9707653^(1/2))/21745142720 - 881023/5600,
             ]) ≈ 0 atol=1e-12


    end

    @testset "Rp" begin
        @test norm(
            SymbolicUtils.substitute(
                Rp,
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
            ) - [0, -(sqrt(41)) / 164, (sqrt(41)) / 82, -(sqrt(41)) / 164],
        ) ≈ 0 atol = 1e-12
    end

end

Ru = doubleIntegrate(Ru, vars)
Rp = doubleIntegrate(Rp, vars)


# TODO -- test against Integrals.jl and see MATLAB as well; substitute all but integration variable
# TODO -- don't repeat all the same numbers
@testset "Integrated Ru and Rp" begin
    
    @testset "Test integrals" begin

    end

    @testset "Ru" begin
        @test norm(
            SymbolicUtils.substitute(
                Ru,
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
            ],
        ) ≈ 0 atol = 1e-12

        @test norm(
            SymbolicUtils.substitute(
                Ru,
                Dict([
                    u1 => 0,
                    u2 => 0,
                    u3 => 1,
                    u4 => 1,
                    u5 => 1,
                    u6 => 1,
                    u7 => 0,
                    u8 => 0,
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
                -(38^(1/2)*(300*38^(1/2) + 180))/13680,
                (38^(1/2)*(60*38^(1/2) - 180))/13680,
                (38^(1/2)*(120*38^(1/2) - 60))/13680,
                (38^(1/2)*(480*38^(1/2) - 60))/13680,
                (38^(1/2)*(480*38^(1/2) + 300))/13680,
                (38^(1/2)*(120*38^(1/2) + 300))/13680,
                (38^(1/2)*(60*38^(1/2) - 60))/13680,
                -(38^(1/2)*(300*38^(1/2) + 60))/13680
            ],
        ) ≈ 0 atol = 1e-12

        @test norm(
            SymbolicUtils.substitute(
                Ru,
                Dict([
                    u1 => 2,
                    u2 => 2,
                    u3 => 1,
                    u4 => 1,
                    u5 => 1,
                    u6 => 1,
                    u7 => 2,
                    u8 => 2,
                    p1 => 1,
                    p2 => 1,
                    p3 => 1,
                    p4 => 0,
                    μ => 0,
                    η => 0,
                    ξ => 0,
                    α => 1,
                    ρ => 1,
                    dx => 1,
                    dy => 1,
                ]),
            ) - [
                    (10^(1/2)*(150*10^(1/2) - 240))/3600,
                    (10^(1/2)*(120*10^(1/2) + 360))/3600,
                    -(10^(1/2)*(150*10^(1/2) + 30))/3600,
                    (10^(1/2)*(150*10^(1/2) - 60))/3600,
                    -(10^(1/2)*(120*10^(1/2) - 300))/3600,
                    -(10^(1/2)*(150*10^(1/2) + 240))/3600,
                    (10^(1/2)*(120*10^(1/2) - 30))/3600,
                    -(10^(1/2)*(120*10^(1/2) + 60))/3600,
            ],
        ) ≈ 0 atol = 1e-12

        @test norm(
            SymbolicUtils.substitute(
                Ru,
                Dict([
                    u1 => 2,
                    u2 => 2,
                    u3 => 0,
                    u4 => 0,
                    u5 => 0,
                    u6 => 0,
                    u7 => 2,
                    u8 => 2,
                    p1 => 0,
                    p2 => 0,
                    p3 => 0,
                    p4 => 0,
                    μ => 0,
                    η => 0,
                    ξ => 0,
                    α => 1,
                    ρ => 1,
                    dx => 1,
                    dy => 1,
                ]),
            ) - [
                -(5^(1/2)*(120*5^(1/2) - 600))/1800,
                -(5^(1/2)*(120*5^(1/2) - 600))/1800,
                -(5^(1/2)*(60*5^(1/2) + 120))/1800,
                -(5^(1/2)*(60*5^(1/2) + 120))/1800,
                -(5^(1/2)*(60*5^(1/2) + 360))/1800,
                -(5^(1/2)*(60*5^(1/2) + 360))/1800,
                -(5^(1/2)*(120*5^(1/2) + 120))/1800,
                -(5^(1/2)*(120*5^(1/2) + 120))/1800
            ],
        ) ≈ 0 atol = 1e-12

        @test norm(
            SymbolicUtils.substitute(
                Ru,
                Dict([
                    u1 => 1,
                    u2 => 1,
                    u3 => 1,
                    u4 => 1,
                    u5 => 1,
                    u6 => 1,
                    u7 => 1,
                    u8 => 1,
                    p1 => 0,
                    p2 => 0,
                    p3 => 0,
                    p4 => 0,
                    μ => 0,
                    η => 0,
                    ξ => 0,
                    α => 1,
                    ρ => 1,
                    dx => 1,
                    dy => 1,
                ]),
            ) - [
                (5^(1/2)*(90*5^(1/2) - 360))/1800,
                (5^(1/2)*(90*5^(1/2) - 360))/1800,
                1/4,
                1/4,
                (5^(1/2)*(90*5^(1/2) + 360))/1800,
                (5^(1/2)*(90*5^(1/2) + 360))/1800,
                1/4,
                1/4
            ],
        ) ≈ 0 atol = 1e-12

        @test norm(
            SymbolicUtils.substitute(
                Ru,
                Dict([
                    u1 => 0,
                    u2 => 0,
                    u3 => 0,
                    u4 => 0,
                    u5 => 0,
                    u6 => 0,
                    u7 => 0,
                    u8 => 0,
                    p1 => 1,
                    p2 => 1,
                    p3 => 1,
                    p4 => 1,
                    μ => 0,
                    η => 0,
                    ξ => 0,
                    α => 1,
                    ρ => 1,
                    dx => 1,
                    dy => 1,
                ]),
            ) - [
                1/2,
                1/2,
               -1/2,
                1/2,
               -1/2,
               -1/2,
                1/2,
               -1/2
            ],
        ) ≈ 0 atol = 1e-12
    end

    @testset "Rp" begin
        @test norm(
            SymbolicUtils.substitute(
                Rp,
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
            ) - [-sqrt(41) / 41, 0, sqrt(41) / 41, 0],
        ) ≈ 0 atol = 1e-12

        @test norm(SymbolicUtils.substitute(
            Rp,
            Dict([
                u1 => 0,
                u2 => 0,
                u3 => 0,
                u4 => 0,
                u5 => 0,
                u6 => 0,
                u7 => 0,
                u8 => 0,
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
        )) == 0.0

        @test norm(SymbolicUtils.substitute(
            Rp,
            Dict([
                u1 => 0,
                u2 => 0,
                u3 => 1,
                u4 => 1,
                u5 => 1,
                u6 => 1,
                u7 => 0,
                u8 => 0,
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
            (sqrt(38)*(18*sqrt(38) - 30))/2736 - (5*sqrt(38))/456,
            (sqrt(38)*(18*38^(1/2) - 6))/2736 - sqrt(38)/456,
            (7*sqrt(38))/456 + (sqrt(38)*(18*sqrt(38) + 42))/2736,
            (sqrt(38)*(18*sqrt(38) - 6))/2736 - sqrt(38)/456,
        ]) == 0.0
    end
end

Re = [Ru; Rp]
Je = formJe(Re, s)

@testset "Testing a few Je components" begin
    @test size(Re) == (12,)
    @test size(Je) == (12,12)


    # TEST Re

    @test norm(SymbolicUtils.substitute(
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
        ])
    ) - [
        ((41)^(1 / 2) * (270 * (41)^(1 / 2) - 360)) / 14760,
        ((41)^(1 / 2) * (270 * (41)^(1 / 2) - 360)) / 14760,
        -1 / 4,
        3 / 4,
        -((41)^(1 / 2) * (90 * (41)^(1 / 2) - 360)) / 14760,
        -((41)^(1 / 2) * (90 * (41)^(1 / 2) - 360)) / 14760,
        3 / 4,
        -1 / 4,
        -((41)^(1/2)) / 41,
        0,
        -((41)^(1/2)) / 41,
        0,
    ]) ≈ 0 atol=1e-12


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
    @test (sje[4,1] - (- (41^(1/2)*((3555*41^(1/2))/41 - 60))/14760 - 3/328)) ≈ 0 atol=1e-12
    @test (sje[12,3] - ((41^(1/2)*(6*41^(1/2) + 6))/2952 - 41^(1/2)/492)) ≈ 0 atol=1e-12
    @test (sje[5,9] - (-(41^(1/2)*(30*41^(1/2) + 150))/14760)) ≈ 0 atol=1e-12

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
    ) - (38 / 9) ≈ 0 atol=1e-12

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
    ) - 1 ≈ 0 atol=1e-12

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
