using Test
using SimpleTopOpt
using SimpleTopOpt.TopFlow


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

@testset "Test ux, px, dudx, and dpdx" begin

    @test size(s) == (12,)
    @test size(p) == (4,)


    println("{TEST -- Testing ux, px, dudx, and dpdx}")
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
            Dict([
                p1 => 1,
                p2 => 1,
                p3 => 1,
                p4 => 1,
                ξ => 1,
                η => 1,
                dx => 1,
                dy => 1,
            ]),
        )) - [0; 0],
    ) == 0

    @test norm(
        (SymbolicUtils.substitute(
            dpdx,
            Dict([
                p1 => 1,
                p2 => 2,
                p3 => 3,
                p4 => 4,
                ξ => 5,
                η => 6,
                dx => 7,
                dy => 8,
            ]),
        )) - [-6 / 7; -3 / 8],
    ) == 0
end


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


end



@testset "Intermediate Ru and Rp" begin
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
            -(√(41)) / 41,
            1 - (√41) / 41,
            2 * (√(41)) / 41,
            2 * (√(41)) / 41,
            1 - (√41) / 41,
            -(√(41)) / 41,
        ],
    ) ≈ 0 atol = 1e-12

    @test norm(
        SymbolicUtils.substitute(
            Ru,
            Dict([
                u1 => 1,
                u2 => 0,
                u3 => 1,
                u4 => 0,
                u5 => 1,
                u6 => 0,
                u7 => 1,
                u8 => 0,
                p1 => 1,
                p2 => 0,
                p3 => 1,
                p4 => 0,
                μ => 1,
                η => 0,
                ξ => 1,
                α => 1,
                ρ => 1,
                dx => 1,
                dy => 1,
            ]),
        ) - [
            (1 / 4) - (√(39) / 78),
            -(√(39) / 78),
            (1 / 4) + (√(39) / 78),
            (1 / 2) + (√(39) / 78),
            (1 / 4) + (√(39) / 78),
            (√(39) / 78) - (1 / 2),
            (1 / 4) - (√(39) / 78),
            -(√(39) / 78),
        ],
    ) ≈ 0 atol = 1e-12


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
        ) - [0, -(sqrt(41)) / 41, 2 * (sqrt(41)) / 41, -(sqrt(41)) / 41],
    ) ≈ 0 atol = 1e-12
end



@testset "Testing simplified Ru and Rp" begin
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
    ) ≈ 0 atol = 1e-12 # Passes


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
    ) ≈ 0 atol = 1e-12 # Passes

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


@testset "Integrated Ru and Rp" begin
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
end




@testset "Testing a few Je components" begin
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
        sje[1, 4] - (-(sqrt(41) * ((3555 * sqrt(41)) / 41 - 60)) / 14760 - 3 / 328)
    ) ≈ 0 atol = 1e-12
    #@test (sje[4,1] - (-(sqrt(41) * ((3555 * sqrt(41))/ 41 - 60))/14760 - 3/328)) ≈ 0 atol=1e-12
    #@test (sje[12,3] - ((sqrt(41)*(6*sqrt(41) + 6))/2952 - sqrt(41)/492)) ≈ 0 atol=1e-12
    #@test (sje[5,9] - (-(41^(1/2)*(30*41^(1/2) + 150))/14760)) ≈ 0 atol=1e-12




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
    ) - 38 / 9 == 0

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

end


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




