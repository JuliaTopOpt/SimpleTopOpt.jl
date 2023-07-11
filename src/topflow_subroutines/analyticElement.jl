using Symbolics
using SymbolicUtils
using LinearAlgebra
using SymbolicNumericIntegration

"""
This script symbolically evaluates and generates code for several
derivatives and functions of interest to TopFlow.
"""

function analyticElement(exportJR::Bool = true, exportPhi::Bool = true; test::Bool = false)

    # Declare variables
    @variables ρ μ α dα ξ η dx dy
    @variables u1 u2 u3 u4 u5 u6 u7 u8 p1 p2 p3 p4

    # Shape functions and matrices
    xv = [-1, 1, 1, -1]
    yv = [-1, -1, 1, 1]
    Np = 1 / 4 * (1 .+ ξ * xv) .* (1 .+ η * yv)
    Nu = zeros(Num, 2, 8)
    Nu[1, 1:2:end-1] = Np
    Nu[2, 2:2:end] = Np

    # Nodal coordinates, interpolation and coordinate transforms
    xc = dx / 2 * xv
    yc = dy / 2 * yv
    x = Np' * xc
    y = Np' * yc

    J = Symbolics.jacobian(vec([x y]), [ξ; η])
    iJ = inv(J)
    detJ = simplify(det(J))

    """
    Test Jacobian
    """
    if !test  # TODO -- take off the not!
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



    # Derivatives of shape functions
    dNpdx =
        iJ[:, 1] * (Symbolics.derivative(Np, ξ))' +
        iJ[:, 2] * (Symbolics.derivative(Np, η))'
    dNudx = zeros(Num, 2, 8, 2)
    for i = 1:2
        dNudx[i, :, :] =
            (
                iJ[:, 1] * Symbolics.derivative(Nu[i, :], ξ)' +
                iJ[:, 2] * Symbolics.derivative(Nu[i, :], η)'
            )'
    end

    """
    Test dNpdx and dNudx
    """
    if !test  # TODO -- take off the not!
        println("Testing dNpdx and dNudx")
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


    # Nodal DOFs
    u = [u1; u2; u3; u4; u5; u6; u7; u8]
    p = [p1; p2; p3; p4]
    s = [u; p]

    ux = Nu * u
    px = Np' * p

    dudx = [(dNudx[:, :, 1] * u) (dNudx[:, :, 2] * u)]
    dpdx = dNpdx * p

    """
    Test ux, px, dudx, and dpdx
    """
    if test  # TODO -- take off the not!
        println("Testing ux, px, dudx, and dpdx")
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

    # Stabilisation parameters
    h = sqrt(dx^2 + dy^2)
    u0 = SymbolicUtils.substitute(ux, Dict([ξ => 0, η => 0]))
    ue = sqrt(u0' * u0)
    τ1 = h / (2 * ue)
    τ3 = ρ * h^2 / (12 * μ)
    τ4 = ρ / α
    τ = (τ1^(-2) + τ3^(-2) + τ4^(-2))^(-1 / 2)

    if test  # TODO -- take off the not!
        println("Testing stabilisation parameters")

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

    Ru = zeros(Num, 8, 1)
    Rp = zeros(Num, 4)

    println("Building Ru and Rp...")
    # Momentum equations
    for g = 1:8
        for i = 1:2
            Ru[g] += α * Nu[i, g] * ux[i]                                        # Brinkman term
            for j = 1:2
                Ru[g] += μ * dNudx[i, g, j] * (dudx[i, j] + dudx[j, i])           # Viscous term
                Ru[g] += ρ * Nu[i, g] * ux[j] * dudx[i, j]                        # Convection term
                Ru[g] += τ * ux[j] * dNudx[i, g, j] * α * ux[i]                   # SUPG Brinkman term
                for k = 1:2
                    Ru[g] += τ * ux[j] * dNudx[i, g, j] * ρ * ux[k] * dudx[i, k]   # SUPG convection term
                end
                Ru[g] += τ * ux[j] * dNudx[i, g, j] * dpdx[i]                     # SUPG pressure term
            end
            Ru[g] -= dNudx[i, g, i] * px                                          # Pressure term
        end
    end

    # Incompressibility equations
    for g = 1:4
        for i = 1:2
            Rp[g] += Np[g] * dudx[i, i]                                          # Divergence term
            Rp[g] += (τ / ρ) * dNpdx[i, g] * α * ux[i]                             # PSPG Brinkman term
            for j = 1:2
                Rp[g] += τ * dNpdx[i, g] * ux[j] * dudx[i, j]                     # PSPG convection term
            end
            Rp[g] += (τ / ρ) * dNpdx[i, g] * dpdx[i]                               # PSPG pressure term
        end
    end

    """
    Ensuring intermediate correctness
    """
    if test
        println("Testing intermediate Ru and Rp")

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
            ) - [
                0,
                -(√(41)) / 164,
                 (√(41)) / 82,
                -(√(41)) / 164,
            ],
        ) ≈ 0 atol = 1e-12

    end


    println("Simplifying detJ*Ru and detJ*Rp...")
    Ru = vec(simplify(detJ * Ru))
    Rp = simplify(detJ * Rp)

    if test
        println("Testing simplified Ru and Rp")
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
                (8 * 29^(1/2))/29 + 1/2,
                (8 * 29^(1/2))/29 + 1/2,
                -(16 * 29^(1/2))/29 - 5/2,
                -(16 * 29^(1/2))/29 - 3/2,
                (24 * 29^(1/2))/29 + 15/2,
                (24 * 29^(1/2))/29 + 15/2,
                -(16 * 29^(1/2))/29 - 3/2,
                -(16 * 29^(1/2))/29 - 5/2,
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
                -(√(41)) / 164,
                1/4 - (√41) / 164,
                (√(41)) / 82,
                (√(41)) / 82,
                1/4 - (√41) / 164,
                -(√(41)) / 164,
            ],
        ) ≈ 0 atol = 1e-12 # Passes

    end

    # FToC
    println("Integrating Ru...")
    F = integrate(Ru, ξ; symbolic = true)
    @assert F[2] == 0.0
    @assert F[3] == 0.0
    F = F[1]

    F = SymbolicUtils.substitute(F, Dict([ξ => 1])) - SymbolicUtils.substitute(F, Dict([ξ => -1]))
    F = integrate(F, η; symbolic = true)
    @assert F[2] == 0.0
    @assert F[3] == 0.0
    F = F[1]

    Ru = SymbolicUtils.substitute(F, Dict([η => 1])) - SymbolicUtils.substitute(F, Dict([η => -1]))


    if test
        println("Testing integrated Ru...")
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
                ((41)^(1/2) * (270*(41)^(1/2) - 360)) / 14760,
                ((41)^(1/2) * (270*(41)^(1/2) - 360)) / 14760,
                -1/4,
                3/4,
                -((41)^(1/2) * (90*(41)^(1/2) - 360)) / 14760,
                ((41)^(1/2) * (90*(41)^(1/2) - 360)) / 14760,
                3/4,
                -1/4
            ],
        ) ≈ 0 atol = 1e-12 
    end

    println("Integrating Rp...")
    G = integrate(Rp, ξ; symbolic = true)
    @assert G[2] == 0.0
    @assert G[3] == 0.0
    G = G[1]

    G = simplify(
        SymbolicUtils.substitute(G, Dict([ξ => 1])) -
        SymbolicUtils.substitute(G, Dict([ξ => -1])),
    )
    G = integrate(G, η; symbolic = true)
    @assert G[2] == 0.0
    @assert G[3] == 0.0
    G = G[1]

    Rp = simplify(
        SymbolicUtils.substitute(G, Dict([η => 1])) -
        SymbolicUtils.substitute(G, Dict([η => -1])),
    )

    if test
        println("Testing integrated Rp")
    end

    Re = [Ru; Rp]

    # Differentiate
    Je = Symbolics.jacobian(Re, s)

    # Export residual and Jacobian
    if exportJR
        JAC = Symbolics.build_function(
            Je,
            dx,
            dy,
            μ,
            ρ,
            α,
            [u1, u2, u3, u4, u5, u6, u7, u8, p1, p2, p3, p4],
        )
        RES = Symbolics.build_function(
            Re,
            dx,
            dy,
            μ,
            ρ,
            α,
            [u1, u2, u3, u4, u5, u6, u7, u8, p1, p2, p3, p4],
        )

        write("subroutines/JAC.jl", string(JAC))
        write("subroutines/RES.jl", string(RES))
    end

    # Optimization part 
    println("Computing Φ...")

    ϕ = (1 / 2) * α * ux' * ux

    # Compute partial derivative of residual wrt. design field
    
    # Export optimization functions

end









if abspath(PROGRAM_FILE) == @__FILE__
    analyticElement(true, true)
end
