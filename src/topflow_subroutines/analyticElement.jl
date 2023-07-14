using Symbolics
using SymbolicUtils
using LinearAlgebra
using SymbolicNumericIntegration

module analyticElement
"""
This script symbolically evaluates and generates code for several
derivatives and functions of interest to TopFlow.
"""

struct Symbols
    """ Contains all symbols involved in the computations """
    ρ
    μ
    α
    dα
    ξ
    η
    dx
    dy
    u1
    u2
    u3
    u4
    u5
    u6
    u7
    u8
    p1
    p2
    p3
    p4

    function Symbols()
        """ Constructor """
        @variables ρ μ α dα ξ η dx dy
        @variables u1 u2 u3 u4 u5 u6 u7 u8 p1 p2 p3 p4

        new(
            ρ,μ,α,dα,ξ,η,dx,dy,u1,u2,u3,u4,u5,u6,u7,u8,p1,p2,p3,p4
        )
    end
end


function generation(exportJR::Bool = true, exportPhi::Bool = true)
    """ Runs through symbolic computations; generates and exports code """
    vars = symbols()
    DIR_PATH = @__DIR__

    # Analytic part
    # Shape functions and matrices
    xv, yv, Np, Nu = shapeFunctionsAndMatrices(vars)

    # Nodal coordinates, interpolation and coordinate transforms
    x, y = nodalCoordsTransforms(xv, yv, Np)

    # Jacobian
    iJ, detJ = jacobianConstruction(vars, x, y)

    # Derivatives of shape functions
    dNpdx, dNudx = shapeFunctionDerivatives(iJ, Np, Nu, vars)

    # Loop over the tensor weak form to form residual
    Ru, Rp = residualFormation(vars, Nu, ux, dNudx, dudx, dpdx, px, Np, dNpdx, detJ)

    # TODO -- convert to pass by reference
    Ru = doubleIntegrate(Ru)
    Rp = doubleIntegrate(Rp)

    Re = [Ru; Rp]

    Je = formJe(Re, s)

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

        println("{exportJR} -- writing out JAC and RES")

        write(DIR_PATH * "subroutines/JAC.jl", string(JAC))
        write(DIR_PATH * "subroutines/RES.jl", string(RES))
    end

    # Optimization part
    # Compute ϕ
    ϕ = computePhi(vars, ux, dudx)
    ϕ = doubleIntegrate(ϕ)

    # Compute partial derivative wrt. design field
    dphidg = computePartialPhiDF(vars, ϕ)
    # Compute partial derivative wrt. state field
    dphids = computePartialPhiSF(ϕ, s)

    # Compute partial derivative of residual wrt. design field
    drdg = computePartialJeDF(vars, Re)

    if exportPhi
        PHI = Symbolics.build_function(
            ϕ,
            dx,
            dy,
            μ,
            α,
            [u1, u2, u3, u4, u5, u6, u7, u8, p1, p2, p3, p4],
        )
        dPHIdg = Symbolics.build_function(
            dphidg,
            dx,
            dy,
            μ,
            α,
            dα,
            [u1, u2, u3, u4, u5, u6, u7, u8, p1, p2, p3, p4],
        )
        dPHIds = Symbolics.build_function(
            dphids[1:8],
            dx,
            dy,
            μ,
            α,
            [u1, u2, u3, u4, u5, u6, u7, u8, p1, p2, p3, p4],
        )

        write(DIR_PATH * "subroutines/PHI.jl", string(PHI))
        write(DIR_PATH * "subroutines/dPHIdg.jl", string(dPHIdg))
        write(DIR_PATH * "/subroutines/dPHIds.jl", string(dPHIds))
    end
    if exportJR
        dRESdg = Symbolics.build_function(
            drdg,
            dx,
            dy,
            μ,
            ρ,
            α,
            dα,
            [u1, u2, u3, u4, u5, u6, u7, u8, p1, p2, p3, p4],
        )

        write(DIR_PATH * "/subroutines/dRESdg", string(dRESdg))
    end

    return
end



function shapeFunctionsAndMatrices(vars::Symbols)
    """ Assemble shape functions """
    xv = [-1, 1, 1, -1]
    yv = [-1, -1, 1, 1]
    Np = 1 / 4 * (1 .+ vars.ξ * xv) .* (1 .+ vars.η * yv)
    Nu = zeros(Num, 2, 8)
    Nu[1, 1:2:end-1] = Np
    Nu[2, 2:2:end] = Np

    return xv, yv, Np, Nu
end

function nodalCoordsTransforms(xv, yv, Np)
    """ Nodal coordinates transform """
    xc = dx / 2 * xv
    yc = dy / 2 * yv
    x = Np' * xc
    y = Np' * yc

    return x, y
end

function jacobianConstruction(vars::Symbols, x, y)
    """ Compute Jacobian wrt ξ, η """
    J = Symbolics.jacobian(vec([x y]), [vars.ξ; vars.η])
    iJ = inv(J)
    detJ = simplify(det(J))

    return iJ, detJ
end

function shapeFunctionDerivatives(iJ, Np, Nu, vars)
    """ Compute shape function derivatives in η, ξ """
    dNpdx =
        iJ[:, 1] * (Symbolics.derivative(Np, vars.ξ))' +
        iJ[:, 2] * (Symbolics.derivative(Np, vars.η))'
    dNudx = zeros(Num, 2, 8, 2)
    for i = 1:2
        dNudx[i, :, :] =
            (
                iJ[:, 1] * Symbolics.derivative(Nu[i, :], vars.ξ)' +
                iJ[:, 2] * Symbolics.derivative(Nu[i, :], vars.η)'
            )'
    end

    return dNpdx, dNudx
end

function nodalDofs(vars, Nu, Np, dNudx, dNpdx)
    """ Assemble nodal dof """
    u = [vars.u1; vars.u2; vars.u3; vars.u4; vars.u5; vars.u6; vars.u7; vars.u8]
    p = [vars.p1; vars.p2; vars.p3; vars.p4]
    s = [u; p]

    ux = Nu * u
    px = Np' * p

    dudx = [(dNudx[:, :, 1] * u) (dNudx[:, :, 2] * u)]
    dpdx = dNpdx * p

    return u, p , s, ux, px, dudx, dpdx
end

function stabilizationParameters(vars, ux)
    # Stabilisation parameters
    h = sqrt(vars.dx^2 + vars.dy^2)
    u0 = SymbolicUtils.substitute(ux, Dict([vars.ξ => 0, vars.η => 0]))
    ue = sqrt(u0' * u0)
    τ1 = h / (2 * ue)
    τ3 = vars.ρ * h^2 / (12 * vars.μ)
    τ4 = vars.ρ / vars.α
    τ = (τ1^(-2) + τ3^(-2) + τ4^(-2))^(-1 / 2)

    return τ, ue, u0, h
end

function residualFormation(vars, Nu, ux, dNudx, dudx, dpdx, px, Np, dNpdx, detJ)
    Ru = zeros(Num, 8, 1)
    Rp = zeros(Num, 4)

    # Momentum equations
    for g = 1:8
        for i = 1:2
            Ru[g] += vars.α * Nu[i, g] * ux[i]                                        # Brinkman term
            for j = 1:2
                Ru[g] += vars.μ * dNudx[i, g, j] * (dudx[i, j] + dudx[j, i])           # Viscous term
                Ru[g] += vars.ρ * Nu[i, g] * ux[j] * dudx[i, j]                        # Convection term
                Ru[g] += vars.τ * ux[j] * dNudx[i, g, j] * vars.α * ux[i]                   # SUPG Brinkman term
                for k = 1:2
                    Ru[g] += vars.τ * ux[j] * dNudx[i, g, j] * vars.ρ * ux[k] * dudx[i, k]   # SUPG convection term
                end
                Ru[g] += vars.τ * ux[j] * dNudx[i, g, j] * dpdx[i]                     # SUPG pressure term
            end
            Ru[g] -= dNudx[i, g, i] * px                                          # Pressure term
        end
    end

    # Incompressibility equations
    for g = 1:4
        for i = 1:2
            Rp[g] += Np[g] * dudx[i, i]                                          # Divergence term
            Rp[g] += (vars.τ / vars.ρ) * dNpdx[i, g] * vars.α * ux[i]                             # PSPG Brinkman term
            for j = 1:2
                Rp[g] += vars.τ * dNpdx[i, g] * ux[j] * dudx[i, j]                     # PSPG convection term
            end
            Rp[g] += (vars.τ / vars.ρ) * dNpdx[i, g] * dpdx[i]                               # PSPG pressure term
        end
    end

    Ru = vec(simplify(detJ * Ru))
    Rp = simplify(detJ * Rp)

    return Ru, Rp
end

function doubleIntegrate(expression)
    """ Integrates over the unit square in ξ and η """
    F = integrate(expression, ξ; symbolic = true)
    @assert F[2] == 0
    @assert F[3] == 0
    F = simplify(
        SymbolicUtils.substitute(F[1], Dict([ξ => 1])) -
        SymbolicUtils.substitute(F[1], Dict([ξ => -1]))
    )

    G = integrate(expression, η; symbolic = true)
    @assert G[2] == 0
    @assert G[3] == 0
    G = simplify(
        SymbolicUtils.substitute(G[1], Dict([η => 1])) -
        SymbolicUtils.substitute(G[1], Dict([η => -1]))
    )

    return G
end

function formJe(Re, s)
    Je = zeros(Num, 12, 12)
    for b = 1:12
        Je[:,b] = Symbolics.derivative(Re, s[b])
    end

    return Je
end


function computePhi(vars::Symbols, ux, dudx)
    ϕ = (1/2) * vars.α * ux' * ux
    for i = 1:2
        for j = 1:2
            ϕ += (1/2) * vars.μ * dudx[i,j] * (dudx[i,j] + dudx[j,i])
        end
    end

    return ϕ
end

function computePartialPhiDF(vars, ϕ)
    return simplify(Symbolics.derivative(ϕ, vars.α) * vars.dα)
end

function computePartialPhiSF(ϕ, s)
    dphids = zeros(Num, 12)
    for a = 1:12
        dphids[a] = simplify(Symbolics.dervative(ϕ, s[a]))
    end

    return dphids
end

function computePartialJeDF(vars, Re)
    return simplify( Symbolics.derivative(Re, vars.α) * dα)
end

