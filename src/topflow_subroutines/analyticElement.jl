module analyticElement

using Symbolics
using SymbolicUtils
using LinearAlgebra
using SymbolicNumericIntegration


"""
This script symbolically evaluates and generates code for several
derivatives and functions of interest to TopFlow.
"""

struct Symbols
    """ Contains all symbols involved in the computations """
    ρ::Num
    μ::Num
    α::Num
    dα::Num
    ξ::Num
    η::Num
    dx::Num
    dy::Num
    u1::Num
    u2::Num
    u3::Num
    u4::Num
    u5::Num
    u6::Num
    u7::Num
    u8::Num
    p1::Num
    p2::Num
    p3::Num
    p4::Num

    function Symbols()
        """ Constructor """
        @variables ρ μ α dα ξ η dx dy
        @variables u1 u2 u3 u4 u5 u6 u7 u8 p1 p2 p3 p4

        new(ρ, μ, α, dα, ξ, η, dx, dy, u1, u2, u3, u4, u5, u6, u7, u8, p1, p2, p3, p4)
    end
end


function generation()
    """ Runs through symbolic computations; generates and exports code """
    vars = Symbols()
    DIR_PATH = @__DIR__

    # Analytic part
    # Shape functions and matrices
    xv, yv, Np, Nu = shapeFunctionsAndMatrices(vars)

    # Nodal coordinates, interpolation and coordinate transforms
    x, y = nodalCoordsTransforms(vars, xv, yv, Np)

    # Jacobian
    iJ, detJ, _ = jacobianConstruction(vars, x, y)

    # Derivatives of shape functions
    dNpdx, dNudx = shapeFunctionDerivatives(iJ, Np, Nu, vars)

    # Nodal DOFs
    s, ux, px, dudx, dpdx = nodalDofs(vars, Nu, Np, dNudx, dNpdx)
    # Stabilization parameters τ
    τ, _ = stabilizationParameters(vars, ux)

    # Loop over the tensor weak form to form residual
    Ru, Rp = residualFormation(vars, τ, Nu, ux, dNudx, dudx, dpdx, px, Np, dNpdx, detJ)


    # TODO -- convert to pass by reference
    Ru = doubleIntegrate(Ru, vars)
    Rp = doubleIntegrate(Rp, vars)

    Re = [Ru; Rp]

    Je = formJe(Re, s)

    # Export residual and Jacobian
    JAC = Symbolics.build_function(
        Je,
        vars.dx,
        vars.dy,
        vars.μ,
        vars.ρ,
        vars.α,
        [
            vars.u1,
            vars.u2,
            vars.u3,
            vars.u4,
            vars.u5,
            vars.u6,
            vars.u7,
            vars.u8,
            vars.p1,
            vars.p2,
            vars.p3,
            vars.p4,
        ],
    )
    RES = Symbolics.build_function(
        Re,
        vars.dx,
        vars.dy,
        vars.μ,
        vars.ρ,
        vars.α,
        [
            vars.u1,
            vars.u2,
            vars.u3,
            vars.u4,
            vars.u5,
            vars.u6,
            vars.u7,
            vars.u8,
            vars.p1,
            vars.p2,
            vars.p3,
            vars.p4,
        ],
    )

    # Optimization part
    # Compute ϕ
    ϕ = computePhi(vars, ux, dudx)
    ϕ = doubleIntegrate(ϕ, vars)

    # Compute partial derivative wrt. design field
    dphidg = computePartialPhiDF(vars, ϕ)
    # Compute partial derivative wrt. state field
    dphids = computePartialPhiSF(ϕ, s)

    # Compute partial derivative of residual wrt. design field
    drdg = computePartialJeDF(vars, Re)

    PHI = Symbolics.build_function(
        ϕ,
        vars.dx,
        vars.dy,
        vars.μ,
        vars.α,
        [
            vars.u1,
            vars.u2,
            vars.u3,
            vars.u4,
            vars.u5,
            vars.u6,
            vars.u7,
            vars.u8,
            vars.p1,
            vars.p2,
            vars.p3,
            vars.p4,
        ],
    )
    dPHIdg = Symbolics.build_function(
        dphidg,
        vars.dx,
        vars.dy,
        vars.μ,
        vars.α,
        vars.dα,
        [
            vars.u1,
            vars.u2,
            vars.u3,
            vars.u4,
            vars.u5,
            vars.u6,
            vars.u7,
            vars.u8,
            vars.p1,
            vars.p2,
            vars.p3,
            vars.p4,
        ],
    )
    dPHIds = Symbolics.build_function(
        dphids[1:8],
        vars.dx,
        vars.dy,
        vars.μ,
        vars.α,
        [
            vars.u1,
            vars.u2,
            vars.u3,
            vars.u4,
            vars.u5,
            vars.u6,
            vars.u7,
            vars.u8,
            vars.p1,
            vars.p2,
            vars.p3,
            vars.p4,
        ],
    )

    dRESdg = Symbolics.build_function(
        drdg,
        vars.dx,
        vars.dy,
        vars.μ,
        vars.ρ,
        vars.α,
        vars.dα,
        [
            vars.u1,
            vars.u2,
            vars.u3,
            vars.u4,
            vars.u5,
            vars.u6,
            vars.u7,
            vars.u8,
            vars.p1,
            vars.p2,
            vars.p3,
            vars.p4,
        ],
    )

    return JAC, RES, PHI, dPHIdg, dPHIds, dRESdg
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

function nodalCoordsTransforms(vars, xv, yv, Np)
    """ Nodal coordinates transform """
    xc = vars.dx / 2 * xv
    yc = vars.dy / 2 * yv
    x = Np' * xc
    y = Np' * yc

    return x, y
end

function jacobianConstruction(vars::Symbols, x, y)
    """ Compute Jacobian wrt ξ, η """
    J = Symbolics.jacobian(vec([x y]), [vars.ξ; vars.η])
    iJ = inv(J)
    detJ = simplify(det(J))

    return iJ, detJ, J
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

    return s, ux, px, dudx, dpdx
end

function stabilizationParameters(vars, ux)
    """ Stabilisation parameter τ calculation """
    h = sqrt(vars.dx^2 + vars.dy^2)
    u0 = SymbolicUtils.substitute(ux, Dict([vars.ξ => 0, vars.η => 0]))
    ue = sqrt(u0' * u0)
    τ1 = h / (2 * ue)
    τ3 = vars.ρ * h^2 / (12 * vars.μ)
    τ4 = vars.ρ / vars.α
    τ = (τ1^(-2) + τ3^(-2) + τ4^(-2))^(-1 / 2)

    return τ, ue
end

function residualFormation(vars, τ, Nu, ux, dNudx, dudx, dpdx, px, Np, dNpdx, detJ)
    """ Forms Ru and Rp/ loop over tensor weak form """
    Ru = zeros(Num, 8, 1)
    Rp = zeros(Num, 4)

    # Momentum equations
    for g = 1:8
        for i = 1:2
            Ru[g] += vars.α * Nu[i, g] * ux[i]                                        # Brinkman term
            for j = 1:2
                Ru[g] += vars.μ * dNudx[i, g, j] * (dudx[i, j] + dudx[j, i])           # Viscous term
                Ru[g] += vars.ρ * Nu[i, g] * ux[j] * dudx[i, j]                        # Convection term
                Ru[g] += τ * ux[j] * dNudx[i, g, j] * vars.α * ux[i]                   # SUPG Brinkman term
                for k = 1:2
                    Ru[g] += τ * ux[j] * dNudx[i, g, j] * vars.ρ * ux[k] * dudx[i, k]   # SUPG convection term
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
            Rp[g] += (τ / vars.ρ) * dNpdx[i, g] * vars.α * ux[i]                             # PSPG Brinkman term
            for j = 1:2
                Rp[g] += τ * dNpdx[i, g] * ux[j] * dudx[i, j]                     # PSPG convection term
            end
            Rp[g] += (τ / vars.ρ) * dNpdx[i, g] * dpdx[i]                               # PSPG pressure term
        end
    end

    Ru = vec(simplify(detJ * Ru))
    Rp = simplify(detJ * Rp)

    return Ru, Rp
end

# TODO -- this function is currently broken
#   the intent is to have G be populated with the 
#   element-wise results f[i]
function doubleIntegrate(expression, vars)
    """ Integrates over the unit square in ξ and η """
    F = integrate.(expression, vars.ξ; symbolic = true)
    # G = {}
    for i in eachindex(F)
        f = F[i]
        @assert f[2] == 0
        @assert f[3] == 0
        F[i] = f[1]
    end
    F =
        simplify.(
            SymbolicUtils.substitute(F[1], Dict([vars.ξ => 1])) -
            SymbolicUtils.substitute(F[1], Dict([vars.ξ => -1])),
        )

    F = integrate.(F, vars.η; symbolic = true)
    for i in eachindex(F)
        f = F[i]
        @assert f[2] == 0
        @assert f[3] == 0
        F[i] = f[1]
    end
    F =
        simplify.(
            SymbolicUtils.substitute(F[1], Dict([vars.η => 1])) -
            SymbolicUtils.substitute(F[1], Dict([vars.η => -1])),
        )

    return F
end

function formJe(Re, s)
    """ Form Jacobian """
    #Je = zeros(Num, 12, 12)
    #for b = 1:12
    #for k = 1:12
    #Je[k, b] = Symbolics.derivative(Re[k], s[b]) 
    #end
    #end

    Je = Symbolics.jacobian(Re, s)

    return Je
end


function computePhi(vars::Symbols, ux, dudx)
    """ Compute objective functional """
    ϕ = (1 / 2) * vars.α * ux' * ux
    for i = 1:2
        for j = 1:2
            ϕ += (1 / 2) * vars.μ * dudx[i, j] * (dudx[i, j] + dudx[j, i])
        end
    end

    return ϕ
end

function computePartialPhiDF(vars, ϕ)
    """ Computes the partial of the objective function wrt. design field """

    return simplify(Symbolics.derivative(ϕ, vars.α) * vars.dα)
end

function computePartialPhiSF(ϕ, s)
    """ Computes the partial of the objective function wrt. state field """
    dphids = zeros(Num, 12)
    for a = 1:12
        dphids[a] = simplify(Symbolics.derivative(ϕ, s[a]))
    end

    return dphids
end

function computePartialJeDF(vars, Re)
    """ Computes the partial of the residual wrt. the design field """
    return simplify(Symbolics.derivative(Re, vars.α) * vars.dα)
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    analyticElement.generation()
end
