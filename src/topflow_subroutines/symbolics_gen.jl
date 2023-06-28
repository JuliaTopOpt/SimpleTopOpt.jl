using Symbolics
using SymbolicUtils
using LinearAlgebra

"""
TODO -- some script description
"""

const buildJR = true
const buildPhi = true
const exportJR = true
const exportPhi = true

# Declare variables
@variables ρ μ α dα ξ η dx dy
@variables u1 u2 u3 u4 u5 u6 u7 u8 p1 p2 p3 p4

∂η = Differential(η)
∂ξ = Differential(ξ)


# Shape functions and matrices
xv = [-1, 1, 1, -1]
yv = [-1, -1, 1, 1]
Np = 1 / 4 * (1 .+ ξ * xv) .* (1 .+ η * yv)
Nu = zeros(Num, 2, 8)
Nu[1, 1:2:end-1] = Np
Nu[2, 2:2:end] = Np

# Nodal coordinates, interpolation and coordinate transforms
xc = dx / 2 * xv;
yc = dy / 2 * yv;
x = Np' * xc;
y = Np' * yc;

J = Symbolics.jacobian(vec([x y]), [ξ; η])
iJ = inv(J)
detJ = det(J)

# Derivatives of shape functions
dNpdx =
    iJ[:, 1] * (Symbolics.derivative(Np, ξ))' + iJ[:, 2] * (Symbolics.derivative(Np, η))'
dNudx = zeros(Num, 2, 8, 2)
for i = 1:2
    dNudx[i, :, :] =
        (
            iJ[:, 1] * Symbolics.derivative(Nu[1, :], ξ)' +
            iJ[:, 2] * Symbolics.derivative(Nu[1, :], η)'
        )'
end


# Nodal DOFs
u = [u1; u2; u3; u4; u5; u6; u7; u8]
p = [p1; p2; p3; p4]
s = [u; p]

ux = Nu * u # CORRECT
px = Np' * p # CORRECT

dudx = [dNudx[:, :, 1] * u dNudx[:, :, 2] * u] # NOTE: THESE PRODUCTS DO NOT WORK OUT AS INTENDED; THE EVEN INDICES ARE MIA
dpdx = dNpdx * p # THIS IS CORRECT

# Stabilisation parameters
h = sqrt(dx^2 + dy^2)
u0 = SymbolicUtils.substitute(ux, Dict([ξ => 0, η => 0]))
ue = sqrt(u0' * u0)
τ1 = h / (2 * ue)
τ3 = ρ * h^2 / (12 * μ)
τ4 = ρ / α
τ = (τ1^(-2) + τ3^(-2) + τ4^(-2))^(-1 / 2)

# Loop over tensor weak form to form residual
if buildJR
    Ru = zeros(Num, 8, 1)
    Rp = zeros(Num, 4)

    # Momentum equations
    for g = 1:8
        for i = 1:2
            Ru[g] += α * Nu[i, g] * ux[i]                                           # Brinkman term
            for j = 1:2
                Ru[g] += μ * dNudx[i, g, j] * (dudx[i, j] + dudx[j, i])            # Viscous term
                Ru[g] += ρ * Nu[i, g] * ux[j] * dudx[i, j]                       # Convection term
                Ru[g] += τ * ux[j] * dNudx[i, g, j] * α * ux[i]                  # SUPG Brinkman term
                for k = 1:2
                    Ru[g] += τ * ux[j] * dNudx[i, g, j] * ρ * ux[k] * dudx[i, k]  # SUPG convection term
                end
                Ru[g] += τ * ux[j] * dNudx[i, g, j] * dpdx[i]                    # SUPG pressure term
            end
            Ru[g] -= dNudx[i, g, i] * px                                         # Pressure term
        end
    end

    # Incompressibility equations
    for g = 1:4
        for i = 1:2
            Rp[g] += Np[1, g] * dudx[i, i]                     # Divergence term
            Rp[g] += τ / ρ * dNpdx[i, g] * α * ux[i]            # PSPG Brinkman term
            for j = 1:2
                Rp[g] += τ * dNpdx[i, g] * ux[j] * dudx[i, j]  # PSPG convection term
            end
            Rp[g] += τ / ρ * dNpdx[i, g] * dpdx[i]              # PSPG pressure term
        end
    end

end

# Export residual and Jacobian


# Optimization part 


# Compute partial derivative of residual wrt. design field


# Export optimization functions
