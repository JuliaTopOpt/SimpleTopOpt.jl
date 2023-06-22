using Symbolics

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

# Shape functions and matrices
xv = [-1 1 1 -1]';
yv = [-1 -1 1 1]';
Np = 1 / 4 * (1 .+ ξ * xv') .* (1 .+ η * yv')
Nu = zeros(Num, 2, 8)
Nu[1, 1:2:end-1] = Np
Nu[2, 2:2:end] = Np

# Nodal coordinates, interpolation and coordinate transforms
xc = dx / 2 * xv;
yc = dy / 2 * yv;
x = Np * xc;
y = Np * yc;

# Declare differentials
∂η = Differential(η)
∂ξ = Differential(ξ)

# ???????????????????
