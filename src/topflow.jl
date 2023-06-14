module TopFlow

using LinearAlgebra
using SparseArrays


export topflow


### Constant parameters
# Physical parameters
const Uin = 1e0
const rho = 1e0
const mu = 1e0

# Continuation strategy
const conit = 50

# Optimization parameters
const mvlim = 0.2
const chlim = 1e-5
const chnum = 5

# Newton Solver parameters
const nltol = 1e-6
const nlmax = 25

"""
    `topflow`

Driving function for the TopFlow module
"""
function topflow(
    problemtype::S = 1,
    volfrac::T = 0.3333,
    Lx::T = 1.0,
    Ly::T = 1.0,
    nely::S = 30,
) where {T<:AbstractFloat,S<:Integer}

end

"""
Struct containing all components required for the finite element analysis
"""
mutable struct FEAContainer{T<:AbstractFloat,S<:Integer}
    dx::T
    dy::T

    nodx::S
    nody::S
    nodtot::S

    neltot::S
    doftot::S

    nodenrs::Matrix{S}
    edofMat::Matrix{S}

    iJ::Matrix{S}
    jJ::Matrix{S}
    iR::Matrix{S}
    jR::Matrix{S}
    jE::Matrix{S}

    """
    Constructor
    """
    function FEAContainer(
        nely::S = 30,
        Lx::T = 1.0,
        Ly::T = 1.0,
    ) where {T<:AbstractFloat,S<:Integer}

        @assert nely > 0.0
        @assert Lx > 0.0
        @assert Ly > 0.0

        nelx = nely * Lx / Ly

        dx = Lx / nelx
        dy = Ly / nely
        nodx = nelx + 1
        nody = nely + 1
        nodtot = nodx * nody
        neltot = nelx * nely
        doftot = 3 * nodtot

        nodenrs = reshape(1:nodtot, nody, nodx)
        edofVecU = reshape(2 * nodenrs[1:end-1, 1:end-1] + 1, neltot, 1)
        edofMatU =
            repeat(edofVecU, 1, 8) + repeat([0 1 2 * nely + [2 3 0 1] -2 -1], neltot, 1)
        edofVecP = reshape(nodenrs[1:end-1, 1:end-1], neltot, 1)
        edofMatP = repeat(edofVecP, 1, 4) + repeat([1 nely + [2 1] 0], neltot, 1)
        edofMat = [edofMatU 2 * nodtot + edofMatP]

        iJ = reshape(kron(edofMat, ones(12, 1))', 144 * neltot, 1)
        jJ = reshape(kron(edofMat, ones(1, 12))', 144 * neltot, 1)
        iR = reshape(edofMat', 12 * neltot, 1)
        jR = ones(12 * neltot, 1)
        jE = repeat(1:neltot, 12, 1)

        new{T,S}(
            dx,
            dy,
            nodx,
            nody,
            nodtot,
            neltot,
            doftot,
            nodenrs,
            edofVecU,
            edofMatU,
            edofVecP,
            edofMatP,
            edofMat,
            iJ,
            jJ,
            iR,
            jR,
            jE,
        )

    end
end

end
