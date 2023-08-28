abstract type FiniteElementDefinition end

"""
Maintains the Topflow finite element formulation
"""
struct TopflowFEA <: FiniteElementDefinition
    neltot::Int64
    doftot::Int64

    nodx::Int64
    nody::Int64
    nodtot::Int64

    edofMat::Matrix{Int64}

    iJ::Matrix{Int64}
    jJ::Matrix{Int64}
    iR::Matrix{Int64}
    jR::Matrix{Int64}
    jE::Matrix{Int64}

    function TopflowFEA(domain::TopflowDomain)
        """
        Assembles most of the finite element matrix problem.
        """

        nodx = domain.nelx + 1
        nody = domain.nely + 1
        nodtot = nodx * nody
        neltot = domain.nelx * domain.nely
        doftot = 3 * nodtot

        nodenrs = reshape(1:nodtot, nody, nodx)
        edofVecU = reshape(2 * nodenrs[1:end-1, 1:end-1] .+ 1, neltot, 1)
        edofMatU =
            repeat(edofVecU, 1, 8) +
            repeat([0 1 (2 * domain.nely .+ [2 3 0 1]) -2 -1], neltot, 1)
        edofVecP = reshape(nodenrs[1:end-1, 1:end-1], neltot, 1)
        edofMatP = repeat(edofVecP, 1, 4) + repeat([1 (domain.nely .+ [2 1]) 0], neltot, 1)

        edofMat = [edofMatU (2 * nodtot .+ edofMatP)]

        iJ = reshape(kron(edofMat, ones(Int64, 12, 1))', 144 * neltot, 1)
        jJ = reshape(kron(edofMat, ones(Int64, 1, 12))', 144 * neltot, 1)
        iR = reshape(edofMat', 12 * neltot, 1)
        jR = ones(12 * neltot, 1)
        jE = repeat((1:neltot)', 12, 1)

        new(neltot, doftot, nodx, nody, nodtot, edofMat, iJ, jJ, iR, jR, jE)
    end
end