module TopFlow

using LinearAlgebra
using SparseArrays

import ..TopflowContainer
export topflow


"""
    `topflow`

Fluidic topology optimization
"""
function topflow(problem_container::T) where {T<:TopflowContainer}
    #TODO: do I need to error check the physical parameters, etc?
    tfdc = problem_container.tfdc
    bc = problem_container.bc

    ### Initialization
    # Solution vector
    S = zeroes(pc.fea.doftot, 1)
    dS = copy(S)
    L = copy(S)
    S[pc.bc.fixedDofs] = bc.DIR[bc.fixedDofs]

    # Design field
    xPhys = pc.volfrac * ones(tfdc.nelx, tfdc.nely)

    # Counters
    loop = 0
    loopcont = 0
    nlittot = 0
    chcnt = 0

    change = Inf
    objOld = Inf


end

end
