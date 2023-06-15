module TopFlow

using LinearAlgebra
using SparseArrays


export topflow

"""
    `topflow`

Fluidic topology optimization
"""
function topflow(problem_container::T) where {T<:TopflowContainer}
    #TODO: do I need to error check the physical parameters, etc?



end

end
