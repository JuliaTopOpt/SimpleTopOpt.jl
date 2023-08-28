module Solutions

export Top88Solution, TophSolution, TopflowSolution

##################################################################################################
# Top88 solution
##################################################################################################
struct Top88Solution
    xPhys::Matrix{Float64}
    converged::Bool = False

    change_hist::Any
    obj_hist::Any
    xPhys_hist::Any
end

##################################################################################################
# TopH solution
##################################################################################################

struct TophSolution
    xPhys::Matrix{Float64}
    converged::Bool = False

    change_hist::Any
    obj_hist::Any
    xPhys_hist::Any
end

##################################################################################################
# Topflow solution
##################################################################################################

struct TopflowSolution
    # TODO -- typehint the below
    xPhys::Any
    loop::Int64
    change_hist::Any
    obj_hist::Any
    xPhys_hist::Any

    converged::Bool

    function TopflowSolution(
        xPhys,
        loop,
        change_hist,
        obj_hist,
        xPhys_hist,
        converged=True,
    ) where {U<:TopflowContainer}

        new{U}(problem_container, xPhys, loop, change_hist, obj_hist, xPhys_hist)
    end
end



end