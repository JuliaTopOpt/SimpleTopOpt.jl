module Solutions

export Top88Solution, TophSolution, TopflowSolution

##################################################################################################
# Top88 solution
##################################################################################################
struct Top88Solution
    design::Matrix{Float64}
    converged::Bool

    loop::Int32
    change_hist::Vector{Float64}
    obj_hist::Vector{Float64}
    # xPhys_hist::Any

    function Top88Solution(
        design::Matrix{Float64},
        converged::Bool = False,
        loop::Int64 = 0,
        change_hist::Vector{Float64} = Vector{Float64}(),
        obj_hist::Vector{Float64} = Vector{Float64}()
    )
        new(design, converged, loop, change_hist, obj_hist)
    end
end

##################################################################################################
# TopH solution
##################################################################################################

struct TophSolution
    design::Matrix{Float64}
    converged::Bool

    loop::Int32
    change_hist::Vector{Float64}
    obj_hist::Vector{Float64}
    # xPhys_hist::Any

    function TophSolution(
        design::Matrix{Float64},
        converged::Bool=False,
        loop::Int64 = 0,
        change_hist::Vector{Float64} = Vector{Float64}(),
        obj_hist::Vector{Float64} = Vector{Float64}()
    )
        new(design, converged, loop, change_hist, obj_hist)
    end
end

##################################################################################################
# Topflow solution
##################################################################################################

struct TopflowSolution
    design::Matrix{Float64}
    converged::Bool

    loop::Int64
    change_hist::Vector{Float64}
    obj_hist::Vector{Float64}
    # xPhys_hist::Any

    function TopflowSolution(
        design::Matrix{Float64},
        converged = False,
        loop = 0,
        change_hist = Vector{Float64}(),
        obj_hist = Vector{Float64}(),
        # xPhys_hist,
    )
        new(design, converged, loop, change_hist, obj_hist)
    end
end



end