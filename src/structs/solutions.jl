module Solutions

export Top88Solution, TophSolution, TopflowSolution

##################################################################################################
# Top88 solution
##################################################################################################
struct Top88Solution
    xPhys::Matrix{Float64}
    converged::Bool = False

    loop::Int32
    change_hist::Vector{Float64}
    obj_hist::Vector{Float64}
    # xPhys_hist::Any

    function Top88Solution(
        xPhys::Matrix{Float64},
        converged::Bool = False,
        loop::Int64 = 0,
        change_hist::Vector{Float64} = Vector{Float64}(),
        obj_hist::Vector{Float64} = Vector{Float64}()
    )
        new(xPhys, converged, loop, change_hist, obj_hist)
    end
end

##################################################################################################
# TopH solution
##################################################################################################

struct TophSolution
    xPhys::Matrix{Float64}
    converged::Bool = False

    loop::Int32
    change_hist::Vector{Float64}
    obj_hist::Vector{Float64}
    # xPhys_hist::Any

    function TophSolution(
        xPhys::Matrix{Float64},
        converged::Bool=False,
        loop::Int64 = 0,
        change_hist::Vector{Float64} = Vector{Float64}(),
        obj_hist::Vector{Float64} = Vector{Float64}()
    )
        new(xPhys, converged, loop, change_hist, obj_hist)
    end
end

##################################################################################################
# Topflow solution
##################################################################################################

struct TopflowSolution
    xPhys::Matrix{Float64}
    converged::Bool = False

    loop::Int64 = 0
    change_hist::Vector{Float64}
    obj_hist::Vector{Float64}
    # xPhys_hist::Any

    function TopflowSolution(
        xPhys,
        converged = False,
        loop = 0,
        change_hist = Vector{Float64}(),
        obj_hist = Vector{Float64}(),
        # xPhys_hist,
    )

        new(xPhys, converged, loop, change_hist, obj_hist, xPhys_hist)
    end
end



end