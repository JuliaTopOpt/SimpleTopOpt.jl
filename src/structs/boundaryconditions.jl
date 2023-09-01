module BoundaryConditionDefinitions

using ..ParameterDefinitions
using ..FEMDefinitions
using ..Domains

export BoundaryConditions, Top88BoundaryConditions, TophBoundaryConditions
export DoublePipeBC, PipeBendBC

abstract type BoundaryConditions end

##################################################################################################
# Top88 BCs
##################################################################################################
abstract type Top88BoundaryConditions <: BoundaryConditions end

# NOTE: only a single set of boundary conditions is implemented currently

##################################################################################################
# Toph BCs
##################################################################################################
abstract type TophBoundaryConditions <: BoundaryConditions end

# NOTE: only a single set of boundary conditions is implemented currently

##################################################################################################
# Topflow BCs
##################################################################################################
abstract type TopflowBoundaryConditions <: BoundaryConditions end

"""@docs
    DoublePipeBC

Specify the boundary conditions for the Topflow DoublePipe problem (type 1)
"""
struct DoublePipeBC <: TopflowBoundaryConditions
    fixedDofs::Matrix{Int64}
    DIR::Matrix{Float64}

    inletLength::Float64

    """
    Constructor
    """
    function DoublePipeBC(domain::TopflowDomain, fea::TopflowFEA, Uin::Float64)
        if mod(domain.nely, 6) > 0
            throw(ArgumentError("Number of elements in y-dir must be divisible by 6."))
        elseif domain.nely < 30
            throw(ArgumentError("Number of elements in y-dir must be at least 30."))
        end

        inletLength = 1 / 6 * domain.nely
        inlet1 = 1 / 6 * domain.nely + 1
        inlet2 = 4 / 6 * domain.nely + 1
        outletLength = 1 / 6 * domain.nely
        outlet1 = 1 / 6 * domain.nely + 1
        outlet2 = 4 / 6 * domain.nely + 1

        nodesInlet = Int.([(inlet1:(inlet1+inletLength))' (inlet2:(inlet2+inletLength))'])
        nodesOutlet =
            (fea.nodx - 1) * fea.nody .+
            Int.([(outlet1:(outlet1+outletLength))' (outlet2:(outlet2+outletLength))'])

        vec_TopBot = Int.(vcat(1:domain.nely+1:fea.nodtot, fea.nody:fea.nody:fea.nodtot))
        vec_LefRig = Int.(vcat(2:domain.nely, (domain.nelx)*fea.nody+2:fea.nodtot-1))

        nodesTopBot = reshape(vec_TopBot, 1, length(vec_TopBot))
        nodesLefRig = reshape(vec_LefRig, 1, length(vec_LefRig))
        fixedDofsTBx = 2 * nodesTopBot .- 1
        fixedDofsTBy = 2 * nodesTopBot
        vec_LRx = 2 * Int.(setdiff(nodesLefRig, [nodesInlet nodesOutlet])) .- 1
        fixedDofsLRx = reshape(vec_LRx, 1, length(vec_LRx))
        fixedDofsLRy = 2 * nodesLefRig
        fixedDofsInX = 2 * nodesInlet .- 1
        fixedDofsInY = 2 * nodesInlet
        fixedDofsOutY = 2 * nodesOutlet
        fixedDofsOutP = 2 * fea.nodtot .+ nodesOutlet


        fixedDofsU = [
            fixedDofsTBx fixedDofsTBy fixedDofsLRx fixedDofsLRy fixedDofsInX fixedDofsInY fixedDofsOutY
        ]
        fixedDofsP = fixedDofsOutP
        fixedDofs = [fixedDofsU fixedDofsP]

        Uinlet = [-4 * x .^ 2 .+ 4 * x for x in (Uin * (0:inletLength) / inletLength)]
        DIRU = zeros(fea.nodtot * 2, 1)
        DIRU[fixedDofsInX] = [Uinlet' Uinlet']
        DIRP = zeros(fea.nodtot, 1)
        DIR = [DIRU; DIRP]

        new(fixedDofs, DIR, inletLength)
    end
end

"""@docs
    PipeBendBC

Specify the boundary conditions for the Topflow PipeBend problem (type 2)
"""
struct PipeBendBC <: TopflowBoundaryConditions
    fixedDofs::Matrix{Int64}
    DIR::Matrix{Float64}
    inletLength::Float64

    """
    Constructor
    """
    function PipeBendBC(domain::TopflowDomain, fea::TopflowFEA, Uin::Float64)
        if (mod(domain.nelx, 10) > 0 || mod(domain.nely, 10) > 0)
            throw(ArgumentError("Number of elements must be divisible by 10."))
        end

        inletLength = 2 / 10 * domain.nely
        inlet1 = 1 / 10 * domain.nely + 1
        outletLength = 2 / 10 * domain.nelx
        outlet1 = 7 / 10 * domain.nelx + 1
        nodesInlet = reshape(Int.(inlet1:(inlet1+inletLength)), 1, Int(inletLength + 1))
        nodesOutlet =
            fea.nody *
            reshape(Int.(outlet1:(outlet1+outletLength)), 1, Int(outletLength + 1))
        vec_TopBot = Int.(vcat(1:domain.nely+1:fea.nodtot, fea.nody:fea.nody:fea.nodtot))
        vec_LefRig = Int.(vcat(2:domain.nely, (fea.nodx-1)*fea.nody+2:fea.nodtot-1))
        nodesTopBot = reshape(vec_TopBot, 1, length(vec_TopBot))
        nodesLefRig = reshape(vec_LefRig, 1, length(vec_LefRig))

        fixedDofsTBx = 2 * nodesTopBot .- 1
        fixedDofsTBy = 2 * Int.(setdiff(nodesTopBot, nodesOutlet))
        fixedDofsTBy = reshape(fixedDofsTBy, 1, length(fixedDofsTBy))
        fixedDofsLRx = 2 * Int.(setdiff(nodesLefRig, nodesInlet)) .- 1
        fixedDofsLRx = reshape(fixedDofsLRx, 1, length(fixedDofsLRx))

        fixedDofsLRy = 2 * nodesLefRig
        fixedDofsInX = 2 * nodesInlet .- 1
        fixedDofsInY = 2 * nodesInlet
        fixedDofsOutX = 2 * nodesOutlet .- 1
        fixedDofsOutP = 2 * fea.nodtot .+ nodesOutlet

        fixedDofsU = [
            fixedDofsTBx fixedDofsTBy fixedDofsLRx fixedDofsLRy fixedDofsInX fixedDofsInY fixedDofsOutX
        ]
        fixedDofsP = fixedDofsOutP
        fixedDofs = [fixedDofsU fixedDofsP]

        DIRU = zeros(fea.nodtot * 2, 1)
        DIRP = zeros(fea.nodtot, 1)
        Uinlet = [-4 * x .^ 2 .+ 4 * x for x in (Uin * (0:inletLength) / inletLength)]
        DIRU[fixedDofsInX] = Uinlet'
        DIR = [DIRU; DIRP]

        new(fixedDofs, DIR, fixedDofsTBy, inletLength)
    end
end

end