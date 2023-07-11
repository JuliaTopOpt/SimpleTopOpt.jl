
abstract type OptimizerContainer end
abstract type ParameterContainer end
abstract type ProblemContainer end
abstract type FiniteElementContainer end

##################################################################################################
# BEGIN HYPERPARAMETER STRUCT DEFINITIONS
##################################################################################################

"""
For the optimality criterion optimizer.
"""
struct OCParameters <: OptimizerContainer
    max_iter::Int64
    mvlim::Float64

    function OCParameters(max_iter::Int64, mvlim::Float64)
        @assert max_iter > 0
        @assert mvlim > 0.0

        new(max_iter, mvlim)
    end
end

"""
Modified SIMP law parameter container
"""
struct SIMPParameters <: ParameterContainer
    rmin::Float64
    penal::Float64

    function SIMPParameters(rmin::Float64, penal::Float64)
        @assert rmin > 0
        @assert penal > 0

        new(rmin, penal)
    end
end

"""
Brinkman penalization parameter container
"""
struct BrinkmanPenalizationParameters <: ParameterContainer
    alphamax::Float64
    alphamin::Float64

    function BrinkmanPenalizationParameters(mu::Float64)
        new(2.5 * mu / (0.01^2), 2.5 * mu / (100^2))
    end
end

##################################################################################################
# BEGIN PROBLEM STRUCT DEFINITIONS
##################################################################################################

"""
Enough to specify a problem to Top88
"""
#struct Top88Container <: ProblemContainer

#end

"""
Enough to specify a problem to Toph
"""
# struct TophContainer{T, S, U} <: ProblemContainer where {T<:AbstractFloat, S<:Integer, U<:OptimizerContainer}
#     nelx::S
#     nely::S

#     volfrac::T

#     simp::SIMPParameters

#     optimizer::U

#     function TophContainer(
#         nelx::S, 
#         nely::S, 
#         volfrac::T, 
#         simp::SIMPParameters,
#         optimizer::U    
#     ) where {T<:AbstractFloat, S<:Integer, U<:OptimizerContainer}

#         @assert nelx > 1
#         @assert nely > 1
#         @assert volfrac > 0.0 && volfrac < 1.0

#         new{T, S, U}(
#             nelx,
#             nely,
#             volfrac,
#             simp,
#             optimizer
#         )

#     end
# end

# BEGIN TOPFLOW STRUCT DEFINITIONS
"""
For defining the domain; still assuming quad elements and rectangular domain
"""
struct TopflowDomain
    Lx::Float64
    Ly::Float64

    dx::Float64
    dy::Float64

    nely::Int64
    nelx::Int64

    function TopflowDomain(Lx::Float64 = 1.0, Ly::Float64 = 1.0, nely::Int64 = 30)
        @assert Lx > 0.0
        @assert Ly > 0.0
        @assert nely > 1

        nelx = nely * Lx / Ly
        dx = Lx / nelx
        dy = Ly / nely

        new(Lx, Ly, dx, dy, nely, nelx)
    end
end

"""
Topflow continuation strategy
"""
struct TopflowContinuation <: ParameterContainer
    ainit::Float64
    qinit::Float64
    qavec::Matrix{Float64}
    qanum::Int64
    conit::Int64

    """
    Constructor
    """
    function TopflowContinuation(
        volfrac::Float64, # = 1/3 by default
        bkman::BrinkmanPenalizationParameters,
        conit::Int64 = 50,
    )
        ainit = bkman.alphamax / 100

        qinit =
            (-volfrac * (bkman.alphamax - bkman.alphamin) - ainit + bkman.alphamax) /
            (volfrac * (ainit - bkman.alphamin))
        qavec = qinit ./ [1 2 10 20]
        qanum = length(qavec)

        new(ainit, qinit, qavec, qanum, conit)

    end
end

"""
Topflow optimization and Newton solver parameters
"""
struct TopflowOptNSParams <: ParameterContainer

    maxiter::Int64
    mvlim::Float64
    # plotdes -- TODO: this is just a plotting utility switch
    chlim::Float64
    chnum::Int64

    # Newton Solver Parameters
    nltol::Float64
    nlmax::Int64
    # plotres -- TODO: this is just a plotting utility switch

    # TODO -- export options?

    """
    Constructor
    """
    function TopflowOptNSParams(
        maxiter::Int64,
        mvlim::Float64,
        chlim::Float64,
        chnum::Int64,
        nltol::Float64,
        nlmax::Int64,
    )
        # TODO: assertions on positivty/ non-negativity, etc.
        @assert maxiter > 0
        @assert mvlim > 0
        @assert chlim > 0
        @assert chnum > 0
        @assert nltol > 0
        @assert nlmax > 0

        new(maxiter, mvlim, chlim, chnum, nltol, nlmax)
    end
end

"""
Will maintain most of the matrices, etc, required for the finite element state evaluation
"""
struct TopflowFEA
    neltot::Int64   # Total number of elements
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

    function TopflowFEA(tfdc::TopflowDomain)
        """
        Assembles most of the finite element matrix problem.
        """

        nodx = tfdc.nelx + 1
        nody = tfdc.nely + 1
        nodtot = nodx * nody
        neltot = tfdc.nelx * tfdc.nely
        doftot = 3 * nodtot

        nodenrs = reshape(1:nodtot, nody, nodx)
        edofVecU = reshape(2 * nodenrs[1:end-1, 1:end-1] .+ 1, neltot, 1)
        edofMatU =
            repeat(edofVecU, 1, 8) +
            repeat([0 1 (2 * tfdc.nely .+ [2 3 0 1]) -2 -1], neltot, 1)
        edofVecP = reshape(nodenrs[1:end-1, 1:end-1], neltot, 1)
        edofMatP = repeat(edofVecP, 1, 4) + repeat([1 (tfdc.nely .+ [2 1]) 0], neltot, 1)

        edofMat = [edofMatU (2 * nodtot .+ edofMatP)]

        iJ = reshape(kron(edofMat, ones(Int64, 12, 1))', 144 * neltot, 1)
        jJ = reshape(kron(edofMat, ones(Int64, 1, 12))', 144 * neltot, 1)
        iR = reshape(edofMat', 12 * neltot, 1)
        jR = ones(12 * neltot, 1)
        jE = repeat((1:neltot)', 12, 1)

        new(neltot, doftot, nodx, nody, nodtot, edofMat, iJ, jJ, iR, jR, jE)
    end
end

abstract type TopflowBoundaryConditions end


struct DoublePipeBC <: TopflowBoundaryConditions
    fixedDofs::Matrix{Int64}
    DIR::Matrix{Float64}

    inletLength::Float64

    """
    Constructor
    """
    function DoublePipeBC(tfdc::TopflowDomain, fea::TopflowFEA, Uin::Float64)
        if mod(tfdc.nely, 6) > 0
            throw(ArgumentError("Number of elements in y-dir must be divisible by 6."))
        elseif tfdc.nely < 30
            throw(ArgumentError("Number of elements in y-dir must be at least 30."))
        end

        inletLength = 1 / 6 * tfdc.nely
        inlet1 = 1 / 6 * tfdc.nely + 1
        inlet2 = 4 / 6 * tfdc.nely + 1
        outletLength = 1 / 6 * tfdc.nely
        outlet1 = 1 / 6 * tfdc.nely + 1
        outlet2 = 4 / 6 * tfdc.nely + 1

        nodesInlet = Int.([(inlet1:(inlet1+inletLength))' (inlet2:(inlet2+inletLength))'])
        nodesOutlet =
            (fea.nodx - 1) * fea.nody .+
            Int.([(outlet1:(outlet1+outletLength))' (outlet2:(outlet2+outletLength))'])

        # TEMP -- there may be a more succinct way of getting this?

        vec_TopBot = Int.(vcat(1:tfdc.nely+1:fea.nodtot, fea.nody:fea.nody:fea.nodtot))
        vec_LefRig = Int.(vcat(2:tfdc.nely, (tfdc.nelx)*fea.nody+2:fea.nodtot-1))

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

abstract type TopflowContainer <: ProblemContainer end

"""
Topflow Problem Type 1 -- the double pipe problem
"""
struct DoublePipeContainer{U<:OptimizerContainer} <: TopflowContainer
    tfdc::TopflowDomain
    tc::TopflowContinuation
    solver_opts::TopflowOptNSParams # TODO: give this a better name
    bkman::BrinkmanPenalizationParameters
    volfrac::Float64
    optimizer::U

    fea::TopflowFEA
    bc::DoublePipeBC

    Uin::Float64
    rho::Float64
    mu::Float64

    Renum::Float64

    function DoublePipeContainer(
        tfdc::TopflowDomain,
        volfrac::Float64,
        optimizer::U,
        Uin::Float64 = 1e0,
        rho::Float64 = 1e0,
        mu::Float64 = 1e0,
    ) where {U<:OptimizerContainer}

        @assert volfrac > 0.0 && volfrac < 1.0

        fea = TopflowFEA(tfdc)
        bc = DoublePipeBC(tfdc, fea, Uin)

        # TODO: fill this in; or have taken as argument
        solver_opts = Nothing       

        bkman = BrinkmanPenalizationParameters(mu)

        tc = TopflowContinuation(volfrac, bkman)

        Renum = Uin * (bc.inletLength * tfdc.Ly / tfdc.nely) * rho / mu

        new{U}(tfdc, tc, solver_opts, bkman, volfrac, optimizer, fea, bc, Renum)
    end
end


"""
BCs for problem 2
"""
struct PipeBendBC <: TopflowBoundaryConditions
    fixedDofs::Matrix{Int64}
    DIR::Matrix{Float64}

    # TEMP
    fixedDofsTBy::Matrix{Int64}

    inletLength::Float64

    """
    Constructor
    """
    function PipeBendBC(tfdc::TopflowDomain, fea::TopflowFEA, Uin::Float64)
        if (mod(tfdc.nelx, 10) > 0 || mod(tfdc.nely, 10) > 0)
            throw(ArgumentError("Number of elements must be divisible by 10."))
        end
        # TODO -- sizes may be problematic; will need to transpose?

        ## DEBUG:
        inletLength = 2 / 10 * tfdc.nely
        inlet1 = 1 / 10 * tfdc.nely + 1
        outletLength = 2 / 10 * tfdc.nelx
        outlet1 = 7 / 10 * tfdc.nelx + 1
        nodesInlet = reshape(Int.(inlet1:(inlet1+inletLength)), 1, Int(inletLength + 1))
        nodesOutlet =
            fea.nody *
            reshape(Int.(outlet1:(outlet1+outletLength)), 1, Int(outletLength + 1))
        vec_TopBot = Int.(vcat(1:tfdc.nely+1:fea.nodtot, fea.nody:fea.nody:fea.nodtot))
        vec_LefRig = Int.(vcat(2:tfdc.nely, (fea.nodx-1)*fea.nody+2:fea.nodtot-1))
        nodesTopBot = reshape(vec_TopBot, 1, length(vec_TopBot))
        nodesLefRig = reshape(vec_LefRig, 1, length(vec_LefRig))

        # TEMP
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

        println("Sizes for problem 2 -- pipe bend")
        println(size(fixedDofsTBx))
        println(size(fixedDofsTBy))
        println(size(fixedDofsLRx))
        println(size(fixedDofsLRy))
        println(size(fixedDofsInX))
        println(size(fixedDofsInY))
        println(size(fixedDofsOutX))
        println(size(fixedDofsOutP))

        for i in eachindex(fixedDofsTBy)
            print(string(i) * ": ")
            println(fixedDofsTBy[i])
        end

        println(sum(fixedDofsTBy))

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



"""
Topflow Problem Type 2 -- the pipe bend problem
"""
struct PipeBendContainer{U<:OptimizerContainer} <: TopflowContainer
    tfdc::TopflowDomain
    tc::TopflowContinuation
    solver_opts::TopflowOptNSParams # TODO Give this a better name

    volfrac::Float64
    optimizer::U

    fea::TopflowFEA
    bc::PipeBendBC

    Uin::Float64
    rho::Float64
    mu::Float64

    Renum::Float64

    function PipeBendContainer(
        tfdc::TopflowDomain,
        tc::TopflowContinuation,
        volfrac::Float64,
        optimizer::U,
        Uin::Float64 = 1e0,
        rho::Float64 = 1e0,
        mu::Float64 = 1e0,
    ) where {U<:OptimizerContainer}

        @assert volfrac > 0.0 && volfrac < 1.0

        fea = TopflowFEA(tfdc)
        bc = PipeBendBC(tfdc, fea, Uin)

        Renum = Uin * (bc.inletLength * tfdc.Ly / tfdc.nely) * rho / mu

        new{U}(tfdc, tc, volfrac, optimizer, fea, bc, Uin, rho, mu, Renum)
    end
end
