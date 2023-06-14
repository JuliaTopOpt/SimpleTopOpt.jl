
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
struct BrinkmanPenalizationParamters <: ParameterContainer
    alphamax::Float64
    alphamin::Float64

    function BrinkmanPenalizationParamters(mu::Float64 = 1e0)
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

struct TopflowContinuation

    function TopflowContinuation(
        mu::Float64,
        xinit::Float64,
        bkman::BrinkmanPenalizationParamters,
    )


    end
end

"""
Will maintain most of the matrices, etc, required for the finite element state evaluation
"""
struct TopflowFEA
    neltot::Int64
    doftot::Int64

    nodx::Int64
    nody::Int64

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







abstract type TopflowContainer <: ProblemContainer end

"""
Topflow Problem Type 1 -- the double pipe problem
"""
struct DoublePipeContainer{U} <: TopflowContainer where {U<:OptimizerContainer}
    tfdc::TopflowDomain
    volfrac::Float64

    optimizer::U

    function DoublePipeContainer(
        tfdc::TopflowDomain,
        volfrac::Float64,
        optimizer::U,
    ) where {U<:OptimizerContainer}

        @assert volfrac > 0.0 && volfrac < 1.0

        new{U}(tfdc, volfrac, optimizer)
    end
end

"""
Topflow Problem Type 2 -- the pipe bend problem
"""
struct PipeBendContainer{U} <: TopflowContainer where {U<:OptimizerContainer}

    # fixed degrees of freedom
    fixedDofs::Int64
    # Reynolds number
    Renum::Float64
    # DIR

    optimizer::U

    function PipeBendContainer(
        fixedDofs::Int64,
        Renum::Float64,
        optimizer::U,
    ) where {U<:OptimizerContainer}

        new{U}(fixedDofs, Renum, optimizer)
    end
end
