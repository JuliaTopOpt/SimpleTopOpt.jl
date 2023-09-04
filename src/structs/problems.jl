module Problems

using ..ParameterDefinitions
using ..Optimizers
using ..Domains
using ..BoundaryConditionDefinitions
using ..FEMDefinitions

export Top88Problem, TophProblem, TopflowProblem, TOProblem, DoublePipeProblem, PipeBendProblem

abstract type TOProblem end

"""
Creates a Top88 problem.
"""
struct Top88Problem{U, T} <: TOProblem where {U <: Optimizer, T <: Filter}
    domain::Top88Domain
    SIMP::ModifiedSIMPParameters
    optimizer::U
    filter::T

    volfrac::Float64
    nu::Float64

    function Top88Problem(
        domain::Top88Domain,
        SIMP::ModifiedSIMPParameters,
        optimizer::U,
        filter::T,
        volfrac::Float64 = 0.4,
        nu::Float64 = 0.3,
    ) where {U <: Optimizer, T <: Filter}

        if !(volfrac > 0.0 && volfrac < 1.0)
            throw(DomainError(volfrac, "The volume fraction must be exclusively between 0 and 1."))
        end

        if SIMP.E0 == 0.0
            throw(DomainError(SIMP.E0, "SIMP parameter E_0 must be non-zero for Top88 problems"))
        elseif SIMP.Emin == 0.0
            throw(DomainError(SIMP.Emin, "SIMP parameter E_min must be non-zero for Top88 problems"))
        end
        
        new{U,T}(domain, SIMP, optimizer, filter, volfrac, nu)
    end
end

"""@docs
Creates a TopH problem
"""
struct TophProblem{U, T} <: TOProblem where {U<:Optimizer, T <: Filter}
    domain::TophDomain
    SIMP::ModifiedSIMPParameters
    optimizer::U
    filter::T

    volfrac::Float64

    function TophProblem(
        domain::TophDomain,
        simp::ModifiedSIMPParameters,
        optimizer::U,
        filter::T,
        volfrac::Float64 = 0.4, 
    ) where {U<:Optimizer, T <: Filter}
        if !(volfrac > 0.0 && volfrac < 1.0)
            throw(DomainError(volfrac, "The volume fraction must be exclusively between 0 and 1."))
        end

        new{U, T}(
            domain,
            simp,
            optimizer,
            filter,
            volfrac,
        )
    end
end

abstract type TopflowProblem <: TOProblem end

"""
Topflow Problem Type 1 -- the double pipe problem
"""
struct DoublePipeProblem{U<:Optimizer} <: TopflowProblem
    domain::TopflowDomain
    continuation::TopflowContinuation
    solver_opts::TopflowNumericals
    bkman::BrinkmanPenalizationParameters
    volfrac::Float64
    optimizer::U

    fea::TopflowFEA
    bc::DoublePipeBC

    physics::ParameterDefinitions.TopflowPhysics

    function DoublePipeProblem(
        domain::TopflowDomain,
        volfrac::Float64,
        optimizer::U,
        solver_opts::TopflowNumericals = TopflowNumericals(),
        physicals::TopflowPhysicals = TopflowPhysicals(),
    ) where {U<:Optimizer}

        @assert volfrac > 0.0 && volfrac < 1.0

        fea = TopflowFEA(domain)
        bc = DoublePipeBC(domain, fea, physicals.Uin)

        bkman = BrinkmanPenalizationParameters(physicals.mu)

        tc = TopflowContinuation(bkman, volfrac)

        physics = ParameterDefinitions.TopflowPhysics(
            physicals.Uin,
            physicals.rho,
            physicals.mu,
            physicals.Uin * (bc.inletLength * domain.Ly / domain.nely) * physicals.rho / physicals.mu
        )       

        new{U}(
            domain,
            tc,
            solver_opts,
            bkman,
            volfrac,
            optimizer,
            fea,
            bc,
            physics,
        )
    end
end


"""
Topflow Problem Type 2 -- the pipe bend problem
"""
struct PipeBendProblem{U<:Optimizer} <: TopflowProblem
    domain::TopflowDomain
    continuation::TopflowContinuation
    solver_opts::TopflowNumericals
    bkman::BrinkmanPenalizationParameters
    volfrac::Float64
    optimizer::U

    fea::TopflowFEA
    bc::PipeBendBC

    physics::ParameterDefinitions.TopflowPhysics

    Renum::Float64

    function PipeBendProblem(
        domain::TopflowDomain,
        volfrac::Float64,
        optimizer::U,
        physicals::TopflowPhysicals = TopflowPhysicals(),
        numericals::TopflowNumericals = TopflowNumericals()
    ) where {U<:Optimizer}

        @assert volfrac > 0.0 && volfrac < 1.0

        fea = TopflowFEA(domain)
        bc = PipeBendBC(domain, fea, physicals.Uin)

        bkman = BrinkmanPenalizationParameters(physicals.mu)
        tc = TopflowContinuation(bkman, volfrac)

        physics = ParameterDefinitions.TopflowPhysics(
            physicals.Uin,
            physicals.rho,
            physicals.mu,
            physicals.Uin * (bc.inletLength * domain.Ly / domain.nely) * physicals.rho / physicals.mu
        )

        new{U}(
            domain,
            tc,
            numericals,
            bkman,
            volfrac,
            optimizer,
            fea,
            bc,
            physics
        )
    end
end

end