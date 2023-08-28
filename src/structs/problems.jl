
abstract type TOProblem end

"""
Creates a Top88 problem.
"""
struct Top88Problem <: TOProblem
    domain::Top88Domain
    SIMP::SIMPParameters

    volfrac::Float64

    nu::Float64
    use_sensitivity::Bool

    function Top88Container(
        t8dc::Top88Domain,
        SIMP::SIMPParameters,
        volfrac::Float64 = 0.5,
        nu::Float64 = 0.3,
        use_sensitivity::Bool = true,
    )

        if SIMP.E0 == 0.0
            throw(DomainError(SIMP.E0, "SIMP parameter E_0 must be non-zero for Top88 problems"))
        elseif SIMP.Emin == 0.0
            throw(DomainError(SIMP.Emin, "SIMP parameter E_min must be non-zero for Top88 problems"))
        end
        

        new(t8dc, SIMP, volfrac, nu, use_sensitivity)
    end
end

"""
Enough to specify a problem to Toph
"""
struct TophProblem{U} <: TOProblem where {U<:Optimizer}
    domain::TophDomain
    volfrac::Float64

    simp::SIMPParameters

    optimizer::U

    function TophContainer(
        thdc::TophDomain,
        volfrac::Float64, 
        simp::SIMPParameters,
        optimizer::U    
    ) where {U<:Optimizer}

        @assert nelx > 1
        @assert nely > 1
        @assert volfrac > 0.0 && volfrac < 1.0

        new{U}(
            thdc,
            volfrac,
            simp,
            optimizer
        )

    end
end


abstract type TopflowProblem <: TOProblem end


"""
Topflow Problem Type 1 -- the double pipe problem
"""
struct DoublePipeProblem{U<:Optimizer} <: TopflowProblem
    domain::TopflowDomain
    tc::TopflowContinuation
    solver_opts::TopflowOptNSParams # TODO: give this a better name
    bkman::BrinkmanPenalizationParameters
    volfrac::Float64
    optimizer::U

    fea::TopflowFEA
    bc::DoublePipeBC

    physicals::TopflowPhysics

    function DoublePipeContainer(
        domain::TopflowDomain,
        volfrac::Float64,
        optimizer::U,
        OptNSParams::TopflowOptNSParams,
        physicals::TopflowPhysics,
    ) where {U<:OptimizerContainer}

        @assert volfrac > 0.0 && volfrac < 1.0

        fea = TopflowFEA(domain)
        bc = DoublePipeBC(domain, fea, Uin)


        bkman = BrinkmanPenalizationParameters(physicals.mu)

        tc = TopflowContinuation(volfrac, bkman)

        physicals.Renum = Uin * (bc.inletLength * domain.Ly / domain.nely) * rho / mu

        new{U}(
            domain,
            tc,
            OptNSParams,
            bkman,
            volfrac,
            optimizer,
            fea,
            bc,
            physicals,
        )
    end
end


"""
Topflow Problem Type 2 -- the pipe bend problem
"""
struct PipeBendProblem{U<:Optimizer} <: TopflowProblem
    domain::TopflowDomain
    tc::TopflowContinuation
    solver_opts::TopflowNonlinearNewt
    bkman::BrinkmanPenalizationParameters
    volfrac::Float64
    optimizer::U

    fea::TopflowFEA
    bc::PipeBendBC

    physics::TopflowPhysics

    Renum::Float64

    function PipeBendContainer(
        domain::TopflowDomain,
        volfrac::Float64,
        optimizer::U,
        physicals::TopflowPhysics = TopflowPhysics(),

    ) where {U<:Optimizer}

        @assert volfrac > 0.0 && volfrac < 1.0

        fea = TopflowFEA(domain)
        bc = PipeBendBC(domain, fea, Uin)

        solver_opts = TopflowOptNSParams()
        bkman = BrinkmanPenalizationParameters(mu)
        tc = TopflowContinuation(volfrac, bkman)

        physicals.Renum = Uin * (bc.inletLength * domain.Ly / domain.nely) * rho / mu

        new{U}(
            domain,
            tc,
            solver_opts,
            bkman,
            volfrac,
            optimizer,
            fea,
            bc,
            physicals
        )
    end
end
