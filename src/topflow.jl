module TopFlow

include("topflow_subroutines/analyticElement.jl")

using LinearAlgebra
using SparseArrays
using FillArrays
# using topflow_MATLAB_functions
using MATLAB
using Statistics

import ..TopflowContainer, ..TopflowOptNSParams, ..TopflowBoundaryConditions, ..TopflowFEA
export topflow

const _dir = @__DIR__

# TODO -- revisit this after plots; instead of circumventing symbolics entirely, thought is to use MATLAB as the integrator
# TODO -- another path = numeric integration over symbolic integration/ burning time with matlab
"""
    call_JAC

The following is a wrapper around the MATLAB interface call
into JAC
"""
function call_JAC(dxv, dyv, muv, rhov, alpha, in6)
    eval_string("addpath('" * _dir * "/topflow_subroutines/matlab_ffi_routines')")

    out = mxcall(:JAC, 1, mxarray(dxv), mxarray(dyv), mxarray(muv), mxarray(rhov), mxarray(alpha), mxarray(in6))

    return out
end

"""
    `call_RES`

The following is a wrapper around the MATLAB interface call
into RES
"""
function call_RES(dxv, dyv, muv, rhov, alpha, in6)
    eval_string("addpath('" * _dir * "/topflow_subroutines/matlab_ffi_routines')")
 
    mdxv = mxarray(dxv)
    mdyv = mxarray(dyv)
    mmuv = mxarray(muv)
    mrhov = mxarray(rhov)
    malpha = mxarray(alpha)
    min6 = mxarray(in6)

    out = mxcall(:RES, 1, mdxv, mdyv, mmuv, mrhov, malpha, min6)

    return out
end

"""
    call_PHI

The following is a wrapper around the MATLAB interface call into PHI
"""
function call_PHI(dxv, dyv, muv, alpha, in5)
    eval_string("addpath('" * _dir * "/topflow_subroutines/matlab_ffi_routines')")
    out = mxcall(:PHI, 1, dxv, dyv, muv, mxarray(alpha), mxarray(in5))

    return out
end

"""
    call_dPHIdg

The following is a wrapper around the MATLAB interface call into dPHIdg
"""
function call_dPHIdg(dxv, dyv, muv, alpha, dalpha, in6)
    eval_string("addpath('" * _dir * "/topflow_subroutines/matlab_ffi_routines')")
    out = mxcall(:dPHIdg, 1, dxv, dyv, muv, mxarray(alpha), mxarray(dalpha), mxarray(in6))

    return out
end

"""
    call_dPHIds

The following is a wrapper around the MATLAB interface call into dPHIds
"""
function call_dPHIds(dxv, dyv, muv, alpha, in5)
    eval_string("addpath('" * _dir * "/topflow_subroutines/matlab_ffi_routines')")
    out = mxcall(:dPHIds, 1, dxv, dyv, muv, mxarray(alpha), mxarray(in5))

    return out
end

"""
    call_dRESdg

The following is a wrapper around the MATLAB interface call into dRESdg
"""
function call_dRESdg(dxv, dyv, muv, rhov, alpha, dalpha, in7)
    eval_string("addpath('" * _dir * "/topflow_subroutines/matlab_ffi_routines')")
    out = mxcall(
        :dRESdg,
        1,
        dxv,
        dyv,
        muv,
        rhov,
        mxarray(alpha),
        mxarray(dalpha),
        mxarray(in7),
    )

    return out
end


"""
    `topflow`

Fluidic topology optimization
"""
function topflow(problem_container::T, writeout::Bool = false) where {T<:TopflowContainer}

    # TODO -- once analyticElement is replaced with pure Julia
    # JAC, RES, PHI, dPHIdg, dPHIds, dRESdg = analyticElement.generation()

    # TODO: do I need to error check the physical parameters, etc?
    tfdc = problem_container.tfdc
    fea = problem_container.fea
    bc = problem_container.bc
    continuation = problem_container.tc
    solver_opts = problem_container.solver_opts
    bkman = continuation.bkman
    vf = problem_container.volfrac
    μ = problem_container.mu
    ρ = problem_container.rho

    ### Boundary conditions ctd
    # Nullspace matrices for imposing boundary conditions
    EN = Diagonal(I, fea.doftot)
    ND = copy(EN)
    ND[bc.fixedDofs, bc.fixedDofs] .= 0.0
    EN -= ND

    # Vectors for free DOFs
    alldofs = 1:fea.doftot
    freedofs = setdiff(alldofs, bc.fixedDofs)

    ### Initialization
    # Solution vector
    S = zeros(fea.doftot, 1)
    dS = copy(S)
    L = copy(S)
    S[bc.fixedDofs] = bc.DIR[bc.fixedDofs]

    # Design Field
    xPhys = vf * ones(tfdc.nely, tfdc.nelx)

    # Counters
    loop = 0
    loopcont = 0
    nlittot = 0
    chcnt = 0

    # Change
    change = Inf
    objOld = Inf

    # Continuation
    qastep = 1
    qa = continuation.qavec[1]

    # Vectorized constants 
    dxv = tfdc.dx * ones(1, fea.neltot)
    dyv = tfdc.dy * ones(1, fea.neltot)
    muv = problem_container.mu * ones(1, fea.neltot)
    rhov = problem_container.rho * ones(1, fea.neltot)

    ### Output
    if writeout
        println("TODO --- fill out problem info string")
    end

    ### Begin optimization loop
    while loop <= solver_opts.maxiter

        # Greyscale indicator
        Md = 100 * (4 * sum(xPhys .* (1 .- xPhys)) / fea.neltot)

        # Material interpolator
        alpha =
            bkman.alphamin .+
            (bkman.alphamax - bkman.alphamin) * (1 .- xPhys[:]) ./ (1 .+ qa * xPhys[:])
        dalpha =
            (qa * (bkman.alphamax - bkman.alphamin) * (xPhys[:] .- 1)) ./
            (xPhys[:] * qa .+ 1) .^ 2 .- (bkman.alphamax - bkman.alphamin) ./ (xPhys[:] * qa .+ 1)

        alpha_T = collect(alpha')
        dalpha_T = collect(dalpha')

        # Newton solve
        newton(xPhys, nlittot, dxv, dyv, muv, rhov, alpha, S, fea, bc, solver_opts, ND, EN)

        # Objective evaluation
        obj = sum(call_PHI(dxv, dyv, muv, alpha_T, S[fea.edofMat']))
        change = abs(objOld - obj) / objOld
        objOld = obj

        # Volume constraint
        V = mean(xPhys)

        # Print results

        # TODO: Here

        # Evaluate current iterate - continue unless considered converged

        if (change < solver_opts.chlim)
            chcnt += 1
        else
            chcnt = 0
        end

        if (qastep == continuation.qanum && ((chcnt == solver_opts.chnum) || (loopcont == continuation.conit)))
            break
        end

        # Adjoint solver
        sR = [call_dPHIds(dxv, dyv, muv, alpha_T, S[fea.edofMat']); zeros(4, fea.neltot)]

        RHS = sparse(fea.iR[:], fea.jR[:], sR[:])

        RHS[bc.fixedDofs] .= 0
        sJ = call_JAC(dxv, dyv, muv, rhov, alpha_T, S[fea.edofMat'])

        # TODO: Check size
        J = sparse(fea.iJ[:], fea.jJ[:], sJ[:])
        J = (ND' * J * ND + EN)
        L = Matrix(J') \ Matrix(RHS)

        # Compute sensitivities
        sR = call_dRESdg(dxv, dyv, muv, rhov, alpha_T, dalpha_T, S[fea.edofMat'])

        # TODO: Check size
        dRdg = sparse(fea.iR[:], fea.jE[:], sR[:])
        dphidg = call_dPHIdg(dxv, dyv, muv, alpha_T, dalpha_T, S[fea.edofMat'])
        sens = reshape(dphidg - L' * dRdg, tfdc.nely, tfdc.nelx)

        # Volume constraint
        dV = ones(tfdc.nely, tfdc.nelx) ./ fea.neltot

        # Optimality criteria update of design variables and physical densities
        # TODO -- put into its own function for modularity
        xnew = xPhys
        xlow = xPhys[:] .- solver_opts.mvlim
        xupp = xPhys[:] .+ solver_opts.mvlim
        ocfac = xPhys[:] .* max.(1e-10, (-sens[:] ./ dV[:])) .^ (1 / 3)
        l1 = 0
        l2 = (1 / (fea.neltot * vf) * sum(ocfac))^3
        while (l2 - l1) / (l1 + l2) > 1e-3
            lmid = 0.5 * (l2 + l1)
            xnew[:] = max.(0, max.(xlow, min.(1, min.(xupp, ocfac / (lmid^(1 / 3))))))
            if mean(xnew[:]) > vf
                l1 = lmid
            else
                l2 = lmid
            end
        end
        xPhys = xnew


        # Continuation update

        if (qastep < continuation.qanum && (loopcont == continuation.conit || chcnt == solver_opts.chnum))
            loopcont = 0
            chcnt = 0
            qastep = qastep + 1
            qa = continuation.qavec[qastep]
        end

    end


    # TODO -- what to return?
    #      -- will need at least xPhys

    return xPhys
end


"""
Nonlinear Newton Solver
"""
function newton(
    xPhys,
    nlittot::Int,
    dxv,
    dyv,
    muv,
    rhov,
    alpha,
    S,
    fea::TopflowFEA,
    bc,
    solver_opts,
    ND,
    EN,
    printout::Bool = false,
)

    fail = -1
    normR = 1
    nlit = 0

    alpha_T = collect(alpha')

    r0 = 0.0
    while fail != 1
        nlit += 1
        nlittot += 1

        # Build residual and Jacobian
        sR = call_RES(dxv, dyv, muv, rhov, alpha_T, S[fea.edofMat'])
        R = sparse(fea.iR[:], fea.jR[:], sR[:])
        R[bc.fixedDofs] .= 0

        if nlit == 1
            r0 = norm(R)
        end
        r1 = norm(R)
        normR = r1 / r0

        # TODO -- PLOTTING LINE HERE

        if normR < solver_opts.nltol
            break
        end

        sJ = call_JAC(dxv, dyv, muv, rhov, alpha_T, S[fea.edofMat'])
        J = sparse(fea.iJ[:], fea.jJ[:], sJ[:])
        J = (ND' * J * ND + EN)

        # Calculate Newton step
        # TODO -- J and R are both sparse matrices, but doing the below cast is expensive
        #         and we completely lose out on having these sparse matrices to begin with
        dS = -Matrix(J) \ Matrix(R)

        # L2-norm line search
        Sp = S + 0.5 * dS
        sR = call_RES(dxv, dyv, muv, rhov, alpha_T, Sp[fea.edofMat'])
        R = sparse(fea.iR[:], fea.jR[:], sR[:])
        R[bc.fixedDofs] .= 0
        r2 = norm(R)

        Sp = S + 1.0*dS;
        sR = call_RES(dxv, dyv, muv, rhov, alpha_T, Sp[fea.edofMat']);
        R = sparse(fea.iR[:], fea.jR[:], sR[:])
        R[bc.fixedDofs] .= 0
        r3 = norm(R)

        # Solution update with "optimal" damping
        lambda = max(0.01, min(1.0, (3 * r1 + r3 - 4 * r2) / (4 * r1 + 4 * r3 - 8 * r2)))
        S = S + lambda * dS

        # if fail, retry from zero solution
        if (nlit == solver_opts.nlmax && fail < 0)
            nlit = 0
            S[freeDofs] = 0.0
            normR = 1
            fail += 1
        end

        if (nlit == solver_opts.nlmax && fail < 1)
            fail += 1
        end

        # TODO -- what to return?

    end

    if printout
        println("<Something informative here>")
    end

    if fail == 1
        println(
            "Newton solver did not converge after retry from zero! Stopping optimization.",
        )
    end

    return xPhys

end


"""
Optimality Criterion optimization update step
"""
function OCUpdate(problem_container::TopflowContainer)
    ## TODO -- how to type annotate where we should accept either problem container but
    ##         only ones using OCParameters

    fea = problem_container.TopflowFEA
    bc = problem_container.bc

    oc = problem_container.optimizer
    vf = problem_container.volfrac


    xlow = xPhys[:] - mvlim
    xupp = xPhys[:] + mvlim
    # TODO -- where are sens and dV coming from?
    ocfac = xPhys[:] .* max(1e-10, (-sens[:] ./ dV[:])) .^ (1 / 3)
    l1 = 0
    l2 = (1 / (fea.neltot * vf) * sum(ocfac))^3
    while (l2 - l1) / (l1 + l2) > 1e-3
        lmid = 0.5 * (l2 + l1)
        xnew[:] = max(0, max(xlow, min(1, min(xupp, ocfac / (lmid^(1 / 3))))))
        if mean(xnew[:]) > vf
            l1 = lmid
        else
            l2 = lmid
        end
    end

    return xPhys
end




end





if abspath(PROGRAM_FILE) == @__FILE__
    Lx = 1.0
    Ly = 1.0
    nely = 30
    volfrac = 1 / 3
    Uin = 1e0
    rho = 1e0
    mu = 1e0
    conit = 50

    tfdc = TopflowDomain(Lx, Ly, nely)
    fea = SimpleTopOpt.TopflowFEA(tfdc)
    optimizer = OCParameters(200, 0.2)
    dpbc = SimpleTopOpt.DoublePipeBC(tfdc, fea, Uin)
    dpc = DoublePipeContainer(tfdc, volfrac, optimizer, Uin, rho, mu)

    t_xPhys = SimpleTopOpt.TopFlow.topflow(dpc)
end
