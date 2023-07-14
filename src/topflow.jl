module TopFlow

using LinearAlgebra
using SparseArrays
using FillArrays

import ..TopflowContainer, ..TopflowOptNSParams, ..TopflowBoundaryConditions, ..TopflowFEA
export topflow

include("topflow_subroutines/analyticElement.jl")

"""
    `topflow`

Fluidic topology optimization
"""
function topflow(problem_container::T, writeout::Bool = false) where {T<:TopflowContainer}


    # TODO -- code generation?

    # TODO: do I need to error check the physical parameters, etc?
    tfdc = problem_container.tfdc
    fea = problem_container.fea
    bc = problem_container.bc
    continuation = problem_container.tc
    solver_opts = problem_container.solver_opts
    bkman = continuation.bkman

    ### Boundary conditions ctd
    # Nullspace matrices for imposing boundary conditions
    EN = Diagonal(I, fea.doftot)
    ND = copy(EN)
    ND[bc.fixedDofs, bc.fixedDofs] = 0.0
    EN -= ND

    # Vectors for free DOFs
    alldofs = 1:fea.doftot
    freedofs = setdiff(alldofs, bc.fixedDofs)

    ### Initialization
    # Solution vector
    S = zeros(fea.doftot, 1)
    dS = copy(S)
    L = copy(S)
    S[bc.fixedDofs] = DIR[bc.fixedDofs]

    # Design Field
    xPhys = volfrac * ones(tfdc.nely, tfdc.nelx)

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
        Md = 100 * (4 * sum(xPhys .* (1 .- xPhys)) / neltot)

        # Material interpolator
        alpha =
            bkman.alpgamin +
            (bkman.alphamax - bkman.alphamin) * (1 .- xPhys) ./ (1 .+ qa * xPhys)
        dalpha =
            (qa * (bkman.alphamax - bkman.alphamin) * (xPhys .- 1)) ./
            (xPhys * qa .+ 1) .^ 2 .- (bkman.alphamax - bkmin.alphamin) ./ (xPhys * qa .+ 1)

        # Newtown solve
        newton(nlittot)       # TODO -- what else does this need, if anything?

        # Objective evaluation

        obj = sum(PHI(dxv, dyv, muv, alpha', S(edofMat')))
        change = abs(objOld - obj) / objOld
        objOld = obj

        # Volume constraint
        V = mean(xPhys)

        # Print results

        # TODO: Here

        # Evaluate current iterate - continue unless considered converged

        if (change < chlim)
            chcnt += 1
        else
            chcnt = 0
        end

        if (qastep == qanum && ((chcnt == chnum) || (loopcont == conit)))
            break
        end

        # Adjoint solver
        sR = [dPHIds(dxv, dyv, muv, alpha', S[edofMat']); zeros(4, neltot)]
        RHS = sparse(iR, jR, sR)
        RHS[fixedDofs] .= 0
        sJ = JAC(dxv, dyv, muv, rhov, alpha', S[edofMat'])

        # TODO: Check size
        J = sparse(iJ, jJ, sJ)
        J = (ND' * J * ND + EN)
        L = J' \ RHS

        # Compute sensitivities
        sR = dRESdg(dxv, dyv, muv, rhov, alpha', dalpha', S[edofMat'])

        # TODO: Check size
        dRdg = sparse(iR, jE, sR)
        dphidg = dPHIdg(dxv, dyv, muv, alpha', dalpha', S[edofMat'])
        sens = reshape(dphidg - L' * dRdg, nely, nelx)

        # Volume constraint
        dV = ones(nely, nelx) ./ neltot

        # Optimality criteria update of design variables and physical densities
        # TODO -- put into its own function for modularity
        xnew = xPhys
        xlow = xPhys(:) - mvlim
        xupp = xPhys(:) + mvlim
        ocfac = xPhys(:) .* max(1e-10, (-sens(:) ./ dV(:))) .^ (1 / 3)
        l1 = 0
        l2 = (1 / (neltot * volfrac) * sum(ocfac))^3
        while (l2 - l1) / (l1 + l2) > 1e-3
            lmid = 0.5 * (l2 + l1)
            xnew(:) = max(0, max(xlow, min(1, min(xupp, ocfac / (lmid^(1 / 3))))))
            if mean(xnew(:)) > volfrac
                l1 = lmid
            else
                l2 = lmid
            end
        end
        xPhys = xnew


        # Continuation update

        if (qastep < qanum && (loopcont == conit || chcnt == chnum))
            loopcont = 0
            chcnt = 0
            qastep = qastep + 1
            qa = qavec(qastep)
        end

    end

    if writeout


    end
end


"""
Nonlinear Newton Solver
"""
function newton(nlittot::Int, printout::Bool = false)

    fail = -1
    normR = 1
    nlit = 0

    while fail != 1
        pc.solver_opts.nlit += 1
        nlittot += 1

        # Build residual and Jacobian
        # TODO -- need RES
        sR = RES(dxv, dyv, muv, rhov, alpha', S(fea.edofMat'))
        R = sparse(iR, jR, sR)
        R(fixedDofs) = 0

        if nlit == 1
            r0 = norm(R)
        end
        r1 = norm(R)
        normR = r1 / r0

        # TODO -- PLOTTING LINE HERE

        if normR < solver_opts.nltol
            break
        end

        sJ = JAC(dxv, dyv, muv, rhov, alpha', S[fea.edofMat'])
        J = sparse(fea.iJ, fea.jJ, sJ)
        J = (ND' * J * ND + EN)

        # Calculate Newton step
        dS = -J \ R

        # L2-norm line search
        Sp = S + 0.5 * dS
        sR = RES(dxv, dyv, muv, rhov, alpha', Sp[fea.edofMat'])
        R = sparse(fea.iR, fea.jR, sR)
        R[bc.fixedDofs] .= 0
        r3 = norm(R)

        # Solution update with "optimal" damping
        lambda = max(0.01, min(1.0, (3 * r1 + r3 - 4 * r2) / (4 * r1 + 4 * r3 - 8 * r2)))
        S = S + lambda * dS

        # if fail, retry from zero solution
        if (nlit == nlmax && fail < 0)
            nlit = 0
            S[freeDofs] = 0.0
            normR = 1
            fail += 1
        end

        if (nlit == nlmax && fail < 1)
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
function OCUpdate(
    xPhys::Matrix{Float64},
    solver_opts::TopflowOptNSParams,
    fea::TopflowFEA,
)

    xlow = xPhys(:) - mvlim
    xupp = xPhys(:) + mvlim
    ocfac = xPhys(:) .* max(1e-10, (-sens(:) ./ dV(:))) .^ (1 / 3)
    l1 = 0
    l2 = (1 / (neltot * volfrac) * sum(ocfac))^3
    while (l2 - l1) / (l1 + l2) > 1e-3
        lmid = 0.5 * (l2 + l1)
        xnew(:) = max(0, max(xlow, min(1, min(xupp, ocfac / (lmid^(1 / 3))))))
        if mean(xnew(:)) > volfrac
            l1 = lmid
        else
            l2 = lmid
        end
    end

    return xPhys

end




end
