module Top88

using LinearAlgebra
using SparseArrays
using Statistics

using ..Structs


export optimize

"""
    optimize

A direct, naive Julia port of Andreassen et al. "Efficient topology optimization in MATLAB
using 88 lines of code." By default, this will reproduce the optimized MBB beam from Sigmund
(2001).

"""
function optimize(
    problem::Top88Problem,
    writeout::Bool = false,
    loop_max::Int = 1000,
)::Top88Solution

    nelx = problem.domain.nelx
    nely = problem.domain.nely

    rmin = problem.filter.rmin
    penal = problem.SIMP.penal
    E0 = problem.SIMP.E0
    Emin = problem.SIMP.Emin

    nu = problem.nu
    volfrac = problem.volfrac

    KE = Top88FEA(nu)

    nodenrs = reshape(1:(1+nelx)*(1+nely), 1 + nely, 1 + nelx)
    edofVec = reshape(2 * nodenrs[1:end-1, 1:end-1] .+ 1, nelx * nely, 1)
    edofMat = zeros(Int64, nelx * nely, 8)

    offsets = [0 1 2 * nely .+ [2 3 0 1] -2 -1]
    for i = 1:8
        for j = 1:nelx*nely
            edofMat[j, i] = edofVec[j] + offsets[i]
        end
    end

    iK = reshape(kron(edofMat, ones(8, 1))', 64 * nelx * nely, 1)
    jK = reshape(kron(edofMat, ones(1, 8))', 64 * nelx * nely, 1)

    F = spzeros(2 * (nely + 1) * (nelx + 1))
    F[2, 1] = -1
    U = spzeros(2 * (nely + 1) * (nelx + 1))

    fixeddofs = union(1:2:2*(nely+1), [2 * (nelx + 1) * (nely + 1)])
    alldofs = 1:2*(nely+1)*(nelx+1)
    freedofs = setdiff(alldofs, fixeddofs)

    # Prepare the filter
    H, Hs = prepare_filter(problem.domain, rmin)

    # Initialize iteration
    x = volfrac * ones(nely, nelx)
    xPhys = x
    loop = 0
    change = 1
    cValues = []

    # Start iteration
    while change > 0.01
        loop += 1
        # FE-Analysis
        sK = [j * ((i + Emin)^penal) for i in ((E0 - Emin) * xPhys[:]') for j in KE[:]]
        K = sparse(iK[:], jK[:], sK)
        K = (K + K') / 2

        KK = cholesky(K[freedofs, freedofs])
        U[freedofs] = KK \ F[freedofs]

        mat = (U[edofMat] * KE) .* U[edofMat]

        # Objective function and sensitivity analysis
        ce = reshape([sum(mat[i, :]) for i = 1:(size(mat)[1])], nely, nelx)
        c = sum(sum((Emin * ones(size(xPhys)) .+ (xPhys .^ penal) * (E0 - Emin)) .* ce))
        push!(cValues, c)
        dc = -penal * (E0 - Emin) * xPhys .^ (penal - 1) .* ce
        dv = ones(nely, nelx)

        # Filtering/ modification of sensitivities
        dc, dv = filter_implementation(problem.filter, x, dc, dv, H, Hs)

        # Optimality criteria update of design variables and physical densities
        xnew = OC(nelx, nely, x, volfrac, dc, dv, H, Hs, problem.filter)
        change = maximum(abs.(x - xnew))
        x = xnew

        if writeout
            println(
                "Loop = ",
                loop,
                ", Change = ",
                change,
                ", c = ",
                c,
                ", structural density = ",
                mean(x),
            )
        end
        loop >= loop_max && break
    end

    converged = change <= 0.01

    return Top88Solution(
        x,
        converged,
        loop,
    )
end


function Top88FEA(nu::Float64)
    A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12]
    A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6]
    B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4]
    B12 = [2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2]
    KE = 1 / (1 - nu^2) / 24 * ([A11 A12; A12' A11] + nu * [B11 B12; B12' B11])

    return KE
end


"""

"""
function filter_implementation(s::SensitivityFilter, x, dc::Matrix{Float64}, dv::Matrix{Float64}, H, Hs)
    dc[:] = H * (x[:] .* dc[:]) ./ Hs ./ max(1e-3, maximum(x[:]))

    return dc, dv
end

"""

"""
function filter_implementation(d::DensityFilter, x, dc::Matrix{Float64}, dv::Matrix{Float64}, H, Hs)
    dc[:] = H * (dc[:] ./ Hs)
    dv[:] = H * (dv[:] ./ Hs)

    return dc, dv
end

"""
Prepare sensitivity/ density filter
"""
function prepare_filter(domain::Top88Domain, rmin::Float64)
    
    nelx = domain.nelx
    nely = domain.nely

    iH = ones(nelx * nely * (2 * (convert(Int64, ceil(rmin) - 1)) + 1)^2)
    jH = ones(size(iH))
    sH = zeros(size(iH))
    k = 0
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (i1 - 1) * nely + j1
            for i2 = max(i1 - (ceil(rmin) - 1), 1):min(i1 + (ceil(rmin) - 1), nelx)
                for j2 = max(j1 - (ceil(rmin) - 1), 1):min(j1 + (ceil(rmin) - 1), nely)
                    e2 = (i2 - 1) * nely + j2
                    k += 1
                    iH[k] = e1
                    jH[k] = e2
                    sH[k] = max(0, rmin - sqrt((i1 - i2)^2 + (j1 - j2)^2))
                end
            end
        end
    end
    H = sparse(iH, jH, sH)
    Hs = [sum(H[i, :]) for i = 1:(size(H)[1])]

    return H, Hs
end

"""
Optimality criteria update
"""
function OC(
    nelx::Int64,
    nely::Int64,
    x,
    volfrac::Float64,
    dc::Matrix{Float64},
    dv::Matrix{Float64},
    H,
    Hs,
    filter::U,
) where U <: Filter
    l1 = 0
    l2 = 1e9
    move = 0.2
    xnew = zeros(nely, nelx)

    while (l2 - l1) / (l1 + l2) > 1e-3
        lmid = 0.5 * (l2 + l1)
        RacBe = sqrt.(-dc ./ dv / lmid)
        XB = x .* RacBe

        for i = 1:nelx
            for j = 1:nely
                xji = x[j, i]
                xnew[j, i] = max(0.000, max(xji - move, min(1, min(xji + move, XB[j, i]))))
            end
        end
        
        xPhys = OC_update_filter(filter, xnew, H, Hs)

        if sum(xPhys[:]) > volfrac * nelx * nely
            l1 = lmid
        else
            l2 = lmid
        end
    end

    return xnew
end

# NOTE -- the below function signatures are gross. Can this be avoided.
function OC_update_filter(s::SensitivityFilter, xnew, H, Hs)
    return xnew
end

function OC_update_filter(d::DensityFilter, xnew, H, Hs)
    return (H * xnew[:]) ./ Hs
end


end
