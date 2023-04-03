
module Top88

using LinearAlgebra
using SparseArrays
using Statistics

export top88
export prepare_filter
export OC

"""
    top88(nelx, nely, volfrac, penal, rmin, ft)

A direct, naive Julia port of Andreassen et al. "Efficient topology optimization in MATLAB
using 88 lines of code." By default, this will reproduce the optimized MBB beam from Sigmund
(2001).

# Arguments
- `nelx::S`: Number of elements in the horizontal direction
- `nely::S`: Number of elements in the vertical direction
- `volfrac::T`: Prescribed volume fraction
- `penal::T`: The penalization power
- `rmin::T`: Filter radius divided by the element size
- `ft::Bool`: Choose between sensitivity (if true) or density filter (if false)
- `write::Bool`: If true, will write out iteration number, changes, and density

# Returns
- `x`: Final material distribution, presented as a matrix

"""
function top88(
    nelx::S=60,
    nely::S=20,
    volfrac::T=0.5,
    penal::T=3.0,
    rmin::T=2.0,
    ft::Bool=true,
    write::Bool=false,
    loop_max::Int=1000
) where {S <: Integer, T <: AbstractFloat}
    # Physical parameters
    E0 = 1; Emin = 1e-9; nu = 0.3;

    # Prepare finite element analysis
    A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12]
    A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6]
    B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4]
    B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2]
    KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11])

    nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx)
    edofVec = reshape(2*nodenrs[1:end-1,1:end-1].+1, nelx*nely, 1)
    edofMat = zeros(Int64, nelx*nely, 8)

    offsets = [0 1 2*nely.+[2 3 0 1] -2 -1]
    for i = 1:8
        for j = 1:nelx*nely
            edofMat[j,i]= edofVec[j] + offsets[i]
        end
    end

    iK = reshape(kron(edofMat,ones(8,1))', 64*nelx*nely,1)
    jK = reshape(kron(edofMat,ones(1,8))', 64*nelx*nely,1)
    
    # Loads and supports
    # OLD: F = sparse([2], [1], [-1], 2*(nely+1)*(nelx+1), 1)
    F = spzeros(2*(nely+1)*(nelx+1))
    F[2,1] = -1
    U = spzeros(2*(nely+1)*(nelx+1))

    fixeddofs = union(1:2:2*(nely+1), [2*(nelx+1)*(nely+1)])
    alldofs = 1:2*(nely+1)*(nelx+1)
    freedofs = setdiff(alldofs, fixeddofs)

    # Prepare the filter
    H, Hs = prepare_filter(nelx, nely, rmin)
    
    # Initialize iteration
    x = volfrac*ones(nely,nelx)
    xPhys = x
    loop = 0
    change = 1
    cValues = []

    # Start iteration
    while change > 0.01
        loop += 1
        # FE-Analysis
        sK = [j*((i+Emin)^penal) for i in ((E0-Emin)*xPhys[:]') for j in KE[:]]
        K = sparse(iK[:], jK[:], sK)
        K = (K+K')/2

        KK = cholesky(K[freedofs,freedofs])
        U[freedofs] = KK \ F[freedofs]
        
        # OLD: edM = [convert(Int64,i) for i in edofMat]
        mat = (U[edofMat]*KE).*U[edofMat]

        # Objective function and sensitivity analysis
        ce = reshape([sum(mat[i,:]) for i = 1:(size(mat)[1])],nely,nelx)
        c = sum(sum((Emin*ones(size(xPhys)).+(xPhys.^penal)*(E0-Emin)).*ce))
        push!(cValues,c)
        dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce
        dv = ones(nely,nelx)
        
        # Filtering/ modification of sensitivities
        if ft
            dc[:] = H*(x[:].*dc[:])./Hs./max(1e-3,maximum(x[:]))
        else
            dc[:] = H*(dc[:]./Hs)
            dv[:] = H*(dv[:]./Hs)
        end

        # Optimality criteria update of design variables and physical densities
        # TODO -- probable issue with xPhys not pass by reference?
        xnew = OC(nelx, nely, x, volfrac, dc, dv, xPhys, ft)

        change = maximum(abs.(x-xnew))
        x = xnew

        write && println("Loop = ", loop, ", Change = ", change ,", c = ", c, ", structural density = ", mean(x))
        loop >= loop_max && break       
    end

    return x
end


"""
Prepare filter
"""
function prepare_filter(nelx::S, nely::S, rmin::T) where {S <: Integer, T <: AbstractFloat}
    iH = ones(nelx*nely*(2*(convert(Int64,ceil(rmin)-1))+1)^2)
    jH = ones(size(iH))
    sH = zeros(size(iH))
    k = 0
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (i1-1)*nely+j1
            for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                    e2 = (i2-1)*nely+j2
                    k += 1
                    iH[k] = e1
                    jH[k] = e2
                    sH[k] = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2))
                end
            end
        end
    end
    H = sparse(iH,jH,sH)
    Hs = [sum(H[i,:]) for i = 1:(size(H)[1])]

    return H, Hs
end

"""
Optimality criteria update
"""
function OC(nelx, nely, x, volfrac, dc, dv, xPhys, ft)
    l1 = 0; l2 = 1e9; move = 0.2
    xnew = zeros(nely, nelx)

    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1)
        RacBe = sqrt.(-dc./dv/lmid)
        XB = x.*RacBe

        for i = 1:nelx
            for j = 1:nely
                xji = x[j,i]
                xnew[j,i]= max(0.000,max(xji-move,min(1,min(xji+move,XB[j,i]))))
            end
        end  

        if ft
            xPhys = xnew
        else
            xPhys[:] = (H*xnew[:])./Hs
        end

        if sum(xPhys[:]) > volfrac*nelx*nely
            l1 = lmid
        else 
            l2 = lmid 
        end
    end

    return xnew
end

end 
