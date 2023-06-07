module TopH

using LinearAlgebra
using SparseArrays
using Statistics

export toph
export OC
export check
export FE

"""
    toph(nelx, nely, volfrac, penal, rmin, write, loop_max)

A direct, naive Julia port of the `toph` code listing from "Topology Optimization"
by Martin Bendsøe and Ole Sigmund.

# Arguments
- `nelx::S`: Number of elements in the horizontal direction
- `nely::S`: Number of elements in the vertical direction
- `volfrac::T`: Prescribed volume fraction
- `penal::T`: The penalization power
- `rmin::T`: Filter radius divided by the element size
- `write::Bool`: If true, will write out iteration number, changes, and density
    for each iteration. Defaults to false.
- `loop_max::Int`: Explicitly set the maximum number of iterations. Defaults to 1000.

# Returns
- `Matrix{T}`: Final material distribution, represented as a matrix.
"""
function toph(
    nelx::S,
    nely::S,
    volfrac::T,
    penal::T,
    rmin::T,
    write::Bool=false, 
    loop_max::Int=100
) where {S <: Integer, T <: AbstractFloat}
    # Initialization
    x = volfrac * ones(nely,nelx)
    loop = 0
    change = 1.
    dc = zeros(nely,nelx)

    # Start iteration
    while change > 0.01
        loop += 1
        xold = x
        c = 0.

        # FE Analysis
        U = FE(nelx,nely,x,penal)

        KE = [ 2/3 -1/6 -1/3 -1/6
              -1/6  2/3 -1/6 -1/3
              -1/3 -1/6  2/3 -1/6
              -1/6 -1/3 -1/6  2/3 ]

        # Objective function/ sensitivity analysis
        for ely = 1:nely
            for elx = 1:nelx
                n1 = (nely+1)*(elx-1)+ely
                n2 = (nely+1)* elx   +ely
                Ue = U[[n1; n2; n2+1; n1+1]]

                c += (0.001+0.999*x[ely,elx]^penal)*Ue'*KE*Ue
                dc[ely,elx] = -0.999*penal*(x[ely,elx])^(penal-1)*Ue'*KE*Ue
            end
        end

        # Sensitivity filtering 
        dc = check(nelx,nely,rmin,x,dc)
        # Design update by optimality criteria method
        x  = OC(nelx,nely,x,volfrac,dc)

        # Print out results if desired
        if write
            change = maximum(abs.(x-xold))
            println("Change = ", change, " c = ", c)
        end

        loop >= loop_max && break
    end

    return x
end

"""
    OC(nelx, nely, x, volfrac, dc)

Optimality criteria update

# Arguments
- `nelx::S`: Number of elements in the horizontal direction
- `nely::S`: Number of elements in the vertical direction
- `x::Matrix{T}`: Current material distribution
- `volfrac::T`: Prescribed volume fraction
- `dc::Matrix{T}`: Sensitivity filter

# Returns
- `Matrix{T}`: Updated material distribution

"""
function OC(
    nelx::S,
    nely::S,
    x::Matrix{T},
    volfrac::T,
    dc::Matrix{T}
) where {S <: Integer, T <: AbstractFloat}
    l1 = 0; l2 = 100000; move = 0.2
    xnew = zeros(nely,nelx)

    while (l2-l1) > 1e-4
        lmid = 0.5*(l2+l1)
        RacBe = sqrt.(-dc/lmid)
        XB = x .* RacBe

        for i = 1:nelx
            for j = 1:nely
                xji = x[j,i]
                xnew[j,i]= max(0.001,max(xji-move,min(1,min(xji+move,XB[j,i]))))
            end
        end

        if (sum(sum(xnew)) - volfrac*nelx*nely) > 0
            l1 = lmid
        else
            l2 = lmid
        end
    end

    return xnew
end

"""
    check(nelx, nely, rmin, x, dc)

Mesh independency filter

# Arguments
- `nelx::S`: Number of elements in the horizontal direction
- `nely::S`: Number of elements in the vertical direction
- `rmin::T`: Sensitivity filter radius divided by element size
- `x::Matrix{T}`: Current material distribution
- `dc::Matrix{T}`: Compliance derivatives

# Returns
- `Matrix{T}`: Updated dc
"""
function check(nelx::S,
    nely::S,
    rmin::T,
    x::Matrix{T},
    dc::Matrix{T}
) where {S <: Integer, T <: AbstractFloat}
    dcn=zeros(nely,nelx)

    for i = 1:nelx
      for j = 1:nely
        sum=0.0

        for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
          for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
            l = Int64(l); k = Int64(k)
            fac = rmin-sqrt((i-k)^2+(j-l)^2)
            sum = sum+max(0,fac)
            dcn[j,i] += max(0,fac)*x[l,k]*dc[l,k]
          end
        end

        dcn[j,i] = dcn[j,i]/(x[j,i]*sum)
      end
    end

    return dcn
end

"""
    FE(nelx, nely, x, penal)

Finite element implementation

# Arguments
- `nelx::S`: Number of elements in the horizontal direction
- `nely::S`: Number of elements in the vertical direction
- `x::Matrix{T}`: Current material distribution
- `penal::T`: The penalization power

# Returns
- `Matrix{T}`: Differential equation solution U
"""
function FE(
    nelx::S,
    nely::S,
    x::Matrix{T},
    penal::T
) where {S <: Integer, T <: AbstractFloat}
    KE = [ 2/3 -1/6 -1/3 -1/6
          -1/6  2/3 -1/6 -1/3
          -1/3 -1/6  2/3 -1/6
          -1/6 -1/3 -1/6  2/3 ]

    K = spzeros((nelx+1)*(nely+1), (nelx+1)*(nely+1))
    U = zeros((nely+1)*(nelx+1))
    F = zeros((nely+1)*(nelx+1))
    for elx = 1:nelx
        for ely = 1:nely 
            n1 = (nely+1)*(elx-1)+ely
            n2 = (nely+1)* elx   +ely
            edof = [n1; n2; n2+1; n1+1]
            K[edof,edof] += (0.001+0.999*x[ely,elx]^penal)*KE
        end
    end

    F .= 0.01
    fixeddofs = Int64(nely/2+1-(nely/20)):Int64(nely/2+1+(nely/20))
    alldofs = 1:(nely+1)*(nelx+1)
    freedofs = setdiff(alldofs,fixeddofs)

    U[freedofs] = K[freedofs, freedofs] \ F[freedofs]
    U[fixeddofs] .= 0

    return U
end

end
