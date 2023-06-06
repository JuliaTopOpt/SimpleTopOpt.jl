module TopFlow

using LinearAlgebra
using SparseArrays


export topflow


"""
  `topflow`

Driving function for the TopFlow module
"""
function topflow(
  problemtype::S=1,
  volfrac::T=0.3333,
  Lx::T=1.0,
  Ly::T=1.0,
  nely::S=30,
  Uin::T=1e0,
  rho::T=1e0,
  mu::T=1e0,
) where {T<:AbstractFloat, S<:Integer}

  @assert Lx > 0.0
  @assert Ly > 0.0
  @assert nely > 0

  if (problemtype > 2 || problemtype < 1)
    throw(ArgumentError("Problem type argument `problemtype` must be 1 or 2"))
  end

  # Assemble ProblemParam

  # Brinkman and Continuation

  # Assemble SolvOptParam

  ### PREPARE FINITE ELEMENT ANALYSIS
  # Create FEAContainer

  fea_struct = FEAContainer()

  # Define Boundary Conditions

  ### Initialization



  return 0
end


"""
Struct containing components required for problem specification
"""
mutable struct ProblemParam{T<:AbstractFloat, S<:Integer}
    # Problem type (1 = Double Pipe, 2 = Pipe Bend)
    probtype::S
    # Volume fraction
    volfrac::T
    
    ### Domain
    # Size
    Lx::T
    Ly::T
    
    ### Physical parameters
    # Properties
    Uin::T
    rho::T
    mu::T


    function ProblemParam{T, S}(
        problemtype::S=1,
        volfrac::T=0.3333,
        Lx::T=1.0,
        Ly::T=1.0,
        Uin::T=1e0,
        rho::T=1e0,
        mu::T=1e0,
    ) where {T<:AbstractFloat, S<:Integer}

        new(problemtype, volfrac, Lx, Ly, Uin, rho, mu)
    end
end

"""
Contains everything for the optimizer/ solver overall
"""
mutable struct SolverParam{T<:AbstractFloat, S<:Integer}
    mvlim::T


    function SolverParam{T, S}(
        mu::T,
    ) where {T<:AbstractFloat, S<:Integer}

    end
end





"""

Struct containing all components required for the finite element analysis

"""
mutable struct FEAContainer{T<:AbstractFloat, S<:Integer}

  # TODO: fill out the fields
  dx::T
  dy::T

  nodx::S
  nody::S
  nodtot::S

  neltot::S
  doftot::S

  nodenrs::Matrix{S}
  edofVecU::Matrix{S}
  edofMatU::Matrix{S}
  edofVecP::Matrix{S}
  edofMatP::Matrix{S}
  edofMat::Matrix{S}


  """
  Constructor
  """
  function FEAContainer(
    nely::S=30,
    Lx::T,
    Ly::T
  ) where {T<:AbstractFloat, S<:Integer}

    @assert nely > 0.0
    @assert Lx > 0.0
    @assert Ly > 0.0

    nelx = nely * Lx / Ly

    dx = Lx/nelx; dy = Ly/nely
    nodx = nelx+1; nody = nely+1; nodtot = nodx*nody
    neltot = nelx*nely; doftot = 3 *nodtot

    nodenrs = reshape(1:nodtot, nody, nodx)
    edofVecU = reshape(2 * nodenrs[1:end-1, 1:end-1]+1, neltot, 1)
    edofMatU = repeat(edofVecU,1,8)+repeat( [0 1 2*nely+[2 3 0 1] -2 -1], neltot, 1)
    edofVecP = reshape(nodenrs[1:end-1,1:end-1], neltot, 1)
    edofMatP = repeat(edofVecP,1,4)+repeat([1 nely+[2 1] 0], neltot, 1)
    edofMat = [edofMatU 2*nodtot+edofMatP]
    
    # iJ = reshape(kron(edofMat, ones(12, 1))', 144*neltot, 1)
    # jJ = reshape(kron(edofMat, ones(1, 12))', 144*neltot, 1)
    iR = reshape(edofMat', 12*neltot, 1); jR = ones(12*neltot, 1)
    jE = repeat(1:neltot, 12, 1)

  end
end



mutable struct ProblemContainer
"""
Results of problems.m
"""

  function ProblemContainer()
    
  end
end





end
