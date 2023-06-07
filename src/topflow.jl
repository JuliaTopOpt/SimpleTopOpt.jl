module TopFlow

using LinearAlgebra
using SparseArrays


export topflow


### Constant parameters
# Physical parameters
const Uin = 1e0
const rho = 1e0
const mu = 1e0

# Continuation strategy
const conit = 50

# Optimization parameters
const mvlim = 0.2
const chlim = 1e-5
const chnum = 5

# Newton Solver parameters
const nltol = 1e-6
const nlmax = 25


"""
    `topflow`

Driving function for the TopFlow module
"""
function topflow(
    problemtype::S = 1,
    volfrac::T = 0.3333,
    Lx::T = 1.0,
    Ly::T = 1.0,
    nely::S = 30,
) where {T<:AbstractFloat,S<:Integer}

    @assert Lx > 0.0
    @assert Ly > 0.0
    @assert nely > 0

    if (problemtype > 2 || problemtype < 1)
        throw(ArgumentError("Problem type argument `problemtype` must be 1 or 2"))
    end

    # Assemble ProblemParam
    PP = ProblemParam(probtype, volfrac, Lx, Ly)
    xinit = PP.volfrac

    # Brinkman and Continuation
    alphamax = 2.5 * mu / (0.01^2)
    alphamin = 2.5 * mu / (100^2)
    ainit = 2.5 * mu / (0.1^2)
    qinit =
        (-xinit * (alphamax - alphamin) - ainit + alphamax) / (xinit * (ainit - alphamin))
    qavec = qinit ./ [1 2 10 20]
    qanum = length(qavec)
    conit = 50

    # Assemble SolvOptParam
    SOP = SolvOptParam(qanum, conit, mvlim)

    ### PREPARE FINITE ELEMENT ANALYSIS
    # Create FEAContainer
    FC = FEAContainer(nely, Lx, Ly)

    ### DEFINE BOUDNARY CONDITIOSN
    PC =
        (problemtype == 1) ? Problem1Container(nely, nodx, nody, nodtot) :
        Problem2Container(nely, nelx, nody, noxx, nodtot)
    # Nullspace matrices for imposing boundary conditions
    EN = sparse(I, doftot, doftot)
    ND = copy(EN)
    ND[fixedDofs, fixedDofs] = 0.0
    EN -= ND
    # Vectors for free dofs
    alldofs = 1:doftot
    freedofs = setdiff(alldofs, fixedDofs)

    ### INITIALIZATION
    # Solution vector
    S = zeros((doftot, 1))
    dS = copy(S)
    L = copy(S)
    S[fixedDofs] = DIR[fixedDofs]
    # Design field
    xPhys = xinit * ones(nely, nelx)
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
    qa = qavec[1]
    # Vectorized constants
    dxv = FC.dx * ones(1, neltot)
    dyv = FC.dy * ones(1, neltot)
    muv = mu * ones(1, neltot)
    rhov = rho * ones(1, neltot)

    # Print out problem information
    println("=========================================================")
    println(
        "      Problem number:" *
        String(problemtype) *
        " - Reynolds number: " *
        String(PC.Renum),
    )
    println("=========================================================")

    # Start iteration
    # destime = tic; ittime = tic 
    # TODO: ^^ TIMING -- just do benchmarking instead?
    while loop <= SOP.maxiter
        # TODO
        # if plotdes 
        # 
        # end

        # Greyscale indicator
        Md = 100 * full(4 * sum(xPhys .* (1 - xPhys)) / FC.neltot)

        # Material interpolation
        alpha = alphamin + (alphamax - alphamin) * (1 .- xPhys) ./ (1 .+ qa * xPhys)
        dalpha =
            (qa * (alphamax - alphamin) * (xPhys .- 1)) ./ (xPhys * qa .+ 1) .^ 2 -
            (alphamax - alphamin) ./ (xPhys * qa .+ 1)
    end

    return 0
end


"""
Struct containing components required for problem specification
"""
struct ProblemParam{T<:AbstractFloat,S<:Integer}
    # Problem type (1 = Double Pipe, 2 = Pipe Bend)
    probtype::S
    # Volume fraction
    volfrac::T
    # Domains length in x and y respectively.
    Lx::T
    Ly::T
end


struct SolvOptParam{S<:Integer}
    maxiter::S # maxiter = qanum * conit
    plotdes::Bool
    plotres::Bool
end

"""
Struct containing all components required for the finite element analysis
"""
mutable struct FEAContainer{T<:AbstractFloat,S<:Integer}
    dx::T
    dy::T

    nodx::S
    nody::S
    nodtot::S

    neltot::S
    doftot::S

    nodenrs::Matrix{S}
    edofMat::Matrix{S}

    iJ::Matrix{S}
    jJ::Matrix{S}
    iR::Matrix{S}
    jR::Matrix{S}
    jE::Matrix{S}

    """
    Constructor
    """
    function FEAContainer(nely::S = 30, Lx::T, Ly::T) where {T<:AbstractFloat,S<:Integer}

        @assert nely > 0.0
        @assert Lx > 0.0
        @assert Ly > 0.0

        nelx = nely * Lx / Ly

        dx = Lx / nelx
        dy = Ly / nely
        nodx = nelx + 1
        nody = nely + 1
        nodtot = nodx * nody
        neltot = nelx * nely
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

        new(
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



mutable struct Problem1Container{S,T} where {S<:AbstractMatrix,T<:AbstractFloat}
    """
    Specification of Problem 1 (double pipe problem)
    """

    fixedDofs::S
    Renum::T
    DIR::S

    function Problem1Container(nely::S, nodx::S, nody::S, nodtot::S) where {S<:Integer}
        """
        Constructor
        """

        if mod(nely, 6) != 0
            throw(
                ArgumentError(
                    "ERROR: Number of elements in y-dir. must be divisable by 6.",
                ),
            )
        elseif (nely < 30)
            throw(ArgumentError("ERROR: Number of elements in y-dir. must be above 30."))
        end

        inletLength = 1 / 6 * nely
        inlet1 = 1 / 6 * nely + 1
        inlet2 = 4 / 6 * nely + 1
        outletLength = 1 / 6 * nely
        outlet1 = 1 / 6 * nely + 1
        outlet2 = 4 / 6 * nely + 1

        nodesInlet = [(inlet1:(inlet1+inletLength))' (inlet2:(inlet2+inletLength))']
        nodesOutlet =
            (nodx - 1) * nody +
            [(outlet1:(outlet1:outletLength))' (outlet2:(outlet2+outletLength))']
        nodesTopBot = [1:nely+1:nodtot nody:nody:nodtot]
        nodesLefRig = [2:nely (nelx)*nody+2:nodtot-1]

        fixedDofsTBx = 2 * nodesTopBot - 1
        fixedDofsTBy = 2 * nodesTopBotw
        fixedDofsLRx = 2 * setdiff(nodesLefRig, [nodesInlet nodesOutlet]) - 1
        fixedDofsLRy = 2 * nodesLefRig
        fixedDofsInX = 2 * nodesInlet - 1
        fixedDofsInY = 2 * nodesInlet
        fixedDofsOutY = 2 * nodesOutlet
        fixedDofsOutP = 2 * nodtot + nodesOutlet
        fixedDofsU = [
            fixedDofsTBx fixedDofsTBy fixedDofsLRx fixedDofsLRy
            fixedDofsInx fixedDofsInY fixedDofsOutY
        ]
        fixedDofsP = [fixedDofsOutP]
        fixedDofs = [fixedDofsU fixedDofsP]

        ## Dirichlet vectors
        Uinlet = Uin * [-4 * (x^2 - x) for x in ([0:inletLength]' / inletLength)]
        DIRU = zeros(nodtot * 2, 1)
        DIRU[fixedDofsInX] = [Uinlet' Uinlet']
        DIRP = zeros(nodtot, 1)
        DIR = [DIRU; DIRP]

        # Inlet Reynolds Number
        Renum = Uin * (inletLength * Ly / nely) * rho / mu

        new(fixedDofs, Renum, DIR)
    end
end


mutable struct Problem2Container{S,T} where {S<:AbstractMatrix,T<:AbstractFloat}
    """
    Specification of Problem 2 (pipe bend problem)
    """

    fixedDofs::S
    Renum::T
    DIR::S

    function Problem2Container(
        nely::S,
        nelx::S,
        nody::S,
        nodx::S,
        nodtot::S,
    ) where {S<:Integer}
        """
        Constructor
        """

        throw(ArgumentError("TO BE DONE"))
    end
end


function NLNS(TEMP)
    """
    Non-linear Newton Solver
    """

    normR = 1
    nlit = 0
    fail = -1
    # TODO timing
    # nltime = tic

end






end
