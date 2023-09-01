module ParameterDefinitions
"""
This module defines structs maintaining input parameter values
"""

export Parameters, Filter
export ModifiedSIMPParameters, DensityFilter, SensitivityFilter, BrinkmanPenalizationParameters
export TopflowContinuation, TopflowPhysics, TopflowNumericals

abstract type Parameters end

##################################################################################################
# Material interpolation techniques
##################################################################################################

"""@docs
    ModifiedSIMPParameters

Modified SIMP (solid isotropic material power) law parameter container. The fields
`E0` and `Emin` are required for Top88 problems, though default values are given;
these are not used for TopFlow.
"""
@kwdef struct ModifiedSIMPParameters <: Parameters
    penal::Float64 = 3.0
    
    E0::Float64 = 1.0
    Emin::Float64 = 1e-9
end

##################################################################################################
# Topflow specific
##################################################################################################

abstract type Filter <: Parameters end

"""
Implements the sensitivity filter for Toph and Top88 problems
"""
@kwdef struct SensitivityFilter <: Filter
    rmin::Float64 = 1.2
end


"""
Implements the density filter for Toph and Top88 problems
"""
@kwdef struct DensityFilter <: Filter
    rmin::Float64 = 1.2
end

##################################################################################################
# Topflow specific
##################################################################################################

"""@docs
    BrinkmanPenalizationParameters

Defines the Brinkman penalization technique for Topflow problems
"""
@kwdef struct BrinkmanPenalizationParameters <: Parameters
    alphamax::Float64 = 25000
    alphamin::Float64 = 0.00025

    function BrinkmanPenalizationParameters(mu::Float64)
        new(2.5 * mu / (0.01^2), 2.5 * mu / (100^2))
    end
end

"""@docs
    TopflowContinuation

Defines the TopFlow continuation strategy and maintains the relevant parameters.
"""
struct TopflowContinuation <: Parameters
    ainit::Float64
    qinit::Float64
    qavec::Matrix{Float64}
    qanum::Int64
    conit::Int64

    bkman::BrinkmanPenalizationParameters

    function TopflowContinuation(
        bkman::BrinkmanPenalizationParameters,
        volfrac::Float64 = (1/3),
        conit::Int64 = 50,
    )
        ainit = bkman.alphamax / 100

        qinit =
            (-volfrac * (bkman.alphamax - bkman.alphamin) - ainit + bkman.alphamax) /
            (volfrac * (ainit - bkman.alphamin))
        qavec = qinit ./ [1 2 10 20]
        qanum = length(qavec)

        new(ainit, qinit, qavec, qanum, conit, bkman)
    end
end

"""@docs
Topflow optimization and Newton solver parameters
"""
@kwdef struct TopflowNumericals <: Parameters
    # Main optimization loop 
    maxiter::Int64 = 200
    mvlim::Float64 = 0.2
    chlim::Float64 = 1e-3
    chnum::Int64 = 5

    # Newton Solver Parameters
    nltol::Float64 = 1e-6
    nlmax::Int64 = 25
end

"""
Physical parameters for Topflow
"""
@kwdef mutable struct TopflowPhysics <: Parameters
    Uin::Float64 = 1e0
    rho::Float64 = 1e0
    mu::Float64 = 1e0

    Renum::Float64 = undef
end


end