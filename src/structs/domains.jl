module Domains

abstract type AbstractDomain end

##################################################################################################
# Top88 domains
##################################################################################################

abstract type AbstractTop88Domain <: AbstractDomain end

"""@docs
    Top88Domain

Domain parameters for Top88 problems
"""
@kwdef struct Top88Domain <: AbstractTop88Domain
    nelx::Int64 = 60
    nely::Int64 = 30
end

##################################################################################################
# Toph domains
##################################################################################################

abstract type AbstractTophDomain <: AbstractDomain end

"""@docs
    TophDomain

Domain parameters for TopH problems
"""
@kwdef struct TophDomain <: AbstractTophDomain
    nelx::Int64 = 30
    nely::Int64 = 30
end

##################################################################################################
# Topflow domains
##################################################################################################

abstract type AbstractTopflowDomain <: AbstractDomain end

"""@docs
    TopflowDomain

Domain parameters for Topflow problems
"""
@kwdef struct TopflowDomain <: AbstractTopflowDomain
    Lx::Float64 = 1.0
    Ly::Float64 = 1.0

    dx::Float64 = (1 / 30)
    dy::Float64 = (1 / 30)

    nely::Int64 = 30
    nelx::Int64 = 30

    """@docs
        TopflowDomain

    Constructor for defining a domain for Topflow problems
    """
    function TopflowDomain(Lx::Float64, Ly::Float64, nely::Int64 = 30)
        @assert Lx > 0.0
        @assert Ly > 0.0
        @assert nely > 1

        nelx = nely * Lx / Ly
        dx = Lx / nelx
        dy = Ly / nely

        new(Lx, Ly, dx, dy, nely, nelx)
    end
end


end