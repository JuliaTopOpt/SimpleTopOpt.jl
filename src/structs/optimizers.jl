"""
This module defines the optimization techniques available for use
and structs to maintain relevant parameters.
"""

abstract type Optimizer end

"""@docs
    OptimalityCriteria

Defines an Optimizer for the Optimality Criteria method
"""
@kwdef struct OptimalityCriteria <: Optimizer
    max_iter::Int32 = 200
    mvlim::Float64 = 0.2
end