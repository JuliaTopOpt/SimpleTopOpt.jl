module Structs

using Reexport

include("structs/parameters.jl")
include("structs/optimizers.jl")
include("structs/domains.jl")
include("structs/femdefinitions.jl")
include("structs/boundaryconditions.jl")
include("structs/problems.jl")
include("structs/solutions.jl")

@reexport using .ParameterDefinitions
@reexport using .Optimizers
@reexport using .Domains
@reexport using .FEMDefinitions
@reexport using .BoundaryConditionDefinitions
@reexport using .Problems
@reexport using .Solutions

end