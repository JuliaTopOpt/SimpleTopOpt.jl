module Structs

include("structs/parameters.jl")
include("structs/optimizers.jl")
include("structs/domains.jl")
include("structs/boundaryconditions.jl")
include("structs/femdefinitions.jl")
include("structs/problems.jl")
include("structs/solutions.jl")

@reexport using .Parameters
@reexport using .Optimizers
@reexport using .Domains
@reexport using .BoundaryConditions
@reexport using .FEMDefinitions
@reexport using .Problems
@reexport using .Solutions

end