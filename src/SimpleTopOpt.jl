module SimpleTopOpt

using Reexport

include("structs/parameters.jl")
include("structs/domains.jl")
include("structs/optimizers.jl")
include("structs/boundaryconditions.jl")
include("structs/femdefinitions.jl")
include("structs/problems.jl")

include("top88.jl")
include("toph.jl")
include("topflow.jl")

@reexport using .Structs
@reexport using .Top88
@reexport using .TopH
@reexport using .TopFlow

end