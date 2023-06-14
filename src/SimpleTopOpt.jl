module SimpleTopOpt

include("structs.jl")
include("top88.jl")
include("toph.jl")
include("topflow.jl")

using .Top88
using .TopH
using .TopFlow

export top88
export topH

export OCParameters

export BrinkmanPenalizationParamters, SIMPParameters

export TopflowDomain, DoublePipeContainer, PipeBendContainer


end
