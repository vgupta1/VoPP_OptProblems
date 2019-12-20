###  
# Optimization problems for VoPP Upper and Lower Bounds
###
#Throughout we assume mu = 1, c = 0, and hbar s.t. E[hbar(Vc)] = 0.
#Methods with aleading underscore are not true bounds but approximations for debugging/configuration

module vopp

using JuMP, Gurobi, Ipopt, LinearAlgebra, QuadGK, Roots

include("LambertW.jl")
include("Helpers.jl")
include("closedFormUB.jl")
#include("debugging.jl")
include("upperBounds.jl")
include("UpperBoundsSymmetric.jl")
include("UpperBoundsUnimodal.jl")
include("lowerBoundsSymmetric.jl")
include("lowerBoundsUnimodal.jl")


end #ends module
