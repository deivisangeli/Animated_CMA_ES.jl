__precompile__()


module Animated_CMA_ES

using Plots, LinearAlgebra, StaticArrays


export cma_es, cma_es_animated


include("CMA_ES.jl")
include("cma_es_animated.jl")

end # module
