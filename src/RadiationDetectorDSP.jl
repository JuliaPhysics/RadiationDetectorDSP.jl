# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

__precompile__(true)

module RadiationDetectorDSP

using ArgCheck
using ArraysOfArrays
using Distributions
using DSP
using ElasticArrays
using LinearAlgebra
using ParallelProcessingTools
using Parameters
using RadiationDetectorSignals
using Random
using RecipesBase
using SIMD
using StaticArrays
using Statistics
using StatsBase
using TypedTables

include("util.jl")
include("samples.jl")
include("legacy/filters.jl")
include("legacy/trapezoidal_filter.jl")
include("legacy/generators.jl")
include("zac.jl")

end # module
