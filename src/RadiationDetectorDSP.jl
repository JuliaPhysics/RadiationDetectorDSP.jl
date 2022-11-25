# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

__precompile__(true)

module RadiationDetectorDSP

using LinearAlgebra
using Random
using Statistics

using Adapt
using ArgCheck
using ArraysOfArrays
using CompositionsBase
using DocStringExtensions
using ElasticArrays
using FFTW
using FillArrays
using FunctionChains
using IntervalSets
using KernelAbstractions
using InverseFunctions
using Parameters
using RadiationDetectorSignals
using StatsBase
using TypedTables
using UnPack
using Unitful

import ChainRulesCore
import DSP
import SIMD


include("array_utils.jl")
include("samples.jl")
include("legacy/filters.jl")
include("legacy/trapezoidal_filter.jl")
include("legacy/generators.jl")
include("filter.jl")
include("functional.jl")
include("global_ops.jl")
include("truncate_filter.jl")
include("convolution_filter.jl")
include("biquad_filter.jl")
include("first_order_iir.jl")
include("circuit_filters.jl")
include("trapezoidal_filter.jl")
include("sg_filter.jl")
include("zac_filter.jl")
include("fixed_pickoff.jl")

end # module
