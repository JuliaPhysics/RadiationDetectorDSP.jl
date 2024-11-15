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
using StructArrays
using TypedTables
using UnPack
using Unitful

import ChainRulesCore
import DSP
import GPUArraysCore
import SIMD

using Unitful: RealOrRealQuantity as RealQuantity


include("ka_compat.jl")
include("fast_indexing.jl")
include("computedevs.jl")
include("samples.jl")
include("math_utils.jl")
include("transpose.jl")
include("array_utils.jl")
include("legacy/filters.jl")
include("legacy/trapezoidal_filter.jl")
include("legacy/generators.jl")
include("filter.jl")
include("global_ops.jl")
include("signalstats.jl")
include("lsqfit.jl")
include("intersect.jl")
include("signal_estimator.jl")
include("truncate_filter.jl")
include("convolution_filter.jl")
include("biquad_filter.jl")
include("first_order_iir.jl")
include("circuit_filters.jl")
include("trapezoidal_filter.jl")
include("sg_filter.jl")
include("zac_filter.jl")
include("cusp_filter.jl")
include("gaussian_filter.jl")


@static if !isdefined(Base, :get_extension)
    using Requires
end

function __init__()
    @static if !isdefined(Base, :get_extension)
        @require CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba" include("../ext/RadiationDetectorDSPCUDAExt.jl")
        @require JLArrays = "27aeb0d3-9eb9-45fb-866b-73c2ecf80fcb" include("../ext/RadiationDetectorDSPJLArraysExt.jl")
    end
end

end # module
