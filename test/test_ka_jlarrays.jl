# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

isdefined(@__MODULE__, :gen_test_waveforms) || include("test_utils.jl")

using JLArrays

import GPUArraysCore
import KernelAbstractions

GPUArraysCore.allowscalar(false)


@testset "KernelAbstractions filtering with JLArrays" begin
    @test RadiationDetectorDSP._ka_get_backend(JLArray(rand(Float32, 3))) isa KernelAbstractions.CPU

    n_samples, n_waveforms = 40, 8
    X = rand(Float32, n_samples, n_waveforms)
    inputs_cpu = ArrayOfSimilarVectors(X)
    inputs_gpu = ArrayOfSimilarVectors(JLArray(X))
    si = SamplingInfo{Float32}(1:n_samples)

    filters = [
        BiquadFilter((0.2f0, 0.15f0, 0.3f0), (-0.8f0, 0.4f0)),
        FirstOrderIIR((0.2f0, 0.15f0), (-0.8f0,)),
        TrapezoidalChargeFilter(5, 3),
        ConvolutionFilter(DirectConvolution(), Float32[0.25, 0.5, 0.25]),
    ]

    for flt in filters
        fi = fltinstance(flt, si)
        Y_cpu = bc_rdfilt(fi, inputs_cpu)
        Y_gpu = with_allowscalar(() -> bc_rdfilt(fi, inputs_gpu))
        @test flatview(Y_gpu) isa JLArray
        @test Array(flatview(Y_gpu)) ≈ flatview(Y_cpu)
    end
end
