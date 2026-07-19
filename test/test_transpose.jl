# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

isdefined(@__MODULE__, :gen_test_waveforms) || include("test_utils.jl")

using LinearAlgebra
using JLArrays

using RadiationDetectorDSP: _lazy_transpose, _nonlazy_transpose, _row_major, _col_major, CPUNormAdaptor

import GPUArraysCore

GPUArraysCore.allowscalar(false)


@testset "transpose" begin
    A = rand(Float32, 7, 5)

    @testset "_nonlazy_transpose" begin
        At = _nonlazy_transpose(A)
        @test At isa Matrix
        @test At == transpose(A)

        A_gpu = JLArray(A)
        At_gpu = with_allowscalar(() -> _nonlazy_transpose(A_gpu))
        @test At_gpu isa JLArray
        @test size(At_gpu) == (5, 7)
        @test Array(At_gpu) == transpose(A)
        # Input must not be mutated:
        @test Array(A_gpu) == A

        # Tile-size edge cases of the KA transpose kernel:
        for sz in ((16, 16), (17, 33), (1, 40))
            B = rand(Float32, sz...)
            @test Array(with_allowscalar(() -> _nonlazy_transpose(JLArray(B)))) == transpose(B)
        end
    end

    @testset "_row_major / _col_major" begin
        @test _row_major(A) isa Transpose
        @test _row_major(A) == A
        @test _row_major(transpose(A)) === transpose(A)

        @test _col_major(A) === A
        @test _col_major(transpose(A)) isa Matrix
        @test _col_major(transpose(A)) == transpose(A)
    end

    @testset "CPUNormAdaptor" begin
        @test adapt(CPUNormAdaptor(), JLArray(A)) isa Array
        @test adapt(CPUNormAdaptor(), JLArray(A)) == A

        At = with_allowscalar(() -> adapt(CPUNormAdaptor(), transpose(JLArray(A))))
        @test At == transpose(A)
    end
end
