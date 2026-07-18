# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

isdefined(@__MODULE__, :gen_test_waveforms) || include("test_utils.jl")


@testset "global_ops" begin
    @testset "shift_waveform" begin
        wfs_x = gen_test_waveforms()
        A = randn(length(wfs_x))
        @test @inferred(shift_waveform(wfs_x[1], 0.2)) ≈ RDWaveform(wfs_x[1].time, wfs_x[1].signal .+ 0.2)
        @test @inferred(shift_waveform(wfs_x[1].signal, 0.2)) ≈ wfs_x[1].signal .+ 0.2
        @test @inferred(broadcast(shift_waveform, wfs_x, A)) ≈ ArrayOfRDWaveforms(map((a,b) -> shift_waveform(a,b), wfs_x, A))
        @test @inferred(broadcast(shift_waveform, wfs_x, A[1])) ≈ ArrayOfRDWaveforms(map((a,b) -> shift_waveform(a,b), wfs_x, fill(A[1], length(wfs_x))))
        @test broadcast(shift_waveform, wfs_x, A).time isa Fill
        @test broadcast(shift_waveform, wfs_x, A).signal isa ArrayOfSimilarVectors
        @test @inferred(broadcast(shift_waveform, wfs_x.signal, A)) ≈ map((a,b) -> a .+ b, wfs_x.signal, A)
        @test broadcast(shift_waveform, wfs_x.signal, A) isa ArrayOfSimilarVectors
        @test all(broadcast(shift_waveform, [wfs_x[1], wfs_x[2]], 0.2) .≈ [shift_waveform(wfs_x[1], 0.2), shift_waveform(wfs_x[2], 0.2)])
    end

    @testset "multiply_waveform" begin
        wfs_x = gen_test_waveforms()
        A = randn(length(wfs_x))
        @test @inferred(multiply_waveform(wfs_x[1], 2.5)) ≈ RDWaveform(wfs_x[1].time, wfs_x[1].signal .* 2.5)
        @test @inferred(multiply_waveform(wfs_x[1].signal, 2.5)) ≈ wfs_x[1].signal .* 2.5
        @test @inferred(broadcast(multiply_waveform, wfs_x, A)) ≈ ArrayOfRDWaveforms(map((a,b) -> multiply_waveform(a,b), wfs_x, A))
        @test broadcast(multiply_waveform, wfs_x, A).time isa Fill
        @test broadcast(multiply_waveform, wfs_x, A).signal isa ArrayOfSimilarVectors
        @test @inferred(broadcast(multiply_waveform, wfs_x.signal, A)) ≈ map((a,b) -> a .* b, wfs_x.signal, A)
        @test broadcast(multiply_waveform, wfs_x.signal, A) isa ArrayOfSimilarVectors
    end

    @testset "reverse_waveform" begin
        wfs_x = gen_test_waveforms()
        @test @inferred(reverse_waveform(wfs_x[1])) ≈ RDWaveform(wfs_x[1].time, reverse(wfs_x[1].signal))
        @test @inferred(reverse_waveform(wfs_x[1].signal)) ≈ reverse(wfs_x[1].signal)
        @test @inferred(broadcast(reverse_waveform, wfs_x)) ≈ ArrayOfRDWaveforms(map(reverse_waveform, wfs_x))
        @test broadcast(reverse_waveform, wfs_x).time isa Fill
        @test broadcast(reverse_waveform, wfs_x).signal isa ArrayOfSimilarVectors
        @test @inferred(broadcast(reverse_waveform, wfs_x.signal)) ≈ map(reverse, wfs_x.signal)
        @test broadcast(reverse_waveform, wfs_x.signal) isa ArrayOfSimilarVectors
    end
end
