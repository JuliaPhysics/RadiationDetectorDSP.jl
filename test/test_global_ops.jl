# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

include("test_utils.jl")


@testset "global_ops" begin
    @testset "shift_waveform" begin
        wfs_x = gen_test_waveforms()
        A = randn(length(wfs_x))
        @test @inferred(shift_waveform(wfs_x[1], 0.2)) ≈ RDWaveform(wfs_x[1].time, wfs_x[1].signal .+ 0.2)
        @test @inferred(broadcast(shift_waveform, wfs_x, A)) ≈ ArrayOfRDWaveforms(map((a,b) -> shift_waveform(a,b), wfs_x, A))
        @test @inferred(broadcast(shift_waveform, wfs_x, A[1])) ≈ ArrayOfRDWaveforms(map((a,b) -> shift_waveform(a,b), wfs_x, fill(A[1], length(wfs_x))))
        @test broadcast(shift_waveform, wfs_x, A).time isa Fill
        @test broadcast(shift_waveform, wfs_x, A).signal isa ArrayOfSimilarVectors 
        # ToDo: Add more detailed tests for shift_waveform
    end
end
