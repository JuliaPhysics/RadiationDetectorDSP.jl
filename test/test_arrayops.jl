# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

isdefined(@__MODULE__, :gen_test_waveforms) || include("test_utils.jl")

using Random

using RadiationDetectorDSP: bc_rdfilt_arrayops, rand_waveforms


@testset "rand_waveforms" begin
    n_samples, n_waveforms = 256, 20
    r = rand_waveforms(Xoshiro(3), Float32, n_samples, n_waveforms; noise_sigma = 0)
    X = r.waveforms
    @test X isa Matrix{Float32}
    @test size(X) == (n_samples, n_waveforms)
    @test length(r.amplitudes) == length(r.tau_rises) == length(r.tau_falls) == n_waveforms
    # Flat baseline before the trigger, pulse height below the collected charge:
    @test X[1:(r.i_trigger - 1), :] ≈ repeat(r.baseline_levels', r.i_trigger - 1)
    pulse_heights = vec(maximum(X, dims = 1)) .- r.baseline_levels
    @test all(0.4 .* r.amplitudes .< pulse_heights .< r.amplitudes)
end


@testset "bc_rdfilt_arrayops" begin
    X = rand_waveforms(Xoshiro(4), Float32, 300, 12).waveforms
    signals = ArrayOfSimilarVectors(X)

    filters = [
        TrapezoidalChargeFilter(40, 10),
        TrapezoidalChargeFilter(40, 10, 25),
        InvCRFilter(30f0),
        CRFilter(30f0),
        RCFilter(5f0),
        BiquadFilter((0.2f0, 0.15f0, 0.3f0), (-0.8f0, 0.4f0)),
        BiquadFilter((0.5f0, 0.1f0, 0f0), (-1.5f0, 0.56f0)),
        ConvolutionFilter(DirectConvolution(), Float32[0.25, 0.5, 0.15, 0.1]),
        ConvolutionFilter(FFTConvolution(), Float32[0.25, 0.5, 0.15, 0.1]),
        SavitzkyGolayFilter(7, 3, 0),
    ]

    for flt in filters
        Y_ref = flatview(flt.(signals))
        Y = flatview(@inferred(bc_rdfilt_arrayops(flt, signals)))
        @test size(Y) == size(Y_ref)
        @test Y ≈ Y_ref rtol = 1e-4
    end

    # Integer samples
    signals_int = ArrayOfSimilarVectors(round.(Int32, X))
    for flt in (TrapezoidalChargeFilter(40, 10), CRFilter(30f0))
        @test flatview(bc_rdfilt_arrayops(flt, signals_int)) ≈ flatview(flt.(signals_int)) rtol = 1e-4
    end
end
