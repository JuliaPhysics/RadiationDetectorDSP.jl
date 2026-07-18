# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using Random, Statistics
using Unitful
using RadiationDetectorSignals, ArraysOfArrays, FillArrays

include("testutils_rand_waveforms.jl")


@testset "DSP chain on rand detector waveforms" begin
    n_samples, n_waveforms = 512, 20
    tau_fall = 64.0

    r = rand_detector_waveforms(
        Xoshiro(42), Float64, n_samples, n_waveforms;
        tau_fall = (tau_fall, tau_fall)
    )
    wfs = r.waveforms

    @test wfs isa ArrayOfRDWaveforms
    @test size(wfs) == (n_waveforms,)
    @test all(length.(wfs.signal) .== n_samples)

    baseline = signalstats.(wfs, 5, r.i_trigger - 20)
    @test all(abs.(getproperty.(baseline, :mean) .- r.baseline_levels) .< 4 .* r.noise_levels)

    wfs_blsub = shift_waveform.(wfs, .-getproperty.(baseline, :mean))
    wfs_pz = InvCRFilter(tau_fall).(wfs_blsub)

    # Pole-zero corrected waveforms should have a flat tail:
    tail_stats = signalstats.(wfs_pz, n_samples - 100, n_samples - 1)
    @test all(abs.(getproperty.(tail_stats, :slope)) .< 0.02 .* r.noise_levels)

    # Energy reconstruction via trapezoidal filter:
    wfs_trap = TrapezoidalChargeFilter(100, 60).(wfs_pz)
    E = maximum.(map(wf -> wf.signal, wfs_trap))
    @test all(abs.(E .- r.amplitudes) ./ r.amplitudes .< 0.05)
    @test median(abs.(E .- r.amplitudes) ./ r.amplitudes) < 0.02

    # Trigger location via threshold intersect:
    trig = Intersect(mintot = 4).(wfs_pz, 0.5 .* r.amplitudes)
    @test all(getproperty.(trig, :multiplicity) .== 1)
    delay = getproperty.(trig, :x) .- r.i_trigger
    @test all(0 .< delay .< r.n_collect .+ 2 .* r.tau_rises)

    # Signal estimation on the shaped waveform (should be close to the
    # sample values in the flat-ish tail region):
    sigest = SignalEstimator(PolynomialDNI(degree = 3, length = 7))
    pos = n_samples - 150
    y_est = sigest.(wfs_pz, pos)
    y_ref = map(wf -> wf.signal[pos], wfs_pz)
    @test all(abs.(y_est .- y_ref) .< 4 .* r.noise_levels)
end


@testset "DSP chain on rand detector waveforms (unitful, Float32)" begin
    n_samples, n_waveforms = 512, 10
    tau_fall = 64.0
    dt = 16u"ns"

    r = rand_detector_waveforms(
        Xoshiro(7), Float32, n_samples, n_waveforms;
        tau_fall = (tau_fall, tau_fall),
        time_axis = (0:(n_samples - 1)) .* 16
    )
    wfs = ArrayOfRDWaveforms((Fill((0:(n_samples - 1)) .* dt, n_waveforms), r.waveforms.signal))

    baseline = signalstats.(wfs, 5 * dt, (r.i_trigger - 20) * dt)
    wfs_blsub = shift_waveform.(wfs, .-getproperty.(baseline, :mean))
    wfs_pz = InvCRFilter(tau_fall * dt).(wfs_blsub)
    wfs_trap = TrapezoidalChargeFilter(100 * dt, 60 * dt).(wfs_pz)

    E = maximum.(map(wf -> wf.signal, wfs_trap))
    @test all(abs.(E .- r.amplitudes) ./ r.amplitudes .< 0.05)

    trig = Intersect(mintot = 4 * dt).(wfs_pz, 0.5 .* r.amplitudes)
    # Noise may re-cross the threshold on small-amplitude waveforms:
    @test all(getproperty.(trig, :multiplicity) .>= 1)
    delay = getproperty.(trig, :x) .- (r.i_trigger * dt)
    @test all(0u"ns" .< delay .< (r.n_collect .+ 2 .* r.tau_rises) .* dt)
end
