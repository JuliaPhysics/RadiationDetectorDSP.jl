# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).
using Test

include("test_utils.jl")

@testset "signalstats" begin
    # Flat signal: mean=5, sigma=0, slope=0
    xs = 0.0:1.0:9.0
    ys = fill(5.0, 10)
    wf = RDWaveform(xs, ys)

    stats = signalstats(wf, 0.0, 9.0)
    @test stats.mean ≈ 5.0
    @test stats.sigma ≈ 0.0 atol=1e-10
    @test stats.slope ≈ 0.0 atol=1e-10
    @test stats.offset ≈ 5.0

    # Linear signal y = 2x + 1: slope=2, offset=1, sigma small
    ys_linear = 2.0 .* collect(xs) .+ 1.0
    wf_linear = RDWaveform(xs, ys_linear)

    stats_linear = signalstats(wf_linear, 0.0, 9.0)
    @test stats_linear.slope ≈ 2.0 atol=1e-10
    @test stats_linear.offset ≈ 1.0 atol=1e-10
    @test stats_linear.slope_residual_sigma ≈ 0.0 atol=1e-10

    # Subrange: only use x in [2.0, 5.0]
    stats_sub = signalstats(wf_linear, 2.0, 5.0)
    @test stats_sub.slope ≈ 2.0 atol=1e-10
    @test stats_sub.mean ≈ 2.0 * 3.5 + 1.0 atol=1e-10  # mean x in [2,5] = 3.5
end