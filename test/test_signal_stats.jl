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

    # Large mean compared to the spread, on a Float32 time axis far from zero
    ts = range(60000f0u"ns", step = 16f0u"ns", length = 4000)
    t_ns = ustrip.(u"ns", collect(ts))
    ys_big = 2e4 .- 0.01 .* t_ns .+ 170 .* sin.(1.7 .* (1:length(ts)) .^ 2)
    tb, yb = big.(t_ns), big.(ys_big)
    mt, my = sum(tb) / length(tb), sum(yb) / length(yb)
    var_t, var_y = sum((tb .- mt) .^ 2) / length(tb), sum((yb .- my) .^ 2) / length(yb)
    cov_ty = sum((tb .- mt) .* (yb .- my)) / length(tb)
    ref_slope = cov_ty / var_t
    ref = (mean = my, sigma = sqrt(var_y), slope = ref_slope, offset = my - ref_slope * mt,
           slope_residual_sigma = sqrt(var_y - cov_ty^2 / var_t))
    for (T, rtol) in ((Float64, 1e-10), (Float32, 1e-4))
        stats_big = signalstats(RDWaveform(ts, T.(ys_big)), first(ts), last(ts))
        @test stats_big.mean ≈ ref.mean rtol = rtol
        @test stats_big.sigma ≈ ref.sigma rtol = rtol
        @test ustrip(u"ns^-1", stats_big.slope) ≈ ref.slope rtol = rtol
        @test stats_big.offset ≈ ref.offset rtol = rtol
        @test stats_big.slope_residual_sigma ≈ ref.slope_residual_sigma rtol = rtol
    end
end
