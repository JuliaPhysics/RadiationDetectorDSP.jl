# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using Random


@testset "legacy charge_trapflt!" begin
    navg, ngap = 5, 3
    n = 100
    x = vcat(fill(0.0, 30), fill(1.0, n - 30)) .+ 0.01 .* randn(Xoshiro(1), n)

    y_legacy = charge_trapflt!(copy(x), navg, ngap)
    @test length(y_legacy) == n

    # Away from the padded tail, the legacy in-place filter matches the
    # current TrapezoidalChargeFilter:
    y_ref = TrapezoidalChargeFilter(navg, ngap)(x)
    @test y_legacy[1:(n - 2 * navg - ngap)] ≈ y_ref[1:(n - 2 * navg - ngap)]

    @test_throws ArgumentError charge_trapflt!(copy(x), 0, 3)
    @test_throws ArgumentError charge_trapflt!(copy(x), 5, -1)
    @test_throws ArgumentError charge_trapflt!(rand(10), 5, 3)
end


@testset "legacy pulse generators" begin
    smpls = zeros(Float64, 20)
    @test add_rect_pulse!(smpls, 5, 10, 2.0) === smpls
    @test smpls == vcat(zeros(4), fill(2.0, 10), zeros(6))

    pulse = gen_rect_pulse(20, 5, 10, 2.0)
    @test pulse == smpls
end
