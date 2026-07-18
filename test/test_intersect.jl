# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using Unitful
using RadiationDetectorSignals


@testset "Intersect" begin
    step_signal = vcat(fill(0.0, 10), fill(2.0, 10))

    @testset "plain samples" begin
        r = @inferred Intersect(mintot = 2)(step_signal, 1.0)
        @test r.multiplicity == 1
        @test r.x ≈ 10.5

        # Intersection exactly at threshold crossing between two samples:
        r2 = Intersect(mintot = 1)(vcat(fill(0.0, 5), fill(1.0, 5)), 0.5)
        @test r2.x ≈ 5.5
    end

    @testset "no intersection" begin
        r = @inferred Intersect(mintot = 2)(fill(0.0, 20), 1.0)
        @test r.multiplicity == 0
    end

    @testset "min time-over-threshold" begin
        # Single-sample spike must not count for mintot = 3:
        spiky = vcat(fill(0.0, 5), [2.0], fill(0.0, 8), fill(2.0, 6))
        r = Intersect(mintot = 3)(spiky, 1.0)
        @test r.multiplicity == 1
        @test 14 < r.x < 15

        # With mintot = 1 the spike counts, too:
        r2 = Intersect(mintot = 1)(spiky, 1.0)
        @test r2.multiplicity == 2
        @test 5 < r2.x < 6
    end

    @testset "waveform with units" begin
        wf = RDWaveform((0:19) * 1.0u"ns", step_signal)
        r = @inferred Intersect(mintot = 2u"ns")(wf, 1.0)
        @test r.multiplicity == 1
        @test r.x ≈ 9.5u"ns"
    end
end
