# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using Unitful
using RadiationDetectorSignals, ArraysOfArrays, FillArrays


@testset "SignalEstimator" begin
    poly(x) = 2 + 3 * x + 0.5 * x^2
    xs = 0.0:1.0:99.0
    ys = poly.(xs)
    wf = RDWaveform(xs, ys)

    sigest = SignalEstimator(PolynomialDNI(degree = 2, length = 7))

    @testset "exact polynomial recovery" begin
        for x in (42.3, 17.0, 63.7)
            @test @inferred(sigest(wf, x)) ≈ poly(x) rtol = 1e-6
            @test sigest(ys, x + 1) ≈ poly(x) rtol = 1e-6
        end
    end

    @testset "default polynomial degree" begin
        sigest_default = SignalEstimator(PolynomialDNI(length = 7))
        @test sigest_default(wf, 42.3) ≈ poly(42.3) rtol = 1e-6
    end

    @testset "out of range" begin
        @test isnan(sigest(wf, -10.0))
        @test isnan(sigest(wf, 150.0))
    end

    @testset "denoising" begin
        noise = 0.1 .* sin.(100 .* xs)
        wf_noisy = RDWaveform(xs, ys .+ noise)
        @test sigest(wf_noisy, 42.3) ≈ poly(42.3) rtol = 1e-2
    end

    @testset "broadcast over waveforms" begin
        n_wf = 5
        Y = collect(hcat((ys .+ i for i in 1:n_wf)...))
        wfs = ArrayOfRDWaveforms((Fill(xs, n_wf), ArrayOfSimilarVectors(Y)))

        y_bc = sigest.(wfs, 42.3)
        @test y_bc ≈ [poly(42.3) + i for i in 1:n_wf] rtol = 1e-6

        x_pos = [40.1, 41.2, 42.3, 43.4, 44.5]
        y_bc2 = sigest.(wfs, x_pos)
        @test y_bc2 ≈ [poly(x_pos[i]) + i for i in 1:n_wf] rtol = 1e-6

        y_bc3 = sigest.(Ref(wf), x_pos)
        @test y_bc3 ≈ poly.(x_pos) rtol = 1e-6
    end

    @testset "waveform with units" begin
        wfu = RDWaveform(xs * u"ns", ys)
        sigest_u = SignalEstimator(PolynomialDNI(degree = 2, length = 7u"ns"))
        @test sigest_u(wfu, 42.3u"ns") ≈ poly(42.3) rtol = 1e-6
    end
end
