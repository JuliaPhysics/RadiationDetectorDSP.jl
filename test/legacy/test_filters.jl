# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using Random, LinearAlgebra, Statistics
using DSP


@testset "Filters" begin
    # Runs Biquad filter in Direct Form 1, for comparison (DSP.jl uses
    # Transposed Direct Form 2):
    df1_filt(f::Biquad, input::Vector{<:Real}) =
        filt([f.b0, f.b1, f.b2], [one(f.a1), f.a1, f.a2], input)

    plot(args...) = nothing
    plot!(args...) = nothing
    hline!(args...) = nothing

    @testset "test_csa_response" begin
        t_drift, τ_rise, τ_decay = 20, 20, 500
        current_signal = vcat(fill(0.0, 10), fill(1.0 / t_drift, t_drift), fill(0.0, 200))
        plot(cumsum(current_signal))
        output = filt(simple_csa_response_filter(τ_rise, τ_decay, 1), current_signal)
        plot!(output)
        output_deconv = filt(inv_cr_filter(τ_decay), output)
        plot!(output_deconv)
        tail = output_deconv[150:end]
        # Tail of reco should be flat:
        @test var(tail) < 1e-5
    end

    @testset "rc_filter" begin
        RC = 20
        X = vcat(fill(0.0, 10), fill(1.0, 40))
        plot(X)
        output = filt(rc_filter(RC), X)
        # output = df1_filt(rc_filter(RC), X)
        plot!(output)
        output_deconv = filt(inv_rc_filter(RC), output)
        plot!(output_deconv)
        hline!([1 - exp(-1)])
        @test X ≈ output_deconv
    end

    @testset "cr_filter" begin
        RC = 10
        X = vcat(fill(0.0, 10), fill(1.0, 30))
        plot(X)
        output = filt(cr_filter(RC), X)
        # output = df1_filt(cr_filter(RC), X)
        plot!(output)
        output_deconv = filt(inv_cr_filter(RC), output)
        # output_deconv = df1_filt(inv_cr_filter(RC), output)
        plot!(output_deconv)
        hline!([exp(-1)])
        @test X ≈ output_deconv
    end

    @testset "crmod_filter" begin
        RC = 10
        X = vcat(fill(0.0, 10), fill(1.0, 30))
        plot(X)
        output = filt(crmod_filter(RC), X)
        # output = df1_filt(crmod_filter(RC), X)
        plot!(output)
        output_deconv = filt(inv_crmod_filter(RC), output)
        # output_deconv = df1_filt(inv_crmod_filter(RC), output)
        plot!(output_deconv)
        hline!([exp(-1)])
        @test X ≈ output_deconv
    end
end
