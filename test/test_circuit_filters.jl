# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using Statistics


@testset "circuit_filters" begin
    plot(args...) = nothing
    plot!(args...) = nothing
    hline!(args...) = nothing

    @testset "SimpleCSAFilter" begin
        t_drift = 40
        current_signal = vcat(fill(0.0, 10), fill(1.0 / t_drift, t_drift), fill(0.0, 200))
        plot(cumsum(current_signal))
        flt = SimpleCSAFilter(tau_rise = 20, tau_decay = 500)
        output = flt(current_signal)
        plot!(output)
        output_deconv = inverse(CRFilter(τ_decay))(output)
        plot!(output_deconv)
        tail = output_deconv[150:end]
        # Tail of reco should be flat:
        @test var(tail) < 1e-5
        output_full_deconv = inverse(flt)(output)
        plot!(output_full_deconv)
        @test output_full_deconv ≈ current_signal
    end

    @testset "rc_filter" begin
        RC = 20
        x = vcat(fill(0.0, 10), fill(1.0, 40))
        plot(x)
        flt = RCFilter(RC)
        output = flt(x)
        # output = df1_filt(rc_filter(RC), X)
        plot!(output)
        output_deconv = inverse(flt)(output)
        plot!(output_deconv)
        hline!([1 - exp(-1)])
        @test x ≈ output_deconv
    end

    @testset "cr_filter" begin
        CR = 10
        x = vcat(fill(0.0, 10), fill(1.0, 30))
        plot(x)
        flt = CRFilter(CR)
        output = flt(x)
        # output = df1_filt(cr_filter(RC), X)
        plot!(output)
        output_deconv = inverse(flt)(output)
        # output_deconv = df1_filt(inv_cr_filter(RC), output)
        plot!(output_deconv)
        hline!([exp(-1)])
        @test x ≈ output_deconv
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
