# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test


@testset "Filters" begin
    plot(args...) = nothing
    plot!(args...) = nothing
    hline!(args...) = nothing

    @testset "SimpleCSAFilter" begin
        t_drift = 500
        current_signal = vcat(fill(0.0, 10), fill(1.0 / t_drift, t_drift), fill(0.0, 200))
        flt = SimpleCSAFilter(tau_rise = 20, tau_decay = 20)
        plot(cumsum(current_signal))
        output = flt(current_signal)
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
