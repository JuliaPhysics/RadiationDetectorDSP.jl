# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using Statistics


@testset "circuit_filters" begin
    plot(args...) = nothing
    plot!(args...) = nothing
    hline!(args...) = nothing

    t_drift = 40
    current_signal = vcat(fill(0.0, 10), fill(1.0 / t_drift, t_drift), fill(0.0, 200))

    step_signal = vcat(fill(0.0, 10), fill(1.0, 30))

    @testset "RCFilter" begin
        x = step_signal
        plot(x)
        flt = RCFilter(rc = 20)
        output = flt(x)
        plot!(output)
        plot!(inverse(flt)(output))
        hline!([1 - exp(-1)])
        @test inverse(flt) isa InvRCFilter
        @test inverse(inverse(flt)) == flt
        InverseFunctions.test_inverse(flt, x)
    end

    @testset "CRFilter" begin
        x = vcat(fill(0.0, 10), fill(1.0, 30))
        plot(x)
        flt = CRFilter(cr = 10)
        output = flt(x)
        plot!(output)
        plot!(inverse(flt)(output))
        hline!([exp(-1)])
        @test inverse(flt) isa InvCRFilter
        @test inverse(inverse(flt)) == flt
        InverseFunctions.test_inverse(flt, x)
    end

    @testset "ModCRFilter" begin
        x = vcat(fill(0.0, 10), fill(1.0, 30))
        plot(x)
        flt = ModCRFilter(cr = 10)
        output = flt(x)
        plot!(output)
        plot!(inverse(flt)(output))
        hline!([exp(-1)])
        @test inverse(flt) isa InvModCRFilter
        @test inverse(inverse(flt)) == flt
        InverseFunctions.test_inverse(flt, x)
    end

    @testset "IntegratorFilter" begin
        x = current_signal
        plot(x)
        flt = IntegratorFilter(gain = 2.0)
        output = flt(x)
        plot!(output)
        plot!(inverse(flt)(output))
        @test inverse(flt) isa DifferentiatorFilter
        @test inverse(inverse(flt)) == flt
        InverseFunctions.test_inverse(flt, x)
    end

    @testset "SimpleCSAFilter" begin
        x = current_signal
        plot(cumsum(x))
        flt = SimpleCSAFilter(tau_rise = 20, tau_decay = 500)
        output = flt(x)
        plot!(output)
        output_deconv = inverse(CRFilter(τ_decay))(output)
        plot!(output_deconv)
        tail = output_deconv[150:end]
        # Tail of reco should be flat:
        @test var(tail) < 1e-5
        InverseFunctions.test_inverse(flt, x)
    end
end
