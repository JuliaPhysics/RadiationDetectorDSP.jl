# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using InverseFunctions
using RadiationDetectorSignals, Unitful
using Statistics


@testset "circuit_filters" begin
    plot(args...) = nothing
    plot!(args...) = nothing
    hline!(args...) = nothing

    t_drift = 40
    current_signal = vcat(fill(0.0, 10), fill(1.0 / t_drift, t_drift), fill(0.0, 200))
    current_wf = RDWaveform(15u"ns"*(0:249), current_signal)

    step_signal = vcat(fill(0.0, 10), fill(1.0, 30))
    step_wf = RDWaveform(15u"ns"*(0:39), step_signal)

    @testset "RCFilter" begin
        x = current_wf
        plot(x)
        flt = RCFilter(rc = 20 * 15u"ns")
        output = flt(x)
        plot!(output)
        plot!(inverse(flt)(output))
        hline!([1 - exp(-1)])
        @test inverse(flt) isa InvRCFilter
        @test inverse(inverse(flt)) == flt
        InverseFunctions.test_inverse(flt, x)
    end

    @testset "CRFilter" begin
        x = step_wf
        plot(x)
        flt = CRFilter(cr = 15u"ns" * 10)
        output = flt(x)
        plot!(output)
        plot!(inverse(flt)(output))
        hline!([exp(-1)])
        @test inverse(flt) isa InvCRFilter
        @test inverse(inverse(flt)) == flt
        InverseFunctions.test_inverse(flt, x)
    end

    @testset "ModCRFilter" begin
        x = step_wf
        plot(x)
        flt = ModCRFilter(cr = 15u"ns" * 10)
        output = flt(x)
        plot!(output)
        plot!(inverse(flt)(output))
        hline!([exp(-1)])
        @test inverse(flt) isa InvModCRFilter
        @test inverse(inverse(flt)) == flt
        InverseFunctions.test_inverse(flt, x)
    end

    @testset "IntegratorFilter" begin
        x = current_wf
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
        x = current_wf
        plot(RDWaveform(x.time, cumsum(x.signal)))
        flt = SimpleCSAFilter(tau_rise = 15u"ns" * 20, tau_decay = 15u"ns" * 500)
        output = flt(x)
        plot!(output)
        output_deconv = inverse(CRFilter(cr = 15u"ns" * 500))(output)
        plot!(output_deconv)
        tail = output_deconv.signal[150:end]
        # Tail of reco should be flat:
        @test var(tail) < 1e-5
    end
end
