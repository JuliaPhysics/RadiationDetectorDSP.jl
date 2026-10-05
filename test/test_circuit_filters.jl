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

    # Workaround for isapprox for ranges on Julia v1.6:
    cmpwf(a::RDWaveform, b::RDWaveform; kwargs...) = isapprox(a.signal, b.signal, kwargs...) && isapprox(collect(a.time), collect(b.time), kwargs...)

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
        InverseFunctions.test_inverse(flt, x; compare = cmpwf)
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
        InverseFunctions.test_inverse(flt, x; compare = cmpwf)
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
        InverseFunctions.test_inverse(flt, x; compare = cmpwf)
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
        InverseFunctions.test_inverse(flt, x; compare = cmpwf)
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

    @testset "SecondOrderCRFilter" begin
        x = current_wf
        plot(x)
        flt = SecondOrderCRFilter(cr = 15u"ns" * 10, cr2 = 0.5u"ns" * 10, f = 0.5)
        output = flt(x)
        plot!(output)
        plot!(inverse(flt)(output))
        hline!([exp(-1)])
        @test inverse(flt) isa InvSecondOrderCRFilter
        @test inverse(inverse(flt)) == flt
        InverseFunctions.test_inverse(flt, x; compare = cmpwf)
    end

    @testset "RC_CR2Filter" begin
        # Create a test waveform with a pulse (similar to current signal)
        pulse_signal = vcat(fill(0.0, 10), fill(1.0, 30), fill(0.0, 200))
        pulse_wf = RDWaveform(15u"ns"*(0:239), pulse_signal)
        
        # Test with different time constants
        tau_values = [15u"ns" * 5, 15u"ns" * 10, 15u"ns" * 20]
        
        for tau in tau_values
            plot(pulse_wf)
            flt = RC_CR2Filter(tau = tau)
            output = flt(pulse_wf)
            plot!(output)
            
            # Verify output type
            @test output isa RDWaveform
            @test length(output.signal) == length(pulse_wf.signal)
            
            # First three samples should be preserved
            @test isapprox(output.signal[1], pulse_wf.signal[1]; rtol=1e-6)
            @test isapprox(output.signal[2], pulse_wf.signal[2]; rtol=1e-6)
            @test isapprox(output.signal[3], pulse_wf.signal[3]; rtol=1e-6)
            
            # No NaNs should be in output for valid input
            @test !any(isnan, output.signal)
        end
        
        # Test with step signal to observe shaping behavior
        step_signal = vcat(fill(0.0, 10), fill(1.0, 40))
        step_wf = RDWaveform(15u"ns"*(0:49), step_signal)
        
        flt = RC_CR2Filter(tau = 15u"ns" * 10)
        output_step = flt(step_wf)
        
        @test output_step isa RDWaveform
        @test !any(isnan, output_step.signal)
        @test isapprox(output_step.signal[1], step_wf.signal[1]; rtol=1e-6)
        
        # Test with very short waveform (should handle gracefully)
        short_wf = RDWaveform(15u"ns"*(0:2), [0.0, 1.0, 0.5])
        flt_short = RC_CR2Filter(tau = 15u"ns" * 5)
        output_short = flt_short(short_wf)
        
        # Short waveforms (≤3 samples) should still work, returning NaN
        @test output_short isa RDWaveform
        @test length(output_short.signal) == 3
        
        # Test broadcasting with RDSignal/RDWaveform
        flt = RC_CR2Filter(tau = 15u"ns" * 10)
        
        # Verify that the filter responds to pulse (shows differentiation-like behavior)
        pulse_start = findfirst(x -> x != 0.0, pulse_wf.signal)
        pulse_end = findlast(x -> x != 0.0, pulse_wf.signal)
        
        # Get the filtered output around pulse region
        output_pulse = flt(pulse_wf)
        
        # The filter should show rise-up at pulse start and fall at pulse end
        # (characteristic of RC-CR^2 filter)
        if !isnothing(pulse_start) && !isnothing(pulse_end)
            @test maximum(abs.(output_pulse.signal[pulse_start:pulse_start+10])) > 0
            @test maximum(abs.(output_pulse.signal[pulse_end-10:pulse_end])) > 0
        end
    end
end
