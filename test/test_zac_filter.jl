# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using Unitful
using RadiationDetectorSignals, ArraysOfArrays, FillArrays
using InverseFunctions
using Adapt

using RadiationDetectorDSP: bc_rdfilt, bc_rdfilt!


@testset "ZACChargeFilter" begin
    step_signal = RCFilter(10)(vcat(fill(0.0, 300), fill(1.0, 300)))
    step_wf = RDWaveform(15u"ns"*(eachindex(step_signal)), step_signal)

    flt = ZACChargeFilter(tau = 15u"ns"*10, toplen = 15u"ns"*20, length = 15u"ns"*100)

    @test flt(step_wf) isa RDWaveform
end
