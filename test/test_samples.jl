# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using Unitful
using RadiationDetectorSignals, ArraysOfArrays, FillArrays

using RadiationDetectorDSP: dspfloattype, _floattype, _smpltype, _smpllen


@testset "samples" begin
    @testset "dspfloattype" begin
        @test dspfloattype(Float64) == Float64
        @test dspfloattype(Float32) == Float32
        @test dspfloattype(Int16) == Float32
        @test dspfloattype(UInt16) == Float32
        @test dspfloattype(Int32) == Float32
        @test dspfloattype(Int64) == Float64

        for T in (Float64, Float32, Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64)
            @test _floattype(T) == dspfloattype(T)
        end
    end

    @testset "SamplingInfo" begin
        si = @inferred SamplingInfo{Float32}(0:99)
        @test _smpltype(si) == Float32
        @test _smpllen(si) == 100
    end

    @testset "smplinfo" begin
        smpls = rand(Float32, 20)
        si = @inferred smplinfo(smpls)
        @test _smpltype(si) == Float32
        @test si.axis == eachindex(smpls)

        wf = RDWaveform((0:19) * 1.5u"ns", smpls)
        si_wf = @inferred smplinfo(wf)
        @test _smpltype(si_wf) == Float32
        @test si_wf.axis == wf.time
    end

    @testset "elsmplinfo" begin
        smpls_vec = [rand(Float32, 20) for _ in 1:5]
        si = elsmplinfo(smpls_vec)
        @test _smpltype(si) == Float32
        @test _smpllen(si) == 20

        si_aosa = elsmplinfo(ArrayOfSimilarVectors(rand(Float32, 20, 5)))
        @test _smpltype(si_aosa) == Float32
        @test _smpllen(si_aosa) == 20

        @test_throws DimensionMismatch elsmplinfo([rand(Float32, 20), rand(Float32, 21)])

        t = (0:19) * 1.5u"ns"
        wfs = ArrayOfRDWaveforms((Fill(t, 5), ArrayOfSimilarVectors(rand(Float32, 20, 5))))
        si_wfs = elsmplinfo(wfs)
        @test _smpltype(si_wfs) == Float32
        @test si_wfs.axis == t

        wfs_plainvec = [RDWaveform(t, rand(Float32, 20)) for _ in 1:5]
        si_pv = elsmplinfo(wfs_plainvec)
        @test si_pv.axis == t

        t_other = (0:19) * 2.0u"ns"
        @test_throws ArgumentError elsmplinfo([RDWaveform(t, rand(Float32, 20)), RDWaveform(t_other, rand(Float32, 20))])
    end
end
