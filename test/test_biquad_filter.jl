# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using Random, LinearAlgebra, InverseFunctions
using Unitful
using RadiationDetectorSignals, ArraysOfArrays, FillArrays


@testset "circuit_filters" begin
    flt = BiquadFilter((0.2, 0.15, 0.3), (-0.8, 0.4))

    wfs = ArrayOfRDWaveforms((
        Fill(1.5u"ns" .* (1:47), 10),
        ArrayOfSimilarArrays([vcat(fill(0.5f0, 14 - i), zeros(Float32, 9 + i),ones(Float32, 24)) for i in 1:10])
    ))
    wf_x = wfs[1]
    x = wf_x.signal
    # From DSP.filt(DSP.Biquad(flt.b_012..., flt.a_12...), x):
    y_ref = [0.1, 0.255, 0.489, 0.6142, 0.62076, 0.575928, 0.5374384, 0.52457952, 0.529688256, 0.5389187968, 0.54425973504, 0.544840269312, 0.5431683214336, 0.44159854942208004, 0.286011510964224, 0.05216978900254721, -0.07266877318365184, -0.07900293414794035, -0.03413483804489155, 0.0042933032232628995, 0.01708857779656694, 0.011953540947948394, 0.0027274016397319393, 0.1974005049326062, 0.5068294432901922, 0.9765033526591114, 1.2284709048112121, 1.242175382785325, 1.152351944303775, 1.0750114023288901, 1.0490683441416022, 1.0592501143817257, 1.0777727538487396, 1.0885181573263014, 1.0897054243215452, 1.0863570765267156, 1.0832034914927544, 1.0820199625835174, 1.0823345734697123, 1.0830596737423628, 1.0835139096060054, 1.0835872581878592, 1.083464242707885, 1.0833364908911645, 1.0832834956297777, 1.0832922001473564, 1.0833203618659741]
    wf_y_ref = RDWaveform(wf_x.time, y_ref)

    fi = @inferred fltinstance(flt, smplinfo(x))
    y = similar(x)

    @test @inferred(rdfilt!(y, fi, x)) ≈ y_ref
    @test y ≈ y_ref
    fill!(y, 0)
    @test @inferred(rdfilt!(y, flt, x)) ≈ y_ref
    @test y ≈ y_ref

    @test @inferred(rdfilt(fi, x)) ≈ y_ref
    @test @inferred(rdfilt(flt, x)) ≈ y_ref

    @test @inferred(fi(x)) ≈ y_ref
    @test @inferred(flt(x)) ≈ y_ref

    @test @inferred(rdfilt(fi, wf_x)) ≈ wf_y_ref
    @test @inferred(rdfilt(flt, x)) ≈ y_ref

    @test isapprox(@inferred(inverse(flt)(y)), x, rtol = 1e-3)
end
