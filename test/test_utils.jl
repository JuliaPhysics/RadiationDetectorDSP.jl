# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using Unitful
using RadiationDetectorSignals, ArraysOfArrays, FillArrays
using InverseFunctions
using Adapt

using RadiationDetectorDSP: bc_rdfilt, bc_rdfilt!


function gen_test_waveforms()
    noise = Float32[0.010744692, 0.011901839, 0.05194631, -0.011829155, -0.0061071543, -0.00019944388, 0.018830482, -0.061652776, 0.017121963, 0.040664993, 0.010738007, -0.06324256, 0.11602211, -0.011405745, 0.0582322, -0.0011262955, 0.10291999, -0.033106547, 0.07659402, -0.014946247, 0.076609485, 0.0036761805, -0.00330089, 0.037618753, 0.003061756, -0.13259014, -0.04753995, 0.06222269, -0.0016896317, -0.09834757, 0.029373525, 0.03482066, -0.10341741, -0.019071044, 0.05447747, 0.016070966, 0.1408137, 0.012808924, 0.0141442865, -0.041680623, -0.027131805, 0.046791397, -0.03010413, 0.11516451, 0.10181626, -0.010448066, -0.024512103]

    ArrayOfRDWaveforms((
        Fill(1.5u"ns" .* (-3:43), 10),
        ArrayOfSimilarArrays([vcat(fill(0.5f0, 14 - i), zeros(Float32, 9 + i),ones(Float32, 24)) .+ noise for i in 1:10])
    ))
end


function test_filter(flt, uflt, wfs_x, wfs_y_ref)
    @test @inferred(adapt(Array, flt)) isa getfield(parentmodule(typeof(flt)), nameof(typeof(flt)))
    @test @inferred(adapt(Array, uflt)) isa getfield(parentmodule(typeof(uflt)), nameof(typeof(uflt)))

    wf_x = wfs_x[1]
    x = wf_x.signal

    wf_y_ref = wfs_y_ref[1]
    y_ref = wf_y_ref.signal

    X = wfs_x.signal
    Y_ref = wfs_y_ref.signal

    fi = @inferred fltinstance(flt, smplinfo(x))
    y = similar(y_ref)
    fill!(y, NaN)

    @test @inferred(rdfilt!(y, fi, x)) === y
    @test y ≈ y_ref
    #fill!(y, NaN)
    #@test @inferred(rdfilt!(y, flt, x)) === y
    #@test y ≈ y_ref


    @test @inferred(rdfilt(fi, x)) ≈ y_ref
    #@test @inferred(rdfilt(flt, x)) ≈ y_ref

    #@test @inferred(fi(x)) ≈ y_ref
    @test @inferred(flt(x)) ≈ y_ref

    @test @inferred(rdfilt(fi, wf_x)) ≈ wf_y_ref
    #@test @inferred(rdfilt(flt, wf_x)) ≈ wf_y_ref

    #@test @inferred(fi(wf_x)) ≈ wf_y_ref
    @test @inferred(uflt(wf_x)) ≈ wf_y_ref

    Y = similar(Y_ref)
    fill!.(Y, NaN)

    @test @inferred(bc_rdfilt!(Y, fi, X)) === Y
    @test Y ≈ Y_ref
    #@test @inferred(bc_rdfilt!(Y, flt, X)) === Y
    #@test Y ≈ Y_ref

    @test @inferred(bc_rdfilt(fi, X)) ≈ Y_ref
    #@test @inferred(bc_rdfilt(flt, X)) ≈ Y_ref
    #@test @inferred(broadcast(fi, X)) ≈ Y_ref
    @test @inferred(broadcast(flt, X)) ≈ Y_ref
    #@test typeof(broadcast(fi, X)) == typeof(Y_ref)
    @test typeof(broadcast(flt, X)) == typeof(Y_ref)

    @test @inferred(bc_rdfilt(fi, wfs_x)) ≈ wfs_y_ref
    #@test @inferred(bc_rdfilt(flt, wfs_x)) ≈ wfs_y_ref
    #@test @inferred(broadcast(fi, wfs_x)) ≈ wfs_y_ref
    @test @inferred(broadcast(uflt, wfs_x)) ≈ wfs_y_ref
    #@test typeof(broadcast(fi, wfs_x)) == typeof(wfs_y_ref)
    @test typeof(broadcast(uflt, wfs_x)) == typeof(wfs_y_ref)
end
