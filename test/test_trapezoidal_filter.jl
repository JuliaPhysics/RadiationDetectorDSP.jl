# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using Unitful
using RadiationDetectorSignals, ArraysOfArrays, FillArrays
using InverseFunctions
using Adapt
using DSP

using RadiationDetectorDSP: bc_rdfilt, bc_rdfilt!


@testset "ConvolutionFilter" begin
    flt = TrapezoidalChargeFilter(6,4)
    uflt = TrapezoidalChargeFilter(6 * 1.5u"ns", 4 * 1.5u"ns")
    
    @test adapt(Array, flt) isa TrapezoidalChargeFilter

    noise = Float32[0.010744692, 0.011901839, 0.05194631, -0.011829155, -0.0061071543, -0.00019944388, 0.018830482, -0.061652776, 0.017121963, 0.040664993, 0.010738007, -0.06324256, 0.11602211, -0.011405745, 0.0582322, -0.0011262955, 0.10291999, -0.033106547, 0.07659402, -0.014946247, 0.076609485, 0.0036761805, -0.00330089, 0.037618753, 0.003061756, -0.13259014, -0.04753995, 0.06222269, -0.0016896317, -0.09834757, 0.029373525, 0.03482066, -0.10341741, -0.019071044, 0.05447747, 0.016070966, 0.1408137, 0.012808924, 0.0141442865, -0.041680623, -0.027131805, 0.046791397, -0.03010413, 0.11516451, 0.10181626, -0.010448066, -0.024512103]

    wfs_x = ArrayOfRDWaveforms((
        Fill(1.5u"ns" .* (1:47), 10),
        ArrayOfSimilarArrays([vcat(fill(0.5f0, 14 - i), zeros(Float32, 9 + i),ones(Float32, 24)) .+ noise for i in 1:10])
    ))
    wf_x = wfs_x[1]
    x = wf_x.signal
    
    ref_flt = ConvolutionFilter(flt)
    y_ref = ref_flt(x)
    wf_y_ref = RDWaveform(wf_x.time[16:end], y_ref)

    X = wfs_x.signal
    Y_ref = ref_flt.(X)
    wfs_y_ref = ArrayOfRDWaveforms((map(t -> t[16:end], wfs_x.time), Y_ref))

    fi = @inferred fltinstance(flt, smplinfo(x))
    y = similar(x, length(x) - 15)
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

    Y = nestedview(similar(flatview(X), (innersize(X, 1) - 15, size(X, 1))))
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
