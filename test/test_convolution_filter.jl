# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using Unitful
using RadiationDetectorSignals, ArraysOfArrays, FillArrays
using InverseFunctions
using DSP

using RadiationDetectorDSP: bc_rdfilt, bc_rdfilt!


@testset "BiquadFilter" begin
    # Savitzky-Golay filter -2:4 (length 7, evaluate at 3rd position), 3rd order, 1st deriv:
#    flt = ConvolutionFilter(DirectConvolution(), [-0.115079365079365, -0.18253968253968259, -0.0753968253968254, 0.09523809523809519, 0.21825396825396814, 0.18253968253968253, -0.12301587301587294])
    flt = ConvolutionFilter(FFTConvolution(), [-0.115079365079365, -0.18253968253968259, -0.0753968253968254, 0.09523809523809519, 0.21825396825396814, 0.18253968253968253, -0.12301587301587294])
    
    noise = Float32[0.010744692, 0.011901839, 0.05194631, -0.011829155, -0.0061071543, -0.00019944388, 0.018830482, -0.061652776, 0.017121963, 0.040664993, 0.010738007, -0.06324256, 0.11602211, -0.011405745, 0.0582322, -0.0011262955, 0.10291999, -0.033106547, 0.07659402, -0.014946247, 0.076609485, 0.0036761805, -0.00330089, 0.037618753, 0.003061756, -0.13259014, -0.04753995, 0.06222269, -0.0016896317, -0.09834757, 0.029373525, 0.03482066, -0.10341741, -0.019071044, 0.05447747, 0.016070966, 0.1408137, 0.012808924, 0.0141442865, -0.041680623, -0.027131805, 0.046791397, -0.03010413, 0.11516451, 0.10181626, -0.010448066, -0.024512103]

    wfs_x = ArrayOfRDWaveforms((
        Fill(1.5u"ns" .* (1:47), 10),
        ArrayOfSimilarArrays([vcat(fill(0.5f0, 14 - i), zeros(Float32, 9 + i),ones(Float32, 24)) .+ noise for i in 1:10])
    ))
    wf_x = wfs_x[1]
    x = wf_x.signal
    # From DSP.filt(DSP.FIRFilter(flt.coeffs), x):
    y_ref = Float32[-0.058776174, -0.15214051, -0.19546856, -0.1468839, -0.027337808, 0.07304379, 0.009391554, 0.008527434, -0.0020374788, -0.0010664854, -0.010996632, -0.006111667, -0.008577821, 0.06304957, 0.13707972, 0.17107229, 0.13423079, 0.046884775, -0.07579133, 0.011823545, 7.777956f-5, 0.005714143, -0.009477851, -0.096994996, -0.29966956, -0.34484184, -0.2542256, -0.040562287, 0.110750146, -0.030599957, -0.01432999, 0.019098327, 0.010801469, -0.008182973, -0.0030080022, 0.004854781, -0.044890568, -0.050525106, 0.008078488, 0.030457651, 0.03784929, 0.03058027, -0.018897993, -0.02196053, -0.04128293, -0.018546097, 0.013328026][7:end]
    wf_y_ref = RDWaveform(wf_x.time[7:end], y_ref)

    X = wfs_x.signal
    Y_ref = ArrayOfSimilarArrays([
        Float32[-0.058776174, -0.15214051, -0.19546856, -0.1468839, -0.027337808, 0.07304379, 0.009391554, 0.008527434, -0.0020374788, -0.0010664854, -0.010996632, -0.006111667, -0.008577821, 0.06304957, 0.13707972, 0.17107229, 0.13423079, 0.046884775, -0.07579133, 0.011823545, 7.777956f-5, 0.005714143, -0.009477851, -0.096994996, -0.29966956, -0.34484184, -0.2542256, -0.040562287, 0.110750146, -0.030599957, -0.01432999, 0.019098327, 0.010801469, -0.008182973, -0.0030080022, 0.004854781, -0.044890568, -0.050525106, 0.008078488, 0.030457651, 0.03784929, 0.03058027, -0.018897993, -0.02196053, -0.04128293, -0.018546097, 0.013328026][7:end],
        Float32[-0.058776174, -0.15214051, -0.19546856, -0.1468839, -0.027337808, 0.07304379, 0.009391554, 0.008527434, -0.0020374788, -0.0010664854, -0.010996632, -0.006111667, 0.048961863, 0.1543194, 0.17477812, 0.12345324, 0.025103811, -0.044385068, -0.014283393, 0.011823545, 7.777956f-5, 0.005714143, -0.009477851, -0.096994996, -0.29966956, -0.34484184, -0.2542256, -0.040562287, 0.110750146, -0.030599957, -0.01432999, 0.019098327, 0.010801469, -0.008182973, -0.0030080022, 0.004854781, -0.044890568, -0.050525106, 0.008078488, 0.030457651, 0.03784929, 0.03058027, -0.018897993, -0.02196053, -0.04128293, -0.018546097, 0.013328026][7:end],
        Float32[-0.058776174, -0.15214051, -0.19546856, -0.1468839, -0.027337808, 0.07304379, 0.009391554, 0.008527434, -0.0020374788, -0.0010664854, -0.010996632, 0.051428016, 0.1402317, 0.19201782, 0.12715907, 0.014326249, -0.06616603, 0.01712287, -0.014283393, 0.011823545, 7.777956f-5, 0.005714143, -0.009477851, -0.096994996, -0.29966956, -0.34484184, -0.2542256, -0.040562287, 0.110750146, -0.030599957, -0.01432999, 0.019098327, 0.010801469, -0.008182973, -0.0030080022, 0.004854781, -0.044890568, -0.050525106, 0.008078488, 0.030457651, 0.03784929, 0.03058027, -0.018897993, -0.02196053, -0.04128293, -0.018546097, 0.013328026][7:end],
        Float32[-0.058776174, -0.15214051, -0.19546856, -0.1468839, -0.027337808, 0.07304379, 0.009391554, 0.008527434, -0.0020374788, -0.0010664854, 0.04654305, 0.14269786, 0.17793012, 0.14439878, 0.01803209, -0.07694359, -0.004658094, 0.01712287, -0.014283393, 0.011823545, 7.777956f-5, 0.005714143, -0.009477851, -0.096994996, -0.29966956, -0.34484184, -0.2542256, -0.040562287, 0.110750146, -0.030599957, -0.01432999, 0.019098327, 0.010801469, -0.008182973, -0.0030080022, 0.004854781, -0.044890568, -0.050525106, 0.008078488, 0.030457651, 0.03784929, 0.03058027, -0.018897993, -0.02196053, -0.04128293, -0.018546097, 0.013328026][7:end],
        Float32[-0.058776174, -0.15214051, -0.19546856, -0.1468839, -0.027337808, 0.07304379, 0.009391554, 0.008527434, -0.0020374788, 0.056473196, 0.13781288, 0.18039627, 0.13031107, 0.035271794, -0.07323775, -0.015435659, -0.004658094, 0.01712287, -0.014283393, 0.011823545, 7.777956f-5, 0.005714143, -0.009477851, -0.096994996, -0.29966956, -0.34484184, -0.2542256, -0.040562287, 0.110750146, -0.030599957, -0.01432999, 0.019098327, 0.010801469, -0.008182973, -0.0030080022, 0.004854781, -0.044890568, -0.050525106, 0.008078488, 0.030457651, 0.03784929, 0.03058027, -0.018897993, -0.02196053, -0.04128293, -0.018546097, 0.013328026][7:end],
        Float32[-0.058776174, -0.15214051, -0.19546856, -0.1468839, -0.027337808, 0.07304379, 0.009391554, 0.008527434, 0.055502206, 0.14774303, 0.1755113, 0.13277723, 0.021184087, -0.05599805, -0.011729809, -0.015435659, -0.004658094, 0.01712287, -0.014283393, 0.011823545, 7.777956f-5, 0.005714143, -0.009477851, -0.096994996, -0.29966956, -0.34484184, -0.2542256, -0.040562287, 0.110750146, -0.030599957, -0.01432999, 0.019098327, 0.010801469, -0.008182973, -0.0030080022, 0.004854781, -0.044890568, -0.050525106, 0.008078488, 0.030457651, 0.03784929, 0.03058027, -0.018897993, -0.02196053, -0.04128293, -0.018546097, 0.013328026][7:end],
        Float32[-0.058776174, -0.15214051, -0.19546856, -0.1468839, -0.027337808, 0.07304379, 0.009391554, 0.066067114, 0.14677204, 0.18544145, 0.12789226, 0.023650238, -0.07008576, 0.0055098864, -0.011729809, -0.015435659, -0.004658094, 0.01712287, -0.014283393, 0.011823545, 7.777956f-5, 0.005714143, -0.009477851, -0.096994996, -0.29966956, -0.34484184, -0.2542256, -0.040562287, 0.110750146, -0.030599957, -0.01432999, 0.019098327, 0.010801469, -0.008182973, -0.0030080022, 0.004854781, -0.044890568, -0.050525106, 0.008078488, 0.030457651, 0.03784929, 0.03058027, -0.018897993, -0.02196053, -0.04128293, -0.018546097, 0.013328026][7:end],
        Float32[-0.058776174, -0.15214051, -0.19546856, -0.1468839, -0.027337808, 0.07304379, 0.06693123, 0.15733695, 0.18447046, 0.1378224, 0.01876527, -0.0676196, -0.008577819, 0.0055098864, -0.011729809, -0.015435659, -0.004658094, 0.01712287, -0.014283393, 0.011823545, 7.777956f-5, 0.005714143, -0.009477851, -0.096994996, -0.29966956, -0.34484184, -0.2542256, -0.040562287, 0.110750146, -0.030599957, -0.01432999, 0.019098327, 0.010801469, -0.008182973, -0.0030080022, 0.004854781, -0.044890568, -0.050525106, 0.008078488, 0.030457651, 0.03784929, 0.03058027, -0.018897993, -0.02196053, -0.04128293, -0.018546097, 0.013328026][7:end],
        Float32[-0.058776174, -0.15214051, -0.19546856, -0.1468839, -0.027337808, 0.13058348, 0.15820108, 0.19503537, 0.13685142, 0.028695418, -0.07250457, -0.0061116647, -0.008577819, 0.0055098864, -0.011729809, -0.015435659, -0.004658094, 0.01712287, -0.014283393, 0.011823545, 7.777956f-5, 0.005714143, -0.009477851, -0.096994996, -0.29966956, -0.34484184, -0.2542256, -0.040562287, 0.110750146, -0.030599957, -0.01432999, 0.019098327, 0.010801469, -0.008182973, -0.0030080022, 0.004854781, -0.044890568, -0.050525106, 0.008078488, 0.030457651, 0.03784929, 0.03058027, -0.018897993, -0.02196053, -0.04128293, -0.018546097, 0.013328026][7:end],
        Float32[-0.058776174, -0.15214051, -0.19546856, -0.1468839, 0.030201875, 0.22185332, 0.19589949, 0.14741632, 0.027724424, -0.062574424, -0.010996635, -0.0061116647, -0.008577819, 0.0055098864, -0.011729809, -0.015435659, -0.004658094, 0.01712287, -0.014283393, 0.011823545, 7.777956f-5, 0.005714143, -0.009477851, -0.096994996, -0.29966956, -0.34484184, -0.2542256, -0.040562287, 0.110750146, -0.030599957, -0.01432999, 0.019098327, 0.010801469, -0.008182973, -0.0030080022, 0.004854781, -0.044890568, -0.050525106, 0.008078488, 0.030457651, 0.03784929, 0.03058027, -0.018897993, -0.02196053, -0.04128293, -0.018546097, 0.013328026][7:end]
    ])
    wfs_y_ref = ArrayOfRDWaveforms((map(t -> t[7:end], wfs_x.time), Y_ref))
 
    fi = @inferred fltinstance(flt, smplinfo(x))
    y = similar(x, length(x) - 6)
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
    @test @inferred(flt(wf_x)) ≈ wf_y_ref

    @test DSP.FIRFilter(flt) isa DSP.FIRFilter # No inference test, FIRFilter ctor is not type stable
    @test DSP.filt(DSP.FIRFilter(flt), x)[7:end] ≈ flt(x)

    Y = nestedview(similar(flatview(X), (innersize(X, 1) - 6, size(X, 1))))
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
    @test @inferred(broadcast(flt, wfs_x)) ≈ wfs_y_ref
    #@test typeof(broadcast(fi, wfs_x)) == typeof(wfs_y_ref)
    @test typeof(broadcast(flt, wfs_x)) == typeof(wfs_y_ref)
end
