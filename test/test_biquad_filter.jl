# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using Unitful
using RadiationDetectorSignals, ArraysOfArrays, FillArrays
using InverseFunctions
using Adapt

using RadiationDetectorDSP: bc_rdfilt, bc_rdfilt!


@testset "BiquadFilter" begin
    flt = BiquadFilter((0.2, 0.15, 0.3), (-0.8, 0.4))
        
    @test adapt(Array, flt) isa BiquadFilter

    wfs_x = ArrayOfRDWaveforms((
        Fill(1.5u"ns" .* (1:47), 10),
        ArrayOfSimilarArrays([vcat(fill(0.5f0, 14 - i), zeros(Float32, 9 + i),ones(Float32, 24)) for i in 1:10])
    ))
    wf_x = wfs_x[1]
    x = wf_x.signal
    # From DSP.filt(DSP.Biquad(flt.b_012..., flt.a_12...), x):
    y_ref = [0.1, 0.255, 0.489, 0.6142, 0.62076, 0.575928, 0.5374384, 0.52457952, 0.529688256, 0.5389187968, 0.54425973504, 0.544840269312, 0.5431683214336, 0.44159854942208004, 0.286011510964224, 0.05216978900254721, -0.07266877318365184, -0.07900293414794035, -0.03413483804489155, 0.0042933032232628995, 0.01708857779656694, 0.011953540947948394, 0.0027274016397319393, 0.1974005049326062, 0.5068294432901922, 0.9765033526591114, 1.2284709048112121, 1.242175382785325, 1.152351944303775, 1.0750114023288901, 1.0490683441416022, 1.0592501143817257, 1.0777727538487396, 1.0885181573263014, 1.0897054243215452, 1.0863570765267156, 1.0832034914927544, 1.0820199625835174, 1.0823345734697123, 1.0830596737423628, 1.0835139096060054, 1.0835872581878592, 1.083464242707885, 1.0833364908911645, 1.0832834956297777, 1.0832922001473564, 1.0833203618659741]
    wf_y_ref = RDWaveform(wf_x.time, y_ref)

    X = wfs_x.signal
    Y_ref = ArrayOfSimilarArrays([
        Float32[0.1, 0.255, 0.489, 0.6142, 0.62076, 0.575928, 0.5374384, 0.5245795, 0.52968824, 0.5389188, 0.5442597, 0.5448403, 0.5431683, 0.44159853, 0.28601152, 0.05216979, -0.072668776, -0.07900293, -0.03413484, 0.004293303, 0.017088577, 0.011953541, 0.0027274017, 0.19740051, 0.50682944, 0.9765034, 1.2284709, 1.2421753, 1.152352, 1.0750114, 1.0490683, 1.0592501, 1.0777727, 1.0885181, 1.0897055, 1.0863571, 1.0832034, 1.0820199, 1.0823345, 1.0830597, 1.0835139, 1.0835873, 1.0834643, 1.0833365, 1.0832835, 1.0832922, 1.0833204],
        Float32[0.1, 0.255, 0.489, 0.6142, 0.62076, 0.575928, 0.5374384, 0.5245795, 0.52968824, 0.5389188, 0.5442597, 0.5448403, 0.4431683, 0.28659856, 0.052011512, -0.07303021, -0.07922877, -0.034170933, 0.004354762, 0.017152183, 0.011979842, 0.0027230002, -0.0026135365, 0.19681998, 0.5085014, 0.9780731, 1.2290579, 1.2420171, 1.1519905, 1.0747856, 1.0490322, 1.0593116, 1.0778364, 1.0885445, 1.089701, 1.086343, 1.083194, 1.082018, 1.0823368, 1.0830623, 1.083515, 1.0835872, 1.0834637, 1.0833361, 1.0832834, 1.0832922, 1.0833205],
        Float32[0.1, 0.255, 0.489, 0.6142, 0.62076, 0.575928, 0.5374384, 0.5245795, 0.52968824, 0.5389188, 0.5442597, 0.44484028, 0.2881683, 0.05259855, -0.07318849, -0.07959021, -0.03439677, 0.0043186657, 0.017213643, 0.012043447, 0.002749301, -0.002617938, -0.003194071, 0.19849192, 0.51007116, 0.97866017, 1.2288997, 1.2416557, 1.1517646, 1.0747495, 1.0490937, 1.0593752, 1.0778626, 1.0885401, 1.089687, 1.0863335, 1.0831921, 1.0820202, 1.0823394, 1.0830634, 1.0835149, 1.0835866, 1.0834633, 1.083336, 1.0832834, 1.0832924, 1.0833205],
        Float32[0.1, 0.255, 0.489, 0.6142, 0.62076, 0.575928, 0.5374384, 0.5245795, 0.52968824, 0.5389188, 0.44425973, 0.28984028, 0.05416832, -0.07260145, -0.07974849, -0.03475821, 0.004092827, 0.017177546, 0.012104906, 0.0028129064, -0.0025916372, -0.0031984723, -0.001522123, 0.2000617, 0.5106582, 0.97850186, 1.2285383, 1.2414298, 1.1517286, 1.074811, 1.0491573, 1.0594015, 1.0778582, 1.088526, 1.0896775, 1.0863316, 1.0831943, 1.0820228, 1.0823405, 1.0830632, 1.0835145, 1.0835862, 1.0834632, 1.0833361, 1.0832835, 1.0832925, 1.0833205],
        Float32[0.1, 0.255, 0.489, 0.6142, 0.62076, 0.575928, 0.5374384, 0.5245795, 0.52968824, 0.4389188, 0.28925973, 0.05584027, -0.07103168, -0.07916145, -0.03491649, 0.003731389, 0.016951706, 0.01206881, 0.002874365, -0.0025280318, -0.0031721715, -0.0015265244, 4.7649017f-5, 0.20064873, 0.5104999, 0.9781405, 1.2283124, 1.2413937, 1.15179, 1.0748745, 1.0491836, 1.0593971, 1.0778443, 1.0885166, 1.0896755, 1.0863339, 1.0831969, 1.082024, 1.0823404, 1.0830628, 1.0835141, 1.0835861, 1.0834633, 1.0833362, 1.0832837, 1.0832925, 1.0833205],
        Float32[0.1, 0.255, 0.489, 0.6142, 0.62076, 0.575928, 0.5374384, 0.5245795, 0.42968825, 0.2839188, 0.055259734, -0.06935973, -0.07759168, -0.03432945, 0.003573111, 0.01659027, 0.011842971, 0.002838269, -0.0024665731, -0.0031085662, -0.0015002236, 4.3247524f-5, 0.0006346875, 0.20049044, 0.5101385, 0.97791463, 1.2282763, 1.2414552, 1.1518537, 1.0749009, 1.0491792, 1.059383, 1.0778347, 1.0885146, 1.0896778, 1.0863364, 1.083198, 1.0820239, 1.0823399, 1.0830624, 1.083514, 1.0835862, 1.0834634, 1.0833362, 1.0832837, 1.0832924, 1.0833205],
        Float32[0.1, 0.255, 0.489, 0.6142, 0.62076, 0.575928, 0.5374384, 0.42457953, 0.27468824, 0.049918797, -0.06994026, -0.07591973, -0.032759678, 0.0041601495, 0.016431991, 0.011481533, 0.00261243, -0.002502669, -0.0030471073, -0.0014366182, 6.9548376f-5, 0.000630286, 0.00047640945, 0.20012902, 0.50991267, 0.9778785, 1.2283378, 1.2415189, 1.1518799, 1.0748965, 1.0491651, 1.0593736, 1.0778328, 1.0885168, 1.0896803, 1.0863376, 1.083198, 1.0820233, 1.0823395, 1.0830623, 1.083514, 1.0835863, 1.0834634, 1.0833362, 1.0832837, 1.0832924, 1.0833205],
        Float32[0.1, 0.255, 0.489, 0.6142, 0.62076, 0.575928, 0.4374384, 0.26957953, 0.040688258, -0.0752812, -0.07650027, -0.03108773, 0.0057299216, 0.01701903, 0.011323255, 0.0022509922, -0.0027285083, -0.0030832035, -0.0013751595, 0.00013315381, 0.00065658684, 0.00047200796, 0.000114971626, 0.19990318, 0.50987655, 0.97793996, 1.2284013, 1.2415451, 1.1518755, 1.0748824, 1.0491557, 1.0593716, 1.077835, 1.0885193, 1.0896815, 1.0863374, 1.0831974, 1.0820229, 1.0823394, 1.0830623, 1.0835141, 1.0835863, 1.0834634, 1.0833362, 1.0832835, 1.0832924, 1.0833205],
        Float32[0.1, 0.255, 0.489, 0.6142, 0.62076, 0.475928, 0.2824384, 0.03557952, -0.08451174, -0.0818412, -0.031668264, 0.007401869, 0.018588802, 0.011910293, 0.0020927142, -0.003089946, -0.0033090424, -0.0014112556, 0.00019461253, 0.00072019227, 0.0004983088, 0.00011057013, -0.00011086741, 0.19986708, 0.509938, 0.97800356, 1.2284276, 1.2415407, 1.1518615, 1.074873, 1.0491537, 1.0593739, 1.0778376, 1.0885205, 1.0896814, 1.0863369, 1.083197, 1.0820228, 1.0823394, 1.0830624, 1.0835142, 1.0835863, 1.0834634, 1.0833362, 1.0832835, 1.0832924, 1.0833205],
        Float32[0.1, 0.255, 0.489, 0.6142, 0.52076, 0.320928, 0.0484384, -0.08962048, -0.09107175, -0.037009202, 0.006821335, 0.02026075, 0.013480065, 0.0026797527, -0.003248224, -0.0036704803, -0.0016370947, 0.00015851643, 0.000781651, 0.00056191423, 0.00013687098, -0.00011526891, -0.00014696352, 0.19992854, 0.5100016, 0.9780299, 1.2284232, 1.2415266, 1.151852, 1.074871, 1.049156, 1.0593764, 1.0778388, 1.0885204, 1.0896809, 1.0863365, 1.0831969, 1.0820229, 1.0823395, 1.0830625, 1.0835142, 1.0835863, 1.0834634, 1.0833362, 1.0832835, 1.0832924, 1.0833205]
    ])
    wfs_y_ref = ArrayOfRDWaveforms((deepcopy(wfs_x.time), Y_ref))
 
    fi = @inferred fltinstance(flt, smplinfo(x))
    y = similar(x)
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

    @test isapprox(@inferred(inverse(flt)(y)), x, rtol = 1e-3)

    #@test @inferred(DSP.Biquad(flt)) isa DSP.Biquad
    #@test DSP.filt(DSP.Biquad(flt), x) ≈ flt(x)

    Y = similar(X)
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
