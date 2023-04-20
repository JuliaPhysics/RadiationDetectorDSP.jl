# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

include("test_utils.jl")


@testset "SavitzkyGolayFilter" begin
    # Reference values from scipy.signal.savgol_coeffs:
    @test isapprox(RadiationDetectorDSP.sg_filter_coeffs(-2:2, 2, 0, 3.2), [-0.08571429,  0.34285714,  0.48571429,  0.34285714, -0.08571429], rtol = 10^-5, atol = 10^-10)
    @test isapprox(RadiationDetectorDSP.sg_filter_coeffs(-2:2, 2, 1, 3.2), [-0.0625, -0.03125, 0.0,  0.03125, 0.0625], rtol = 10^-5, atol = 10^-10)
    @test isapprox(RadiationDetectorDSP.sg_filter_coeffs(-3:3, 3, 2, 1), [0.119047619, 0.0, -0.0714285714, -0.0952380952, -0.0714285714, 0.0, 0.119047619], rtol = 10^-5, atol = 10^-10)
    @test isapprox(RadiationDetectorDSP.sg_filter_coeffs(-3:3, 3, 2, 2.6), [0.0176106, 0.0, -0.0105664, -0.0140885, -0.0105664, 0.0, 0.0176106], rtol = 10^-5, atol = 10^-10)
    @test isapprox(RadiationDetectorDSP.sg_filter_coeffs(-3:3, 3, 4, 2.6), [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], rtol = 10^-5, atol = 10^-10)    

    @test @inferred(RadiationDetectorDSP.sg_filter_coeffs(-2:2, 2, 0)) ≈ @inferred(RadiationDetectorDSP.sg_filter_coeffs(-2:2, 2, 0, 1))
    @test @inferred(RadiationDetectorDSP.sg_filter_coeffs(-2:2, 2, 1)) ≈ @inferred(RadiationDetectorDSP.sg_filter_coeffs(-2:2, 2, 1, 1))
    @test @inferred(RadiationDetectorDSP.sg_filter_coeffs(-3:3, 3, 2)) ≈ @inferred(RadiationDetectorDSP.sg_filter_coeffs(-3:3, 3, 2, 1))
    @test @inferred(RadiationDetectorDSP.sg_filter_coeffs(-4:4, 3, 3)) ≈ @inferred(RadiationDetectorDSP.sg_filter_coeffs(-4:4, 3, 3, 1))

    @test SavitzkyGolayFilter(2, 1, 1)(float.(cumsum(1:5))) ≈ [2.5, 3.5, 4.5]

    flt = SavitzkyGolayFilter(7, 3, 2)
    uflt = SavitzkyGolayFilter(10.3u"ns", 3, 2)
    wfs_x = gen_test_waveforms()

    @test isapprox(ConvolutionFilter(flt).coeffs, [0.119047619, 0.0, -0.0714285714, -0.0952380952, -0.0714285714, 0.0, 0.119047619], rtol = 10^-5, atol = 10^-10)
    @test ConvolutionFilter(flt).offset == -3
    @test uflt(wfs_x[1]).time == wfs_x[1].time[begin+3:end-3]

    wfs_y_ref = ArrayOfRDWaveforms((Fill(wfs_x.time[1][begin+3:end-3], length(wfs_x)), ConvolutionFilter(flt).(wfs_x.signal)))

    test_filter(flt, uflt, wfs_x, wfs_y_ref)
end
