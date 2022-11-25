# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

include("test_utils.jl")

using IntervalSets


@testset "FixedPickoff" begin
    flt = FixedPickoff(27)
    uflt = TruncateFilter(19.6u"ns"..52.4u"ns")
    wfs_x = gen_test_waveforms()
    Y_ref = flatview(wfs_x.signal)[27, :]

    test_filter(flt, uflt, wfs_x, Y_ref)
end
