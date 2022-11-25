# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

include("test_utils.jl")

using IntervalSets


@testset "TruncateFilter" begin
    flt = TruncateFilter(13..35)
    uflt = TruncateFilter(13.4u"ns"..46.6u"ns")
    wfs_x = gen_test_waveforms()
    wfs_y_ref = ArrayOfRDWaveforms((map(t -> t[13:35], wfs_x.time), nestedview(flatview(wfs_x.signal)[13:35, :])))

    test_filter(flt, uflt, wfs_x, wfs_y_ref)
end
