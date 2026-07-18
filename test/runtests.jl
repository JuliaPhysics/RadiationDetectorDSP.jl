# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

import Test
import RadiationDetectorDSP
import Documenter

Test.@testset "Package RadiationDetectorDSP" begin
    include("test_array_utils.jl")
    include("test_fast_indexing.jl")
    include("test_samples.jl")
    include("test_transpose.jl")
    include("legacy/test_filters.jl")
    include("legacy/test_trapezoidal_filter.jl")
    include("test_global_ops.jl")
    include("test_truncate_filter.jl")
    include("test_convolution_filter.jl")
    include("test_biquad_filter.jl")
    include("test_first_order_iir.jl")
    include("test_circuit_filters.jl")
    include("test_trapezoidal_filter.jl")
    include("test_sg_filter.jl")
    include("test_zac_filter.jl")
    include("test_cusp_filter.jl")
    include("test_gaussian_filter.jl")
    include("test_signal_stats.jl")
    include("test_intersect.jl")
    include("test_signal_estimator.jl")
    include("test_ka_jlarrays.jl")
    include("test_rand_waveforms.jl")

    # doctests
    Documenter.DocMeta.setdocmeta!(
        RadiationDetectorDSP,
        :DocTestSetup,
        :(using RadiationDetectorDSP);
        recursive=true,
    )
    Documenter.doctest(RadiationDetectorDSP)
end # testset
