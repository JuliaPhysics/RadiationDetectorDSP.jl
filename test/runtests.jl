# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

import Test
import RadiationDetectorDSP
import Documenter

Test.@testset "Package RadiationDetectorDSP" begin
    include("legacy/test_filters.jl")
    include("test_biquad_filter.jl")
    include("test_circuit_filters.jl")

    # doctests
    Documenter.DocMeta.setdocmeta!(
        RadiationDetectorDSP,
        :DocTestSetup,
        :(using RadiationDetectorDSP);
        recursive=true,
    )
    Documenter.doctest(RadiationDetectorDSP)
end # testset
