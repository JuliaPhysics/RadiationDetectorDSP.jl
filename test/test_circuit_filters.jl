# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using Random, LinearAlgebra, InverseFunctions
using DSP


@testset "circuit_filters" begin
    x = vcat(zeros(23),ones(24))
    flt = RCFilter(7)
    fi = fltinstance(flt, x)
    y = similar(x)
    rdfilt!(y, fi, x)

end
