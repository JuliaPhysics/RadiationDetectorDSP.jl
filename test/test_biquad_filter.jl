# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using Random, LinearAlgebra, InverseFunctions
using DSP


@testset "circuit_filters" begin
    x = vcat(zeros(23),ones(24))
    flt = BiquadFilter((0.2, 0.15, 0.3), (-0.8, 0.4))
    # From DSP.filt(DSP.Biquad(flt.b_012..., flt.a_12...), x):
    y_ref = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.51, 0.978, 1.2284, 1.24152, 1.151856, 1.0748768, 1.04915904, 1.059376512, 1.0778375936, 1.08851947008, 1.089680538624, 1.0863366428672, 1.0831970988441602, 1.0820230219284481, 1.0823395780050944, 1.0830624536326963, 1.0835141317041193, 1.0835863239102168, 1.0834634064465258, 1.083336195593134, 1.083283593895897, 1.0832923968794639, 1.0833204799452123]

    fi = @inferred fltinstance(flt, x)
    y = similar(x)
    @test @inferred(rdfilt!(y, fi, x)) ≈ y_ref
    @test y ≈ y_ref
    @test @inferred(rdfilt(fi, x)) ≈ y_ref
    @test @inferred(fi(x)) ≈ y_ref
    @test @inferred(flt(x)) ≈ y_ref

    plot(x); plot!(y)
end
