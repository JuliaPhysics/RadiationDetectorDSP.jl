# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    struct FirstOrderIIR{T<:RealQuantity} <: AbstractRadIIRFilter

A [biquad filter](https://en.wikipedia.org/wiki/Digital_biquad_filter).

Constructors:

* ```$(FUNCTIONNAME)(fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct FirstOrderIIR{T<:RealQuantity} <: AbstractRadIIRFilter
    "Coefficients b_0 to b_1"
    b_01::NTuple{2,T}

    "Coefficient a_1, a_0 equals one implicitly"
    a_1::NTuple{1,T}
end

export FirstOrderIIR


function fltinstance(flt::FirstOrderIIR, si::SamplingInfo{T}) where T
    FirstOrderIIRInstance{T}(flt.b_01, flt.a_1, _smpllen(si))
end


function InverseFunctions.inverse(flt::FirstOrderIIR)
    # In direct form 1:
    # y_i[i] = b0 * x_i[i] + b1 * x_i[i-1]  * x_i[i-2] - a1 * y_i[i-1] - a2 * y_i[i-2]
    # x_i[i] = 1/b0 * y_i[i] + a1/b0 * y_i[i-1] + a2/b0 * y_i[i-2] - b1/b0 * x_i[i-1] /b0 * x_i[i-2]

    b0, b1 = flt.b_01
    a1 = flt.a_1[1]
    inv_b0 = inv(b0)

    FirstOrderIIR((inv_b0, inv_b0 * a1), (inv_b0 * b1,))
end


function Base.:(∘)(f::FirstOrderIIR, g::FirstOrderIIR)
    f_b0, f_b1 = f.b_01
    f_a1 = f.a_1[1]
    g_b0, g_b1 = g.b_01
    g_a1 = g.a_1[1]

    b0 = f_b0 * g_b0
    b1 = f_b0 * g_b1 + f_b1 * g_b0
    b2 = f_b1 * g_b1
    
    a1 = f_a1 + g_a1
    a2 = f_a1 * g_a1
    BiquadFilter((b0, b1, b2), (a1, a2))
end


function BiquadFilter(flt::FirstOrderIIR{T}) where {T<:RealQuantity}
    BiquadFilter((flt.b_01..., T(0)), (flt.a_1..., T(0)))
end

function DSP.Biquad(flt::FirstOrderIIR{T}) where {T<:RealQuantity}
    DSP.Biquad(map(T, flt.b_01)..., T(0), map(T, flt.a_1)..., T(0))
end



struct FirstOrderIIRInstance{T} <: AbstractRadSigFilterInstance{LinearFiltering}
    b_01::NTuple{2,T}
    a_1::NTuple{1,T}
    n::Int
end

@inline function rdfilt!(Y::AbstractVector{T}, fi::FirstOrderIIRInstance{T}, X::AbstractVector{T}) where {T<:Real}
    a1 = fi.a_1[1]
    neg_a1 = -a1
    b0, b1 = fi.b_01
    s1::T = zero(T) # s_init[1]
    s2::T = zero(T)

    #!!! @assert eachindex(X) == eachindex(Y)

    @inbounds @simd for i in eachindex(X)
        x_i = T(X[i])

        z1 = fma(b0, x_i, s1)
        z2 = fma(b1, x_i, s2)

        y_i = z1
        Y[i] = y_i

        s1 = fma(neg_a1, y_i, z2)
    end
    Y
end


flt_output_smpltype(fi::FirstOrderIIRInstance) = flt_input_smpltype(fi)
flt_input_smpltype(fi::FirstOrderIIRInstance{T}) where T = T
flt_output_length(fi::FirstOrderIIRInstance) = flt_input_length(fi)
flt_input_length(fi::FirstOrderIIRInstance) = fi.n
flt_output_time_axis(fi::FirstOrderIIRInstance, time::AbstractVector{<:RealQuantity}) = time
