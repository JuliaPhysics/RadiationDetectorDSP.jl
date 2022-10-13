# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    struct BiquadFilter{T<:RealQuantity} <: AbstractRadIIRFilter

A [biquad filter](https://en.wikipedia.org/wiki/Digital_biquad_filter).

Constructors:

* ```$(FUNCTIONNAME)(fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct BiquadFilter{T<:RealQuantity} <: AbstractRadIIRFilter
    "Coefficients b_0 to b_2"
    b_012::NTuple{3,T}

    "Coefficients a_1 to a_2, a_0 equals one implicitly"
    a_12::NTuple{2,T}
end

export BiquadFilter


Adapt.adapt_structure(to, flt::BiquadFilter) = flt


function fltinstance(flt::BiquadFilter, si::SamplingInfo{T}) where T
    BiquadFilterInstance{T}(flt.b_012, flt.a_12, _smpllen(si))
end


function InverseFunctions.inverse(flt::BiquadFilter)
    # In direct form 1:
    # y_i[i] = b0 * x_i[i] + b1 * x_i[i-1] + b2 * x_i[i-2] - a1 * y_i[i-1] - a2 * y_i[i-2]
    # x_i[i] = 1/b0 * y_i[i] + a1/b0 * y_i[i-1] + a2/b0 * y_i[i-2] - b1/b0 * x_i[i-1] - b2/b0 * x_i[i-2]

    b0, b1, b2 = flt.b_012
    a1, a2 = flt.a_12
    inv_b0 = inv(b0)

    BiquadFilter((inv_b0, inv_b0 * a1, inv_b0 * a2), (inv_b0 * b1, inv_b0 * b2))
end


function DSP.Biquad(flt::BiquadFilter{T}) where {T<:RealQuantity}
    DSP.Biquad(map(T, flt.b_012)..., map(T, flt.a_12)...)
end



struct BiquadFilterInstance{T} <: AbstractRadSigFilterInstance{LinearFiltering}
    b_012::NTuple{3,T}
    a_12::NTuple{2,T}
    n::Int
end

@inline function rdfilt!(Y::AbstractVector{T}, fi::BiquadFilterInstance{T}, X::AbstractVector{T}) where {T<:Real}
    a1, a2 = fi.a_12
    neg_a1, neg_a2 = -a1, -a2
    b0, b1, b2 = fi.b_012
    s1::T = zero(T) # s_init[1]
    s2::T = zero(T) # s_init[2]
    s3::T = zero(T)

    @assert eachindex(X) == eachindex(Y)

    @inbounds @simd for i in eachindex(X)
        x_i = T(X[i])

        z1 = fma(b0, x_i, s1)
        z2 = fma(b1, x_i, s2)
        z3 = fma(b2, x_i, s3)

        y_i = z1
        Y[i] = y_i

        s1 = fma(neg_a1, y_i, z2)
        s2 = fma(neg_a2, y_i, z3)
    end
    Y
end

function bc_rdfilt!(
    outputs::ArrayOfSimilarVectors{<:RealQuantity},
    fi::BiquadFilterInstance,
    inputs::ArrayOfSimilarVectors{<:RealQuantity}
)
    _ka_bc_rdfilt!(outputs, fi, inputs)
end


flt_output_smpltype(fi::BiquadFilterInstance) = flt_input_smpltype(fi)
flt_input_smpltype(fi::BiquadFilterInstance{T}) where T = T
flt_output_length(fi::BiquadFilterInstance) = flt_input_length(fi)
flt_input_length(fi::BiquadFilterInstance) = fi.n
flt_output_time_axis(fi::BiquadFilterInstance, time::AbstractVector{<:RealQuantity}) = time
