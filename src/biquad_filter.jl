# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    abstract type AbstractRadIIRFilter <: AbstractRadSigFilter{LinearFiltering}

Abstract type for IIR filters.
"""
abstract type AbstractRadIIRFilter <: AbstractRadSigFilter{LinearFiltering} end
export AbstractRadIIRFilter



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

## For testing:
#fltparameters(f::BiquadFilter) = (b_012 = f.b_012, a_12 = f.a_12)
#fltparameters(f::DSP.Biquad) = (b_012 = SVec(f.b0, f.b1, f.b2), a_12 = SVec(one(f.a1), f.a1, f.a2))


function DSP.Biquad(flt::BiquadFilter{T}) where {T<:RealQuantity}
    DSP.Biquad(map(T, flt.b_012)..., map(T, flt.a_12)...)
end


function InverseFunctions.inverse(flt::BiquadFilter)
    # In direct form 1:
    # y[i] = b0 * x[i] + b1 * x[i-1] + b2 * x[i-2] - a1 * y[i-1] - a2 * y[i-2]
    # x[i] = 1/b0 * y[i] + a1/b0 * y[i-1] + a2/b0 * y[i-2] - b1/b0 * x[i-1] - b2/b0 * x[i-2]

    b0, b1, b2 = flt.b_012
    a1, a2 = flt.a_12
    inv_b0 = inv(b0)

    BiquadFilter((inv_b0, inv_b0 * a1, inv_b0 * a2), (inv_b0 * b1, inv_b0 * b2))
end



"""
    struct RCFilter{T<:RealQuantity} <: AbstractRadIIRFilter

A simple RC-filter.

Constructors:

* ```$(FUNCTIONNAME)(fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct RCFilter{T<:RealQuantity} <: AbstractRadIIRFilter
    "RC time constant"
    rc::T
end

export RCFilter

InverseFunctions.inverse(flt::RCFilter) = InvRCFilter(flt.rc)

function BiquadFilter(flt::RCFilter)
    RC = float(flt.rc)
    α = 1 / (1 + RC)
    T = typeof(α)
    BiquadFilter((α, T(0), T(0)), (α - T(1), T(0)))
end


struct BiquadFilterInstance{T} <: AbstractRadSigFilter{LinearFiltering}
    b_012::NTuple{3,T}
    a_12::NTuple{2,T}
    n::Int
end

@inline function rdfilt!(Y::AbstractVector{T}, fi::BiquadFilterInstance{T}, X::AbstractVector{T}) where {T<:Real}
    a1, a2 = fi.a_12
    neg_a1, neg_a2 = -a1, -a2
    b0, b1, b2 = fi.b_012
    s1::T = zero(U) # s_init[1]
    s2::T = zero(U) # s_init[2]
    s3::T = zero(U)

    # @assert eachindex(X) == eachindex(Y)

    @inbounds @simd for i in eachindex(X)
        x = U(X[i])

        z1 = fma(b0, x, s1)
        z2 = fma(b1, x, s2)
        z3 = fma(b2, x, s3)

        y = z1
        Y[i] = y

        s1 = fma(neg_a1, y, z2)
        s2 = fma(neg_a2, y, z3)
    end
    Y
end



flt_output_smpltype(fi::BiquadFilterInstance) = flt_input_smpltype(fi)
flt_input_smpltype(fi::BiquadFilterInstance{T}) where T = T
flt_output_length(fi::BiquadFilterInstance) = flt_input_length(fi)
flt_input_length(fi::BiquadFilterInstance) = fi.n
flt_output_time_axis(fi::BiquadFilterInstance, time::AbstractVector{<:RealQuantity}) = time
