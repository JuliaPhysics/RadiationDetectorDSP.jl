# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    abstract type AbstractRadIIRFilter <: AbstractRadSigFilter{LinearFiltering}

Abstract type for IIR filters.
"""
abstract type AbstractRadIIRFilter <: AbstractRadSigFilter{LinearFiltering} end
export AbstractRadIIRFilter



"""
    struct BiquadFilter{T<:Real} <: AbstractRadIIRFilter

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


function DSP.Biquad{T}(flt::BiquadFilter) where {T<:Real}
    DSP.Biquad(map(T, flt.b_012)..., map(T, flt.a_12)...)
end


# For testing:
fltparameters(f::BiquadFilter) = (b_012 = f.b_012, a_12 = f.a_12)
fltparameters(f::DSP.Biquad) = (b_012 = SVec(f.b0, f.b1, f.b2), a_12 = SVec(one(f.a1), f.a1, f.a2))


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
    struct RCFilter{T<:Real} <: AbstractRadIIRFilter

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



"""
    struct InvRCFilter{T<:Real} <: AbstractRadIIRFilter

A simple RC-filter.

Constructors:

* ```$(FUNCTIONNAME)(fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct InvRCFilter{T<:RealQuantity} <: AbstractRadIIRFilter
    "RC time constant"
    rc::T
end

export InvRCFilter

InverseFunctions.inverse(flt::InvRCFilter) = RCFilter(flt.rc)

function BiquadFilter(flt::InvRCFilter)
    RC = float(flt.rc)
    k = 1 + RC
    T = typeof(k)
    BiquadFilter((k, T(1) - k, T(0)), (T(0), T(0)))
end



"""
    struct CRFilter{T<:Real} <: AbstractRadIIRFilter

A simple CR-filter.

Constructors:

* ```$(FUNCTIONNAME)(fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct CRFilter{T<:RealQuantity} <: AbstractRadIIRFilter
    "CR time constant"
    cr::T
end

export CRFilter

InverseFunctions.inverse(flt::CRFilter) = InvCRFilter(flt.cr)

function BiquadFilter(flt::CRFilter)
    CR = float(flt.cr)
    α = CR / (CR + 1)
    T = typeof(α)
    BiquadFilter((α, -α, T(0)), (-α, T(0)))
end



"""
    struct InvCRFilter{T<:Real} <: AbstractRadIIRFilter

Inverse of [`CRFilter`](@ref).

Constructors:

* ```$(FUNCTIONNAME)(fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct InvCRFilter{T<:RealQuantity} <: AbstractRadIIRFilter
    "CR time constant"
    cr::T
end

export InvCRFilter

InverseFunctions.inverse(flt::InvCRFilter) = CRFilter(flt.cr)

function BiquadFilter(flt::InvCRFilter)
    CR = float(flt.cr)
    α = inv(1 + CR)
    k = 1 + inv(CR)
    Biquad((k, k * (α - 1), 0), (-1, 0))
end




# ToDo:

#=
function crmod_filter(RC::Real)
    T = float(typeof(RC))
    α = RC / (RC + 1)
    Biquad(T(1), T(-1), T(0), T(-α), T(0))
end


function inv_crmod_filter(RC::Real)
    T = float(typeof(RC))
    α = 1 / (1 + RC)
    Biquad(T(1), T(α - 1), T(0), T(-1), T(0))
end


function integrator_filter(gain::Real)
    T = float(typeof(gain))
    Biquad(T(gain), T(0), T(0), T(-1), T(0))
end


function differentiator_filter(gain::Real)
    T = float(typeof(gain))
    Biquad(T(gain), T(-gain), T(0), T(0), T(0))
end


function integrator_cr_filter(gain::Real, RC::Real)
    T = float(promote_type(typeof(gain), typeof(RC)))
    α = 1 / (1 + RC)
    Biquad(T(gain), T(-α), T(0), T(α - 1), T(0))
end


function integrator_crmod_filter(gain::Real, RC::Real)
    T = float(promote_type(typeof(gain), typeof(RC)))
    α = 1 / (1 + RC)
    Biquad(T(gain), T(0), T(0), T(α - 1), T(0))
end


function simple_csa_response_filter(τ_rise::Real, τ_decay::Real, gain::Real = one(τ_rise))
    # TODO: Use a single biquad filter

    T = float(promote_type(promote_type(typeof(τ_rise), typeof(τ_decay)), typeof(gain)))
    rc_filter(T(τ_rise)) * integrator_cr_filter(T(gain), T(τ_decay))
end
=#
