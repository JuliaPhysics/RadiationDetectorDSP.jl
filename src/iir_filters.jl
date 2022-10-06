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


# inverse:
# y[i] = b0 * x[i] + b1 * x[i-1] + b2 * x[i-2] - a1 * y[i-1] - a2 * y[i-2]
# x[i] = 1 * y[i] + a1/b0 * y[i-1] + a2/b0 * y[i-2] - b1/b0 * x[i-1] - b2/b0 * x[i-2]




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

function BiquadFilter(flt::RCFilter)
    RC = float(flt.RC)
    T = typeof(RC)
    α = 1 / (1 + RC)
    BiquadFilter((α, 0, 0), (α - 1, 0))
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

Base.inv(flt::CRFilter) = InvCRFilter(flt.cr)

function BiquadFilter(flt::CRFilter)
    CR = float(flt.cr)
    α = CR / (CR + 1)
    BiquadFilter((α, -α, 0), (-α, 0))
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

Base.inv(flt::InvCRFilter) = CRFilter(flt.cr)

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
