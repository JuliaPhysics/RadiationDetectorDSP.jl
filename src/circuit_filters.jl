# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    struct RCFilter{T<:RealQuantity} <: AbstractRadIIRFilter

An RC-filter. It's inverse is [`InvRCFilter`](@ref).

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

fltinstance(flt::RCFilter, fi::SamplingInfo) = fltinstance(BiquadFilter(flt), fi)

InverseFunctions.inverse(flt::RCFilter) = InvRCFilter(flt.rc)

function BiquadFilter(flt::RCFilter)
    RC = float(flt.rc)
    α = 1 / (1 + RC)
    T = typeof(α)
    BiquadFilter((α, T(0), T(0)), (α - T(1), T(0)))
end



"""
    struct InvRCFilter{T<:RealQuantity} <: AbstractRadIIRFilter

Inverse of [`RCFilter`](@ref).

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

fltinstance(flt::InvRCFilter, fi::SamplingInfo) = fltinstance(BiquadFilter(flt), fi)

InverseFunctions.inverse(flt::InvRCFilter) = RCFilter(flt.rc)

function BiquadFilter(flt::InvRCFilter)
    RC = float(flt.rc)
    k = 1 + RC
    T = typeof(k)
    BiquadFilter((k, T(1) - k, T(0)), (T(0), T(0)))
end



"""
    struct CRFilter{T<:RealQuantity} <: AbstractRadIIRFilter

An RC-filter. It's inverse is [`InvCRFilter`](@ref).

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

fltinstance(flt::CRFilter, fi::SamplingInfo) = fltinstance(BiquadFilter(flt), fi)

InverseFunctions.inverse(flt::CRFilter) = InvCRFilter(flt.cr)

function BiquadFilter(flt::CRFilter)
    CR = float(flt.cr)
    α = CR / (CR + 1)
    T = typeof(α)
    BiquadFilter((α, -α, T(0)), (-α, T(0)))
end



"""
    struct InvCRFilter{T<:RealQuantity} <: AbstractRadIIRFilter

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

fltinstance(flt::InvCRFilter, fi::SamplingInfo) = fltinstance(BiquadFilter(flt), fi)

InverseFunctions.inverse(flt::InvCRFilter) = CRFilter(flt.cr)

function BiquadFilter(flt::InvCRFilter)
    CR = float(flt.cr)
    α = inv(1 + CR)
    k = 1 + inv(CR)
    Biquad((k, k * (α - 1), 0), (-1, 0))
end


"""
    struct ModCRFilter{T<:RealQuantity} <: AbstractRadIIRFilter

A modified CR-filter. It's inverse is [`InvModCRFilter`](@ref).

Constructors:

* ```$(FUNCTIONNAME)(fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct ModCRFilter{T<:RealQuantity} <: AbstractRadIIRFilter
    "CR time constant"
    cr::T
end

export ModCRFilter

fltinstance(flt::ModCRFilter, fi::SamplingInfo) = fltinstance(BiquadFilter(flt), fi)

InverseFunctions.inverse(flt::ModCRFilter) = InvModCRFilter(flt.cr)

function BiquadFilter(flt::ModCRFilter)
    CR = float(flt.cr)
    α = CR / (CR + 1)
    T = typeof(α)
    BiquadFilter((T(1), T(-1), T(0)), (-α, T(0)))
end


"""
    struct InvModCRFilter{T<:RealQuantity} <: AbstractRadIIRFilter

Inverse of [`ModCRFilter`](@ref).

Constructors:

* ```$(FUNCTIONNAME)(fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct InvModCRFilter{T<:RealQuantity} <: AbstractRadIIRFilter
    "CR time constant"
    cr::T
end

export InvModCRFilter

fltinstance(flt::InvModCRFilter, fi::SamplingInfo) = fltinstance(BiquadFilter(flt), fi)

InverseFunctions.inverse(flt::InvModCRFilter) = ModCRFilter(flt.cr)

function BiquadFilter(flt::InvModCRFilter)
    CR = float(flt.cr)
    k = CR / (CR + 1)
    T = typeof(k)
    Biquad(T(1), T(-1), T(0), -k, T(0))
end



"""
    struct IntegratorFilter{T<:RealQuantity} <: AbstractRadIIRFilter

An integrator filter. It's inverse is [`DifferentiatorFilter`](@ref).

Constructors:

* ```$(FUNCTIONNAME)(fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct IntegratorFilter{T<:RealQuantity} <: AbstractRadIIRFilter
    "Filter gain"
    gain::T = 1
end

export IntegratorFilter

fltinstance(flt::IntegratorFilter, fi::SamplingInfo) = fltinstance(BiquadFilter(flt), fi)

InverseFunctions.inverse(flt::IntegratorFilter) = DifferentiatorFilter(flt.gain)

function BiquadFilter(flt::IntegratorFilter)
    g = flt.gain
    T = typeof(g)
    BiquadFilter((g, T(0), T(0)), (T(-1), T(0)))
end



"""
    struct DifferentiatorFilter{T<:RealQuantity} <: AbstractRadIIRFilter

An integrator filter. It's inverse is [`IntegratorFilter`](@ref).

Constructors:

* ```$(FUNCTIONNAME)(fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct DifferentiatorFilter{T<:RealQuantity} <: AbstractRadIIRFilter
    "Filter gain"
    gain::T = 1
end

export DifferentiatorFilter

fltinstance(flt::DifferentiatorFilter, fi::SamplingInfo) = fltinstance(BiquadFilter(flt), fi)

InverseFunctions.inverse(flt::DifferentiatorFilter) = IntegratorFilter(flt.gain)

function BiquadFilter(flt::DifferentiatorFilter)
    k = flt.gain
    T = typeof(g)
    BiquadFilter((k, -k, T(0)), (T(0), T(0)))
end



"""
    struct IntegratorCRFilter{T<:RealQuantity} <: AbstractRadIIRFilter

A modified CR-filter. The filter has an inverse.

Constructors:

* ```$(FUNCTIONNAME)(fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct IntegratorCRFilter{T<:RealQuantity} <: AbstractRadIIRFilter
    "Filter gain"
    gain::T = 1
    "CR time constant"
    cr::T
end

export IntegratorCRFilter

fltinstance(flt::IntegratorCRFilter, fi::SamplingInfo) = fltinstance(BiquadFilter(flt), fi)

InverseFunctions.inverse(flt::IntegratorCRFilter) = inverse(BiquadFilter(flt))

function BiquadFilter(flt::IntegratorCRFilter)
    CR = float(flt.cr)
    α = 1 / (1 + CR)
    T = typeof(α)
    g = T(flt.gain)
    BiquadFilter((g, -α, T(0)), (α - 1, T(0)))
end



"""
    struct IntegratorModCRFilter{T<:RealQuantity} <: AbstractRadIIRFilter

A modified CR-filter. The filter has an inverse.

Constructors:

* ```$(FUNCTIONNAME)(fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct IntegratorModCRFilter{T<:RealQuantity} <: AbstractRadIIRFilter
    "Filter gain"
    gain::T = 1
    "CR time constant"
    cr::T
end

export IntegratorModCRFilter

fltinstance(flt::IntegratorModCRFilter, fi::SamplingInfo) = fltinstance(BiquadFilter(flt), fi)

InverseFunctions.inverse(flt::IntegratorModCRFilter) = inverse(BiquadFilter(flt))

function BiquadFilter(flt::IntegratorModCRFilter)
    CR = float(flt.cr)
    α = 1 / (1 + CR)
    T = typeof(α)
    g = T(flt.gain)
    BiquadFilter((g, T(0), T(0)), (α - 1, T(0)))
end



"""
    struct SimpleCSAFilter{T<:RealQuantity} <: AbstractRadIIRFilter

A modified CR-filter. The filter has an inverse.

Constructors:

* ```$(FUNCTIONNAME)(fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct SimpleCSAFilter{T<:RealQuantity} <: AbstractRadIIRFilter
    "Rise time constant"
    tau_rise::T

    "Decay time constant"
    tau_decay::T

    "Gain"
    gain::T = 1
end

export SimpleCSAFilter

fltinstance(flt::SimpleCSAFilter, input) = fltinstance(_lower_flt(flt), input)

InverseFunctions.inverse(flt::SimpleCSAFilter) = inverse(_lower_flt(flt))

function _lower_flt(flt::SimpleCSAFilter)
    flt1 = RCFilter(rc = flt.tau_rise)
    flt2 = IntegratorCRFilter(cr = flt.tau_decay, gain = flt.gain)
    FunctionChain((flt1, flt2))
end
