# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    struct RCFilter <: AbstractRadII>RFilter

A first-order RC lowpass filter.

The inverse filter is [`InvCRFilter`](@ref), but note that this is unstable
in the presence of additional noise and so will typically not be useful to
deconvolve signals in practical applications.

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

function fltinstance(flt::RCFilter, fi::SamplingInfo)
    fltinstance(FirstOrderIIR(RCFilter(ustrip(NoUnits, flt.rc / step(fi.axis)))), fi)
end

InverseFunctions.inverse(flt::RCFilter) = InvRCFilter(flt.rc)

function FirstOrderIIR(flt::RCFilter)
    RC = float(flt.rc)
    α = 1 / (1 + RC)
    T = typeof(α)
    FirstOrderIIR((α, T(0)), (α - T(1),))
end



"""
    struct InvRCFilter <: AbstractRadIIRFilter

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

function fltinstance(flt::InvRCFilter, fi::SamplingInfo)
    fltinstance(FirstOrderIIR(InvRCFilter(ustrip(NoUnits, flt.rc / step(fi.axis)))), fi)
end

InverseFunctions.inverse(flt::InvRCFilter) = RCFilter(flt.rc)

function FirstOrderIIR(flt::InvRCFilter)
    RC = float(flt.rc)
    k = 1 + RC
    T = typeof(k)
    FirstOrderIIR((k, T(1) - k), (T(0),))
end



"""
    struct CRFilter <: AbstractRadIIRFilter

A first-order CR highpass filter.

The inverse filter is [`InvCRFilter`](@ref), this is typically stable even in
the presence of additional noise.

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

function fltinstance(flt::CRFilter, fi::SamplingInfo)
    fltinstance(FirstOrderIIR(CRFilter(ustrip(NoUnits, flt.cr / step(fi.axis)))), fi)
end

InverseFunctions.inverse(flt::CRFilter) = InvCRFilter(flt.cr)

function FirstOrderIIR(flt::CRFilter)
    CR = float(flt.cr)
    α = CR / (CR + 1)
    FirstOrderIIR((α, -α), (-α,))
end



"""
    struct InvCRFilter <: AbstractRadIIRFilter

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

function fltinstance(flt::InvCRFilter, fi::SamplingInfo)
    fltinstance(FirstOrderIIR(InvCRFilter(ustrip(NoUnits, flt.cr / step(fi.axis)))), fi)
end

InverseFunctions.inverse(flt::InvCRFilter) = CRFilter(flt.cr)

function FirstOrderIIR(flt::InvCRFilter)
    CR = float(flt.cr)
    k = 1 + inv(CR) # equivalent to k = -1 / (α - 1)
    T = typeof(k)
    FirstOrderIIR((k, T(-1)), (T(-1),))
end


"""
    struct ModCRFilter <: AbstractRadIIRFilter

A first-order CR highpass filter, modified for full-amplitude step-signal
response.

The resonse of the standard digital [`CRFilter`](@ref) will not recover the
full amplitude of a digital step stignal since a step from one sample to the
still has a finite rise time. This version of a CR filter compensates for
this loss in amplitude, so it effectively treats a step as having

The inverse filter is [`InvModCRFilter`](@ref), this is typically stable even in
the presence of additional noise.

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

function fltinstance(flt::ModCRFilter, fi::SamplingInfo)
    fltinstance(FirstOrderIIR(ModCRFilter(ustrip(NoUnits, flt.cr / step(fi.axis)))), fi)
end

InverseFunctions.inverse(flt::ModCRFilter) = InvModCRFilter(flt.cr)

function FirstOrderIIR(flt::ModCRFilter)
    CR = float(flt.cr)
    k = CR / (CR + 1)
    T = typeof(k)
    FirstOrderIIR((T(1), T(-1)), (-k,))
end


"""
    struct InvModCRFilter <: AbstractRadIIRFilter

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

function fltinstance(flt::InvModCRFilter, fi::SamplingInfo)
    fltinstance(FirstOrderIIR(InvModCRFilter(ustrip(NoUnits, flt.cr / step(fi.axis)))), fi)
end

InverseFunctions.inverse(flt::InvModCRFilter) = ModCRFilter(flt.cr)

function FirstOrderIIR(flt::InvModCRFilter)
    CR = float(flt.cr)
    α = 1 / (1 + CR)
    T = typeof(α)
    FirstOrderIIR((T(1), α - T(1)), (T(-1),))
end



"""
    struct IntegratorFilter <: AbstractRadIIRFilter

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

fltinstance(flt::IntegratorFilter, fi::SamplingInfo) = fltinstance(FirstOrderIIR(flt), fi)

InverseFunctions.inverse(flt::IntegratorFilter) = DifferentiatorFilter(inv(flt.gain))

function FirstOrderIIR(flt::IntegratorFilter)
    g = flt.gain
    T = typeof(g)
    FirstOrderIIR((g, T(0)), (T(-1),))
end



"""
    struct DifferentiatorFilter <: AbstractRadIIRFilter

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

fltinstance(flt::DifferentiatorFilter, fi::SamplingInfo) = fltinstance(FirstOrderIIR(flt), fi)

InverseFunctions.inverse(flt::DifferentiatorFilter) = IntegratorFilter(inv(flt.gain))

function FirstOrderIIR(flt::DifferentiatorFilter)
    g = flt.gain
    T = typeof(g)
    FirstOrderIIR((g, -g), (T(0),))
end



"""
    struct IntegratorCRFilter <: AbstractRadIIRFilter

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

function fltinstance(flt::IntegratorCRFilter, fi::SamplingInfo)
    fltinstance(FirstOrderIIR(IntegratorCRFilter(flt.gain, ustrip(NoUnits, flt.cr / step(fi.axis)))), fi)
end

InverseFunctions.inverse(flt::IntegratorCRFilter) = inverse(FirstOrderIIR(flt))

function FirstOrderIIR(flt::IntegratorCRFilter)
    CR = float(flt.cr)
    α = 1 / (1 + CR)
    T = typeof(α)
    g = T(flt.gain)
    FirstOrderIIR((g, -α), (α - T(1),))
end



"""
    struct IntegratorModCRFilter <: AbstractRadIIRFilter

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

function fltinstance(flt::IntegratorModCRFilter, fi::SamplingInfo)
    fltinstance(FirstOrderIIR(IntegratorCRFilter(flt.gain, ustrip(NoUnits, flt.cr / step(fi.axis)))), fi)
end

InverseFunctions.inverse(flt::IntegratorModCRFilter) = inverse(FirstOrderIIR(flt))

function BiquadFirstOrderIIRFilter(flt::IntegratorModCRFilter)
    CR = float(flt.cr)
    α = 1 / (1 + CR)
    T = typeof(α)
    g = T(flt.gain)
    FirstOrderIIR((g, T(0)), (α - T(1)))
end



"""
    struct SimpleCSAFilter <: AbstractRadIIRFilter

Simulates the current-signal response of a charge-sensitive preamplifier with
resistive reset, the output is a charge signal.

It is equivalent to the composition

```julia
CRFilter(cr = tau_decay) ∘
Integrator(gain = gain) ∘
RCFilter(rc = tau_rise)
```

and maps to a single `BiquadFilter`.

This filter has an inverse, but the inverse is very unstable in the presence
of additional noise if `tau_rise` is not zero (since the inverse of an
RC-filter is unstable under noise). Even if `tau_rise` is zero the inverse
will still amplify noise (since it differentiates), so it should be used very
carefully when deconvolving signals in practical applications.

Constructors:

* ```$(FUNCTIONNAME)(fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct SimpleCSAFilter{T<:RealQuantity,U<:RealQuantity} <: AbstractRadIIRFilter
    "Rise time constant"
    tau_rise::T

    "Decay time constant"
    tau_decay::T

    "Gain"
    gain::U = 1
end

export SimpleCSAFilter

function fltinstance(flt::SimpleCSAFilter, fi::SamplingInfo)
    fltinstance(BiquadFilter(SimpleCSAFilter(
        ustrip(NoUnits, flt.tau_rise / step(fi.axis)),
        ustrip(NoUnits, flt.tau_decay / step(fi.axis)),
        flt.gain,
    )), fi)
end

InverseFunctions.inverse(flt::SimpleCSAFilter) = inverse(BiquadFilter(flt))

function BiquadFilter(flt::SimpleCSAFilter)
    flt1 = RCFilter(rc = flt.tau_rise)
    tau_decay, gain = promote(flt.tau_decay, flt.gain)
    flt2 = IntegratorCRFilter(cr = tau_decay, gain = gain)
    FirstOrderIIR(flt1) ∘ FirstOrderIIR(flt2)
end
