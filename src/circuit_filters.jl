# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    struct RCFilter <: AbstractRadII>RFilter

A first-order RC lowpass filter.

The inverse filter is [`InvCRFilter`](@ref), but note that this is unstable
in the presence of additional noise. As an RC filter attenuates
high-frequency noise, its inverse amplifies such noise and will typically not
be useful to deconvolve signals in practical applications.

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
the presence of additional noise. This is because a CR filter passes
high-frequency noise and so it's inverse passes such noise as well without
amplifying it.

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
the presence of additional noise (see [`CRFilter`](@ref)).

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



"""
    struct SecondOrderCRFilter <: AbstractRadIIRFilter

A scond order CR highpass filter. The filter has an inverse [`InvSecondOrderCRFilter`](@ref).

Constructors:

* ```$(FUNCTIONNAME)(fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct SecondOrderCRFilter{T<:RealQuantity, U<:RealQuantity, V<:Real} <: AbstractRadIIRFilter
    "time constant of the first exponential to be deconvolved"
    cr::T
    "time constant of the second exponential to be deconvolved"
    cr2::U
    "the fraction faktor which the second exponential contributes"
    f::V
end

export SecondOrderCRFilter

function fltinstance(flt::SecondOrderCRFilter, fi::SamplingInfo)
    fltinstance(BiquadFilter(SecondOrderCRFilter(ustrip(NoUnits, flt.cr / step(fi.axis)), ustrip(NoUnits, flt.cr2 / step(fi.axis)), flt.f )), fi)
end

InverseFunctions.inverse(flt::SecondOrderCRFilter) = InvSecondOrderCRFilter(flt.cr, flt.cr2, flt.f)

function BiquadFilter(flt::SecondOrderCRFilter)
    a = exp(-1/float(flt.cr))
    b = exp(-1/float(flt.cr2))
    frac = float(flt.f)
    transfer_denom_1 = frac * b - frac * a - b - 1
    transfer_denom_2 = -(frac * b - frac * a - b)
    transfer_num_1 = -(a + b)
    transfer_num_2 = a * b
    BiquadFilter((1.0, transfer_denom_1, transfer_denom_2), (transfer_num_1, transfer_num_2))
end



"""
    struct InvSecondOrderCRFilter <: AbstractRadIIRFilter

Inverse of [`SecondOrderCRFilter`](@ref).
Apply a double pole-zero cancellation using the provided time constants to the waveform.


Constructors:

* ```$(FUNCTIONNAME)(fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct InvSecondOrderCRFilter{T<:RealQuantity, U<:RealQuantity, V<:Real} <: AbstractRadIIRFilter
    "time constant of the first exponential to be deconvolved"
    cr::T
    "time constant of the second exponential to be deconvolved"
    cr2::U
    "the fraction faktor which the second exponential contributes"
    f::V
end

export InvSecondOrderCRFilter

function fltinstance(flt::InvSecondOrderCRFilter, fi::SamplingInfo)
    fltinstance(BiquadFilter(InvSecondOrderCRFilter(ustrip(NoUnits, flt.cr / step(fi.axis)), ustrip(NoUnits, flt.cr2 / step(fi.axis)), flt.f )), fi)
end

InverseFunctions.inverse(flt::InvSecondOrderCRFilter) = SecondOrderCRFilter(flt.cr, flt.cr2, flt.f)

function BiquadFilter(flt::InvSecondOrderCRFilter)
    a = exp(-1/float(flt.cr))
    b = exp(-1/float(flt.cr2))
    frac = float(flt.f)
    transfer_denom_1 = frac * b - frac * a - b - 1
    transfer_denom_2 = -(frac * b - frac * a - b)
    transfer_num_1 = -(a + b)
    transfer_num_2 = a * b
    BiquadFilter((1.0, transfer_num_1, transfer_num_2), (transfer_denom_1, transfer_denom_2))
end


"""
    struct RC_CR2Filter{T<:RealQuantity} <: AbstractRadIIRFilter

A RC-CR² shaping filter useful for determining pileup and trigger times.
The filter is computed using a matched z-transform to keep the poles/zeroes 
of the analog transfer function in the same location.

Constructors:

* ```$(FUNCTIONNAME)(fields...)```

Fields:

$(TYPEDFIELDS)

YAML Configuration Example
--------------------------

.. code-block:: yaml

    wf_RC_CR2:
      function: rc_cr2
      module: RadiationDetectorDSP
      args:
        - wf_bl
        - "300*ns"
        - wf_RC_CR2
"""
Base.@kwdef struct RC_CR2Filter{T<:RealQuantity} <: AbstractRadIIRFilter
    "RC-CR² time constant"
    tau::T
end

export RC_CR2Filter


struct RC_CR2FilterInstance{T} <: AbstractRadSigFilterInstance{LinearFiltering}
    a::T
    denom_2::T
    denom_3::T
    denom_4::T
    n::Int
end


function fltinstance(flt::RC_CR2Filter, fi::SamplingInfo)
    tau_norm = float(flt.tau / step(fi.axis))
    T = typeof(tau_norm)  # Get the numeric type (Float32 or Float64)
    a = exp(-1 / tau_norm)
    
    denom_2 = -3 * a
    denom_3 = 3 * a^2
    denom_4 = -(a^3)
    
    RC_CR2FilterInstance{T}(a, denom_2, denom_3, denom_4, _smpllen(fi))
end




Adapt.adapt_structure(to, flt::RC_CR2Filter) = flt


@inline function rdfilt!(Y::AbstractVector{T}, fi::RC_CR2FilterInstance{T}, X::AbstractVector{T}) where {T<:Real}
    # Check input validity
    if any(isnan, X) || length(X) <= 3
        fill!(Y, T(NaN))
        return Y
    end
    
    # Initialize first three samples
    Y[1] = X[1]
    Y[2] = X[2]
    Y[3] = X[3]
    
    # Use higher precision buffer to avoid float truncation
    w_tmp = zeros(Float64, 4)
    w_tmp[1] = Float64(X[1])
    w_tmp[2] = Float64(X[2])
    w_tmp[3] = Float64(X[3])
    
    a = Float64(fi.a)
    denom_1 = 1.0
    denom_2 = Float64(fi.denom_2)
    denom_3 = Float64(fi.denom_3)
    denom_4 = Float64(fi.denom_4)
    
    num_1 = 1.0
    num_2 = -2.0
    num_3 = 1.0
    
    @inbounds for i in 4:length(X)
        w_tmp[4] = (
            -denom_2 * w_tmp[3]
            - denom_3 * w_tmp[2]
            - denom_4 * w_tmp[1]
            + num_1 * Float64(X[i])
            + num_2 * Float64(X[i - 1])
            + num_3 * Float64(X[i - 2])
        ) / denom_1
        
        Y[i] = T(w_tmp[4])
        
        # Shuffle the buffers
        w_tmp[1] = w_tmp[2]
        w_tmp[2] = w_tmp[3]
        w_tmp[3] = w_tmp[4]
    end
    
    # Check output for NaNs
    if any(isnan, Y)
        fill!(Y, T(NaN))
    end
    
    Y
end


# RC_CR2FilterInstance methods
adapt_memlayout(::RC_CR2FilterInstance, ::GPU, A::AbstractArray{<:Number}) = _row_major(A)

function bc_rdfilt!(
    outputs::ArrayOfSimilarVectors{<:RealQuantity},
    fi::RC_CR2FilterInstance,
    inputs::ArrayOfSimilarVectors{<:RealQuantity}
)
    _ka_bc_rdfilt!(outputs, fi, inputs)
end

flt_output_smpltype(fi::RC_CR2FilterInstance) = flt_input_smpltype(fi)
flt_input_smpltype(fi::RC_CR2FilterInstance{T}) where T = T

flt_output_length(fi::RC_CR2FilterInstance) = flt_input_length(fi)
flt_input_length(fi::RC_CR2FilterInstance) = fi.n
flt_output_time_axis(fi::RC_CR2FilterInstance, time::AbstractVector{<:RealQuantity}) = time
