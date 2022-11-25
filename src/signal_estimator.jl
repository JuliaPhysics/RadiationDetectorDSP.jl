# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    abstract type DNIMethod

Abstract type for denoising and interpolation methods.
"""    
abstract type DNIMethod end
export DNIMethod



"""
    struct PolynomialDNI <: DNIMethod

Polynomial denoising and interpolation method.

Operates in a similar way as a
[Savitzky-Golay filter](https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter),
but interpolates as well.

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
@with_kw struct PolynomialDNI{
    T <: RealQuantity
} <: DNIMethod
    "polynomial degree"
    degree::Int = 2,

    "length"
    length::T
end

export PolynomialDNI


"""
    struct SignalEstimator <: Function

Estimates a signal at a given position `x`.

Usage:

```julia
(f::SamplesOrWaveform)(input::RDWaveform, x::RealQuantity)
```

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
@with_kw struct SignalEstimator{
    DNI <: DNIMethod
} <: Function
    "denoising and interpolation method"
    dni::DNI
end

export SignalEstimator


function (f::SignalEstimator{<:PolynomialDNI})(input::SamplesOrWaveform, x::RealQuantity)
    X_axis, Y = _get_axis_and_signal(input)
    degree = f.dni.degree
    n = round(Int, ustrip(NoUnits, f.dni.length / step(X_axis)))
    @argcheck n > degree
    lsqfitatpos(degree, n, X_axis, Y, x)
end


function Base.Broadcast.broadcasted(f::SignalEstimator{<:PolynomialDNI}, bc_inputs, bc_xs)
    inputs = Base.materialize(bc_inputs)
    X_axis, Ys = _get_axis_and_signals(inputs)
    X_pos = Base.materialize(bc_xs)
    degree = f.dni.degree
    # ToDo: Remove code duplication
    n = round(Int, ustrip(NoUnits, f.dni.length / step(X_axis)))
    _sigest_poly_bc(degree, n, X_axis, Ys, X_pos)
end

function _sigest_poly_bc(degree::Integer, n::Integer, X_axis::AbstractArray{<:RealQuantity}, Ys::_BC_RQ_Arrays, X_pos::_BC_RQs)
    lsqfitatpos.(degree, n, Ref(X_axis), Ys, X_pos)
end

function _sigest_poly_bc(degree::Integer, n::Integer, X_axis::AbstractRange{<:RealQuantity}, Ys::_BC_RQ_AosAs, X_pos::_BC_RQs)
    result_sz = Base.Broadcast.broadcast_shape(size(Ys), size(X_pos))
    Ys_flat, _ = _kbc_flatview(Ys)
    Y_est = similar(Ys_flat, eltype(Ys_flat), result_sz)
    bc_lsqfitat!(Y_est, degree, n, X_axis, Ys, X_pos)
end
