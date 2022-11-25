# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    struct SavitzkyGolayFilter <: AbstractRadFIRFilter

A [Savitzky-Golay filter](https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter).

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct SavitzkyGolayFilter{T<:RealQuantity} <: AbstractRadFIRFilter
    "filter length"
    length::T

    "Polynomial defgree"
    degree::Int

    "n-th derivative (0 for no derivative)"
    derivative::Int = 0
end

export SavitzkyGolayFilter

Adapt.adapt_structure(to, flt::SavitzkyGolayFilter) = flt

function fltinstance(flt::SavitzkyGolayFilter, fi::SamplingInfo)
    # ToDo: Round up and down to nearest odd length and use weighted sum of coeffs?
    len_orig = round(Int, ustrip(NoUnits, flt.length / step(fi.axis)))
    len_odd = 2 * div(len_orig, 2) + 1
    fltinstance(ConvolutionFilter(SavitzkyGolayFilter(
        len_odd,
        flt.degree,
        flt.derivative
    )), fi)
end

function ConvolutionFilter(flt::SavitzkyGolayFilter)
    rmin = - div(Int(flt.length), 2)
    rmax = flt.length + rmin - 1
    coeffs = RadiationDetectorDSP.sg_filter_coeffs(rmin:rmax, flt.degree, flt.derivative)
    ConvolutionFilter(DirectConvolution(), coeffs, rmin)
end


using LinearAlgebra

function sg_filter_coeffs(x_range::AbstractUnitRange{<:Integer}, degree::Integer, deriv::Integer, delta::Real)
    J = _vandermonde(x_range, degree)
    y = _onehot_vec(degree + 1, deriv + 1, factorial(deriv) / delta ^ deriv)
    J' \ y
end

function sg_filter_coeffs(x_range::AbstractUnitRange{<:Integer}, degree::Integer, deriv::Integer)
    J = _vandermonde(x_range, degree)
    y = _onehot_vec(degree + 1, deriv + 1, factorial(deriv))
    J' \ y
end
