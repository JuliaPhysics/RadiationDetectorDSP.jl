# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    struct Gauss1DFilter <: AbstractRadFIRFilter

One dimensional gaussian filter defined as:
f(x) = beta * exp(-0.5*(x/sigma)^2) / length

where x is in the interval [-alpha*sigma, alpha*sigma]

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
@with_kw struct Gauss1DFilter{
    T <: RealQuantity
} <: AbstractRadFIRFilter
    "standard deviation"
    sigma::T = 1.

    "the amount of standard deviations to cover in the gaussian window"
    alpha::T = 1.

    "scaling factor"
    beta::T = 1.

    "total length of the filter"
    length::T = 100.
end

export Gauss1DFilter

Adapt.adapt_structure(to, flt::Gauss1DFilter) = flt

function fltinstance(flt::Gauss1DFilter, fi::SamplingInfo)
    fltinstance(ConvolutionFilter(Gauss1DFilter()), fi)
end

function ConvolutionFilter(flt::Gauss1DFilter)
    coeffs = gaussian_coeffs(
        flt.length, flt.sigma, flt.alpha, flt.beta)
    ConvolutionFilter(FFTConvolution(), coeffs)
end

"""
    gaussian_coeffs(N::Int, alpha::T, sigma::T, beta::T) where {T}

compute a gaussian kernel, where `N` is the total length of the kernel, 
`alpha` the amount of standard deviations to cover, `sigma` the standard 
deviation and `beta` the total scaling factor.
"""
function gaussian_coeffs(N::Int, alpha::T, sigma::T, beta::T) where {T}
    y = Vector{Float64}(undef, N)
    xᵢ = -sigma*alpha
    Δx = abs(2*xᵢ)/(N-1)
    for i=Base.OneTo(N)
        y[i] = beta*exp(-0.5*(xᵢ/alpha)^2) / N
        xᵢ += Δx
    end
    y
end