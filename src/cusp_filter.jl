# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    struct CUSPChargeFilter <: AbstractRadFIRFilter

CUSP filter.

For the definition the filter and a discussion of the filter properties, see

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
@with_kw struct CUSPChargeFilter{
    T <: RealQuantity, U <: RealQuantity, V <: RealQuantity, W <: Real
} <: AbstractRadFIRFilter
    "equivalent of shaping time (τₛ)"
    sigma::U = 450

    "length of flat top (FT)"
    toplen::T = 10

    "decay constant of the exponential"
    tau::V = 20

    "total length of the filter (L)"
    length::T = 100

    "scaling factor"
    beta::W = 100.
end

export CUSPChargeFilter

Adapt.adapt_structure(to, flt::CUSPChargeFilter) = flt

function fltinstance(flt::CUSPChargeFilter, fi::SamplingInfo)
    fltinstance(ConvolutionFilter(CUSPChargeFilter(
        ustrip(NoUnits, flt.sigma / step(fi.axis)),
        round(Int, ustrip(NoUnits, flt.toplen / step(fi.axis))),
        ustrip(NoUnits, flt.tau / step(fi.axis)),
        round(Int, ustrip(NoUnits, flt.length / step(fi.axis))),
        flt.beta
    )), fi)
end

function ConvolutionFilter(flt::CUSPChargeFilter)
    coeffs = cusp_charge_filter_coeffs(
        flt.length, flt.sigma, flt.toplen, flt.tau, flt.beta)
    ConvolutionFilter(FFTConvolution(), coeffs)
end

"""
    cusp_charge_filter_coeffs(N::Int, sigma::U, FT::Int, tau::V, beta::W
    ) where {U, V, W <: AbstractFloat}

return a vector representing the cusp filter applicaible on a charge 
signal, where `N` is the total length of the filter, `FT` the length of 
the flat top, `sigma` the filter shaping time,`tau` the decay constant 
and `a` the scaling factor.
"""
function cusp_charge_filter_coeffs(N::Int, sigma::U, FT::Int, tau::V, beta::W
    ) where {U, V, W <: AbstractFloat}

    T = promote_type(U, V, W)
    L::Int = ((N - FT) % 2 == 0) ? (N - FT)÷2 : (N - (FT+=1))÷2
    FF = Vector{T}(undef, N-1)
    
    α::T = -exp(-1.0/tau)
    C::T = sinh(L/sigma)
    β::T = beta/N/C         # scaling factor
    Δ::T = (α + 1)*beta/N
    
    # directly compute CUSP convolved with inverse response function
    # FF[i] = α * CUSP[i] + CUSP[i+1]
    for i in Base.OneTo(L)
        FF[i] = β*(α*sinh((i - 1)/sigma) + sinh(i/sigma))
        FF[end-i+1] = β*(α*sinh(i/sigma) + sinh((i - 1)/sigma))
    end
    for i in Base.OneTo(FT-1)
        FF[L+i] = Δ
    end
    FF
end