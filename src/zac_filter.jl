# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    struct ZACChargeFilter <: AbstractRadFIRFilter

Zero area cusp (ZAC) filter.

For the definition the filter and a discussion of the filter properties, see
["Improvement of the energy resolution via an optimized digital signal processing in GERDA Phase I", Eur. Phys. J. C 75, 255 (2015)](https://doi.org/10.1140/epjc/s10052-015-3409-6).

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
@with_kw struct ZACChargeFilter{
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

export ZACChargeFilter

Adapt.adapt_structure(to, flt::ZACChargeFilter) = flt

function fltinstance(flt::ZACChargeFilter, fi::SamplingInfo)
    fltinstance(ConvolutionFilter(ZACChargeFilter(
        ustrip(NoUnits, flt.sigma / step(fi.axis)),
        round(Int, ustrip(NoUnits, flt.toplen / step(fi.axis))),
        ustrip(NoUnits, flt.tau / step(fi.axis)),
        round(Int, ustrip(NoUnits, flt.length / step(fi.axis))),
        flt.beta
    )), fi)
end

function ConvolutionFilter(flt::ZACChargeFilter)
    coeffs = zac_charge_filter_coeffs(
        flt.length, flt.sigma, flt.toplen, flt.tau, flt.beta)
    ConvolutionFilter(FFTConvolution(), coeffs)
end

"""
    zac_charge_filter_coeffs(
        N::Int, sigma::V, FT::Int, tau::T, beta::U
        ) where {V, T, U <: AbstractFloat}

return a vector representing the zac filter applicaible on a charge 
signal, where `N` is the total length of the filter, `FT` the length of 
the flat top, `sigma` the filter shaping time, `tau` the decay constant
and beta an additional scaling factor.
(see Eur. Phys. J. C (2015) 75:255).
"""
function zac_charge_filter_coeffs(N::Int, sigma::U, FT::Int, tau::V, beta::W
    ) where {U, V, W <: AbstractFloat}

    T = promote_type(U, V, W)
    L::Int = ((N - FT) % 2 == 0) ? (N - FT)÷2 : (N - (FT+=1))÷2
    FF = Vector{T}(undef, N-1)
    
    α::T = -exp(-1.0/tau)
    C::T = sinh(L/sigma)
    β::T = beta/N             # scaling factor
    Δ::T = (α + 1)*β
    # sum of cusp filter
    A::T = FT + 
        1.0/((C / (1 - exp(L / sigma))) - (C / (1 - exp((L-1) /sigma))))
    # sum of polynomial part
    A /= (L^3 - L)/3
    
    # directly compute ZAC convolved with inverse response function
    # FF[i] = α * ZAC[i] + ZAC[i+1]
    for i in Base.OneTo(L)
        FF[i] = β*α*(sinh((i - 1)/sigma)/C + A*(i - 1)*((i - 1) - L)) 
        FF[i] += β*(sinh(i/sigma)/C + A*i*(i - L))
        FF[end-i+1] = β*α*(sinh(i/sigma)/C + A*i*(i - L)) 
        FF[end-i+1] += β*(sinh((i - 1)/sigma)/C + A*(i - 1)*((i - 1) - L))
    end
    for i in Base.OneTo(FT-1)
        FF[L+i] = Δ
    end
    FF
end