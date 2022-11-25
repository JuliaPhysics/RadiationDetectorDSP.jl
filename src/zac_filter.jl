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
    T <: RealQuantity
} <: AbstractRadFIRFilter
    "equivalent of shaping time (τₛ)"
    tau::T = 20,

    "length of flat top (FT)"
    toplen::T = 10

    "total length of the filter (L)"
    length::T = 100
end

export ZACChargeFilter

Adapt.adapt_structure(to, flt::ZACChargeFilter) = flt

function fltinstance(flt::ZACChargeFilter, fi::SamplingInfo)
    fltinstance(ConvolutionFilter(ZACChargeFilter(
        ustrip(NoUnits, flt.tau / step(fi.axis)),
        ustrip(NoUnits, flt.toplen / step(fi.axis)),
        round(ustrip(NoUnits, flt.length / step(fi.axis)))
    )), fi)
end

function ConvolutionFilter(flt::ZACChargeFilter)
    coeffs = zac_charge_filter_coeffs(flt.toplen, flt.tau, 1, flt.length)
    ConvolutionFilter(FFTConvolution(), coeffs)
end


"""
    zac_current_filter_integral(L, FT, τₛ)

Compute the integral of the unnormalised zac filter using its analytical
expression, where `L` is the length of the polynomial part, `FT` the length 
of the flat top part and `τₛ` the filter shaping time
"""
@inline zac_current_filter_integral(L, FT, τₛ) = 
    6/L^3 * (τₛ*(cosh(L/τₛ) - 1) + sinh(L/τₛ)*FT/2)

"""
    zac_current_filter_shape(t, τₛ, L, FT, A)
       
value of zac filter evaluated at `t` for the given set of parameters 
(see `zac_charge_filter_coeffs`)
"""
zac_current_filter_shape(t::T, τₛ::T, L::T, FT::T, A::T
) where {T<:AbstractFloat} = begin
    if zero(T) < t < L 
        sinh(t / τₛ) + A*t*(t - L)
    elseif L <= t <= L + FT
        sinh(L / τₛ)
    elseif L + FT < t < 2L + FT
        sinh((2L + FT - t)/τₛ) + A*(L + FT - t)*(2L + FT - t)
    else
        zero(T)
    end
end

"""
    zac_charge_filter_coeffs(FT, τₛ, Δt, T)

return a vector representing the zac filter applicaible on a charge 
signal, where `FT` is the length of the flat top, `τₛ` the filter 
shaping time, `Δt` the sampling time and `T` the total length of the 
filter (see Eur. Phys. J. C (2015) 75:255).
"""
zac_charge_filter_coeffs(FT, τₛ, Δt, T) = begin
    L = (T - FT) / 2
    A = zac_current_filter_integral(L, FT, τₛ)
    FT_t, tau_t, Δt_t, T_t, A_t, L_t  = promote(FT, τₛ, Δt, T, A, L)
    U = typeof(FT_t)
    t = Δt_t:Δt_t:(T_t - Δt_t)
    zac_diff = zeros(U, length(t)-1)
    for i in eachindex(zac_diff)
        fᵢ = zac_current_filter_shape(t[i], tau_t, L_t, FT_t, A_t)
        fᵢ₊₁ = zac_current_filter_shape(t[i+1], tau_t, L_t, FT_t, A_t)
        zac_diff[i] = (fᵢ₊₁ - fᵢ)
    end
    return zac_diff
end

@inline tons(u, dt) = ceil(u / dt) |> Int

"""
    create_zac_filter(Nₜ, FTₜ, τₜ, Δt)

create the zac filter according to Eur. Phys. J. C (2015) 75:255], where 
Nₜ is the length of the filter, FTₜ the length of the flat top, τₜ the the 
exponential decay factor and Δt the sampling time
"""
function create_zac_filter(Nₜ, FTₜ, τₜ, Δt)
    # convert parameters to integers respecting whatever resolution is 
    # given by Δt
    N, τₛ, FT = tons.((Nₜ, τₜ, FTₜ), Δt)

    L = ((N - FT) % 2 == 0) ? (N - FT)÷2 : (N - (FT+=1))÷2
    CUSP = Array{Float64, 1}(undef, N)
    Poly = zeros(Int, N)

    # build the cusp filter (sinh part) and polynomial part
    for i in 1:L
        CUSP[i] = CUSP[N-i+1] = sinh(i/τₛ)
        Poly[i] = Poly[N-i+1] = i*(i - L)
    end
    # build the flat top part
    C = sinh(L/τₛ)
    for i in 1:FT
        CUSP[L+i] = C
    end

    ∫Poly = mean(Poly)
    ∫CUSP = mean(CUSP)

    # build the normalized ZAC filter
    CUSP .+ (Poly .* (-∫CUSP/∫Poly))
end