export zac_current_filter_integral, zac_current_filter_shape,
       zac_charge_filter_coeffs

"""
    zac_current_filter_integral(L, FT, τₛ)

compute the integral of the unnormalised zac filter using its analytical
expression, where L is the length of the polynomial part, FT the length 
of the flat top part and τₛ the 
"""
@inline zac_current_filter_integral(L, FT, τₛ) = 
    6/L^3 * (τₛ*(cosh(L/τₛ) - 1) + sinh(L/τₛ)*FT/2)

"""
    zac_current_filter_shape(t, τₛ, L, FT, A)
       
value of zac filter evaluated at `t` for the given set of parameters
(see Eur. Phys. J. C (2015) 75:255)
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
    zac_charge_filter_coeffs(FT, tau, Δt, L)

return a vector representing the zac filter applicaible on a charge 
signal.
"""
zac_charge_filter_coeffs(FT, tau, Δt, L) = begin
    T = promote_type(typeof.((FT, tau, Δt, L))...)
    A = zac_current_filter_integral(L, FT, τₛ)
    t = zero(T):Δt:L
    zac_diff = zeros(T, length(t)-1)
    for i in eachindex(zac_diff)
        nothing
    end
    return zac_diff
end

"""
    zac_diff_filter_coefficients(t, Δt, τ, τₛ, L, FT)

filter coefficients equivalent to a ZAC filter convolved with a 
differential filter. 
"""
@inline zac_diff_filter_coefficients(t, Δt, τ, τₛ, L, FT) = begin
    a = exp(-Δt/τ)
    if 0 <= t < L
        A = _A(L, FT, τₛ)
        sinh((t+Δt)/τₛ) - a*sinh(t/τₛ) 
        + A*((t + Δt)*Δt + (t + Δt - a*t)*(t - L))
    elseif L <= t <= L+FT
        (1-a)*sinh(L/τₛ)
    elseif L + FT < t <= 2L + FT
        A = _A(L, FT, τₛ)
        t̃ = 2L + FT - t
        -a*(sinh((t̃+Δt)/τₛ) - sinh(t̃/τₛ)/a 
        + A*((t̃ + Δt)*Δt + (t̃ + Δt - t̃/a)*(t̃ - L)))
    else
        0
    end
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