export zac_filter_coefficients, create_zac_filter, 
       zac_diff_filter_coefficients, create_zac_diff_filter

@inline _A(L, FT, τₛ) = 6/L^3 * (τₛ*(cosh(L/τₛ) - 1) + sinh(L/τₛ)*FT/2)

"""
    zac_filter_coefficients(s, τₛ, L, FT)
       
ZAC filter function as defined in Eq. (7) in [Eur. Phys. J. C (2015) 75:255]
"""
@inline zac_filter_coefficients(t, τₛ, L, FT) = begin
    if 0 <= t < L
        sinh(t / τₛ) + _A(L, FT, τₛ)*t*(t - L)
    elseif L <= t <= L + FT
        sinh(L / τₛ)
    elseif L + FT <= t <= 2L + FT
        sinh((2L + FT - t)/τₛ) + _A(L, FT, τₛ)*(L + FT - t)*(2L + FT - t)
    else
        0
    end
end

@inline zac_diff_filter_coefficients(t, Δt, τ, τₛ, L, FT) = begin
    a = exp(-1/τ)
    if 0 <= t < L
        A = _A(L, FT, τₛ)
        # @debug "A = $A"
        sinh((t+Δt)/τₛ) - a*sinh(t/τₛ) + A*((t + Δt)*Δt + (t + Δt - a*t)*(t - L))
    elseif L <= t <= L+FT
        (1-a)*sinh(L/τₛ)
    elseif L + FT < t <= 2L + FT
        A = _A(L, FT, τₛ)
        # @debug "A = $A"
        -a*(sinh((t+Δt)/τₛ) - sinh(t/τₛ)/a + A*((t + Δt)*Δt + (t + Δt - t/a)*(t - L)))
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

    # define length of filter
    L = ((N - FT) % 2 == 0) ? (N - FT)÷2 : (N - (FT+=1))÷2
    @debug "partial length: $L"
    @debug "Total length of Filter: $N"
    @debug "exponential decay: $(τₛ)"
    @debug "Flat top: $FT"
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

    ∫Poly = sum(Poly)
    ∫CUSP = sum(CUSP)
    @debug "∫Poly = $(∫Poly)"
    @debug "∫CUSP = $(∫CUSP)"
    @debug "A = $(-∫CUSP/∫Poly)"
    # build the normalized ZAC filter
    CUSP .+ (Poly .* (-∫CUSP/∫Poly))
end

function create_zac_diff_filter(Nₜ, FTₜ, τₜ, τ₂ₜ, Δt)
    # convert parameters to integers respecting whatever resolution is 
    # given by Δt
    N, τₛ, FT, τ = tons.((Nₜ, τₜ, FTₜ, τ₂ₜ), Δt)

    # define length of filter
    L = ((N - FT) % 2 == 0) ? (N - FT)÷2 : (N - (FT+=1))÷2
    @debug "partial length: $L"
    @debug "Total length of Filter: $N"
    @debug "exponential decay: $(τₛ)"
    @debug "Flat top: $FT"
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

    ∫Poly = sum(Poly)
    ∫CUSP = sum(CUSP)
    @debug "∫Poly = $(∫Poly)"
    @debug "∫CUSP = $(∫CUSP)"
    @debug "A = $(-∫CUSP/∫Poly)"
    # build the normalized ZAC filter
    CUSP .+ (Poly .* (-∫CUSP/∫Poly))
end