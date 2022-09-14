export zac_filter_coefficients, create_zac_filter

"""
    zac_filter_coefficients(s, τₛ, A, L, FT)

ZAC filter function as defined in Eq. (7) in [Eur. Phys. J. C (2015) 75:255]
"""
@inline zac_filter_coefficients(t, τₛ, A, L, FT) = begin
    if 0 <= t < L
        sinh(t / τₛ) + A*t*(t - L)
    elseif L <= t < L + FT
        sinh(L / τₛ)
    elseif L + FT <= t <= 2L + FT
        sinh((2L + FT - t)/τₛ) + A*((L + FT - t)*(2L + FT - t))
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
    @debug "Number of bins: $L"
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
    # build the normalized ZAC filter
    CUSP .+ (Poly .* (-∫CUSP/∫Poly))
end