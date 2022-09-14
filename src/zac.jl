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

"""

    ZACfilter similar to cpp script
"""
@inline tons(u, dt) = ceil(u / dt) |> Int

create_zac_filter(Nₜ, FTₜ, τₜ, Δt) = begin
    # convert parameters to integers respecting nano seconds resolution
    N, τₛ, FT = tons.((Nₜ, τₜ, FTₜ), Δt)

    # define length of filter
    L = ((N - FT) % 2 == 0) ? L = (N - FT)÷2 : (FT += 1; true) && (N - FT)÷2
    CUSP = Array{Float64, 1}(undef, N)
    Poly = zeros(Int, N)

    # normalization constant, such that the ZAC filter is equal to one at the 
    # center
    C = 1/sinh(L/τₛ)

    # build the cusp filter (sinh part) and polynomial part
    for i in 1:L
        y = sinh(i/τₛ)*C
        CUSP[i] = CUSP[N-i+1] = y
        Poly[i] = Poly[N - i + 1] = i*(i - L)
    end
    for i in L+1:L+FT
        CUSP[i] = 1.
    end

    ∫Poly = sum(Poly)
    ∫CUSP = sum(CUSP)

    # build the ZAC filter
    return CUSP .+ (Poly .* (-∫CUSP/∫Poly))
end