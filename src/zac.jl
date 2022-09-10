export zac_filter_coefficients

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