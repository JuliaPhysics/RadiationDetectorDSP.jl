export zac_filter!, zac

"""
    zac_filter(s, τₛ, A, L, FT)

Apply the ZAC filter (Eur. Phys. J. C (2015) 75:255) to the array s and 
save the result in s.

"""
zac_filter!(s::AbstractArray{T}, τₛ::T, A::T, L::T, FT::T) = 
    s .= zac.(s, τₛ, A, L, FT)

@inline zac(t::T, τₛ::T, A::T, L::T, FT::T) where T = begin
    if 0 <= t < L
        sinh(t / τₛ) + A*t*(t - L)
    elseif L <= t < L + FT
        sinh(L / τₛ)
    elseif L + FT <= t <= 2L + FT
        sinh((2L + FT - t)/τₛ) + A*((L + FT - t)*(2L + FT - t))
    end
end