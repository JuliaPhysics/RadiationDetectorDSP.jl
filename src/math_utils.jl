# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


function _smoothstep(x::Real)
    @fastmath begin
        xc = ifelse(x < zero(x), zero(x), ifelse(x > one(x), one(x), x))
        xc = clamp(x, 0, 1)
        x_2 = xc * xc 
        3 * x_2 - 2 * x_2 * xc
    end
end
