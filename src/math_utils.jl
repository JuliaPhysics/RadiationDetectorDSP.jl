# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


_floattype(::Type{T}) where T = float(T)
_floattype(::Type{Int8}) = Float32
_floattype(::Type{UInt8}) = Float32
_floattype(::Type{Int16}) = Float32
_floattype(::Type{UInt16}) = Float32
_floattype(::Type{Int32}) = Float32
_floattype(::Type{UInt32}) = Float32


function _smoothstep(x::Real)
    @fastmath begin
        xc = ifelse(x < zero(x), zero(x), ifelse(x > one(x), one(x), x))
        xc = clamp(x, 0, 1)
        x_2 = xc * xc 
        3 * x_2 - 2 * x_2 * xc
    end
end
