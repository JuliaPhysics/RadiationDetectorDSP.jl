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
        xc = clamp(x, zero(x), one(x))
        xc * xc * (3 - 2 * xc)
    end
end
