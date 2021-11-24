# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

const RealOrSIMD{T<:Real} = Union{T,<:Vec{N,<:T} where N}
