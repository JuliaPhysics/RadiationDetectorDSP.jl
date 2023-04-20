# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

module RadiationDetectorDSPJLArraysExt

isdefined(Base, :get_extension) ? (using JLArrays) : (using ..JLArrays)

import RadiationDetectorDSP
import KernelAbstractions


RadiationDetectorDSP._ka_get_backend(X::Union{JLArray,SubArray{T,N,<:JLArray} where {T,N}}) = KernelAbstractions.CPU()


end # module RadiationDetectorDSPJLArraysExt
