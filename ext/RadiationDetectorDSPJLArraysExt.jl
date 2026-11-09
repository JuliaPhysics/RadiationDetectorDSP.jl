# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

module RadiationDetectorDSPJLArraysExt

using JLArrays

import RadiationDetectorDSP
import KernelAbstractions


# Before KernelAbstractions v0.10 (recognizable by its POCL backend), the JLArrays backend
# is incomplete, but the native KernelAbstractions CPU backend can run kernels on JLArrays:
@static if !isdefined(KernelAbstractions, :POCLBackend)
    RadiationDetectorDSP._ka_get_backend(X::Union{JLArray,SubArray{T,N,<:JLArray} where {T,N}}) = KernelAbstractions.CPU()
end


end # module RadiationDetectorDSPJLArraysExt
