# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

module RadiationDetectorDSPCUDAExt

using CUDA
using RadiationDetectorDSP


function RadiationDetectorDSP._nonlazy_transpose(A::CuArray{T,2}) where {T<:Number}
    A_T = similar(A, reverse(size(A)))
    # CUBLAS.geam! is extremely fast:
    CUDA.CUBLAS.geam!('T', 'T', one(T), A, zero(T), A, A_T)
    return A_T
end


end # module RadiationDetectorDSPCUDAExt
