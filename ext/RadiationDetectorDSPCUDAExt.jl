# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

module RadiationDetectorDSPCUDAExt

isdefined(Base, :get_extension) ? (using CUDA) : (using ..CUDA)
using RadiationDetectorDSP
import Adapt


function RadiationDetectorDSP._nonlazy_transpose(A::CuArray{T,2}) where {T<:Number}
    A_T = similar(A, reverse(size(A)))
    # CUBLAS.geam! is extremely fast:
    CUDA.CUBLAS.geam!('T', 'T', one(T), A, zero(T), A, A_T)
    return A_T
end


RadiationDetectorDSP.deviceof_impl(@nospecialize(TypeHistory::Type), A::CuArray) = CUDA.device(A)
RadiationDetectorDSP.deviceof_impl(@nospecialize(TypeHistory::Type), A::CUDA.CUDA.CUSPARSE.AbstractCuSparseArray) = CUDA.device(A.nzVal)

RadiationDetectorDSP.deviceof_impl(@nospecialize(TypeHistory::Type), e::CUDA.CuEvent) = device(e.ctx)


struct _CuDeviceAdaptor
    devhandle::Int32
end

RadiationDetectorDSP.device_adaptor(dev::CuDevice) = _CuDeviceAdaptor(dev.handle)

function Adapt.adapt_storage(adaptor::_CuDeviceAdaptor, x)
    oldhandle = CUDA.device().handle
    try
        oldhandle != adaptor.devhandle && CUDA.device!(adaptor.devhandle)
        Adapt.adapt(CuArray, x)
    finally
        oldhandle != adaptor.devhandle && CUDA.device!(oldhandle)
    end
end


device_total_memory(dev::CuDevice) = CUDA.totalmem(dev)
device_free_memory(dev::CuDevice) = unsigned(CUDA.device!(CUDA.available_memory, dev))


end # module RadiationDetectorDSPCUDAExt
