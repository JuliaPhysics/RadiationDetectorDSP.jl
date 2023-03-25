function deviceof end
function device_adaptor end
function device_total_memory end
function device_free_memory end


struct ComputingDeviceIndependent end

struct UnknownComputeDeviceOf{T}
    x::T
end


struct MixedComputeSystem end


struct CPUDevice end

struct _CPUDeviceAdaptor end

device_adaptor(dev::CPUDevice) = _CPUDeviceAdaptor()

Adapt.adapt_storage(adaptor::_CPUDeviceAdaptor, x::AbstractArray) = Adapt.adapt(Array, x)

device_total_memory(dev::CPUDevice) = Sys.total_memory()
device_free_memory(dev::CPUDevice) = Sys.free_memory()



function merge_compute_devices end

merge_compute_devices() = ComputingDeviceIndependent()

@inline function merge_compute_devices(a, b, c, ds::Vararg{Any,N}) where N
    a_b = merge_compute_devices(a,b)
    return merge_compute_devices(a_b, c, ds...)
end

@inline merge_compute_devices(a::UnknownComputeDeviceOf, b::UnknownComputeDeviceOf) = a
@inline merge_compute_devices(a::UnknownComputeDeviceOf, b::Any) = a
@inline merge_compute_devices(a::Any, b::UnknownComputeDeviceOf) = b

@inline function merge_compute_devices(a, b)
    return (a === b) ? a : compute_device_mergeresult(
        compute_device_mergerule(a, b),
        compute_device_mergerule(b, a),
    )
end

struct NoCDeviceMergeRule end

@inline compute_device_mergerule(a::Any, b::Any) = NoCDeviceMergeRule()
@inline compute_device_mergerule(a::UnknownComputeDeviceOf, b::Any) = a
@inline compute_device_mergerule(a::UnknownComputeDeviceOf, b::UnknownComputeDeviceOf) = a
@inline compute_device_mergerule(a::ComputingDeviceIndependent, b::Any) = b

@inline compute_device_mergeresult(a_b::NoCDeviceMergeRule, b_a::NoCDeviceMergeRule) = MixedComputeSystem()
@inline compute_device_mergeresult(a_b, b_a::NoCDeviceMergeRule) = a_b
@inline compute_device_mergeresult(a_b::NoCDeviceMergeRule, b_a) = b_a
@inline compute_device_mergeresult(a_b, b_a) = a_b === b_a ? a_b : MixedComputeSystem()


deviceof(x) = deviceof_impl(Union{}, x)


function deviceof_impl end


@inline deviceof_impl(@nospecialize(TypeHistory::Type), ::Array) = CPUDevice()

# Guard against object reference loops:
@inline deviceof_impl(::Type{TypeHistory}, x::T) where {TypeHistory,T<:TypeHistory} = begin
    UnknownComputeDeviceOf(x) 
end

@generated function deviceof_impl(::Type{TypeHistory}, x) where TypeHistory
    if isbitstype(x)
        :(ComputingDeviceIndependent())
    else
        NewTypeHistory = Union{TypeHistory, x}
        impl = :(begin dev_0 = ComputingDeviceIndependent() end)
        append!(impl.args, [:($(Symbol(:dev_, i)) = merge_compute_devices(deviceof_impl($NewTypeHistory, getfield(x, $i)), $(Symbol(:dev_, i-1)))) for i in 1:fieldcount(x)])
        push!(impl.args, :(return $(Symbol(:dev_, fieldcount(x)))))
        impl
    end
end
