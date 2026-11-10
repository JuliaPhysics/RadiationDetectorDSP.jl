# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    abstract type FilteringType

Intended as a type parameter to designate the behavior of a filter
as linear or nonlinear.

Subtypes are [`LinearFiltering`](@ref) and [`NonlinearFiltering`](@ref).
"""
abstract type FilteringType end


"""
    abstract type LinearFiltering <: FilteringType

When used as a type parameter value, marks linear behavior of a filter.
"""
abstract type LinearFiltering <: FilteringType end
export LinearFiltering


"""
    abstract type NonlinearFiltering <: FilteringType

When used as a type parameter value, marks non linear behavior of a filter.
"""
abstract type NonlinearFiltering <: FilteringType end
export NonlinearFiltering



"""
    abstract type AbstractRadSigFilter{FT<:FilteringType} <: Function

Abstract type for signal filters.

Filters are callable as `(flt::AbstractRadSigFilter)(input)` and come with
specialized broadcasting.

Subtypes of `AbstractRadSigFilter` must implement

```
fltinstance(flt::AbstractRadSigFilter, si::SamplingInfo)::[`AbstractRadSigFilterInstance`](@ref)
```

Invertible filters should also implement

* `InverseFunctions.inverse(flt::SomeFilter)`

Note that while a filter may have an inverse, it may, depending on the
filter parameters, be very unstable in the presence of additional noise.
Filters with a high-pass characteristic pass high-frequency noise, so their
inverses pass such noise as well without amplifying it (substantially).
Filters with a low-pass characteristic, on the other hand, attenuate
high-frequency noise, so their inverses amplify such noise and are
typically not useful to deconvolve signals in practical applications.
"""
abstract type AbstractRadSigFilter{FT<:FilteringType} <: Function end
export AbstractRadSigFilter

(flt::AbstractRadSigFilter)(input) = rdfilt(fltinstance(flt, smplinfo(input)), input)

# ToDo: Support mutating broadcasts (needs a bc_rdfilt! that takes chained filters into account)
function Base.Broadcast.broadcasted(flt::AbstractRadSigFilter, inputs)
    X = Base.materialize(inputs)
    fi = fltinstance(flt, elsmplinfo(X))
    bc_rdfilt(fi, X)
end



"""
    abstract type AbstractRadSigFilterInstance{FT<:FilteringType}

Abstract type for signal filter instances. Filter instances are specialized to
a specific length and numerical type of input and output.

Filter instances are callable as `(fi::SomeFilterInstance)(input)` and come
with specialized broadcasting.

Subtypes of `AbstractRadSigFilterInstance` must implement

* `RadiationDetectorDSP.rdfilt!(output, fi::SomeFilterInstance, input)`

* `RadiationDetectorDSP.flt_output_smpltype(fi::SomeFilterInstance)`

* `RadiationDetectorDSP.flt_input_smpltype(fi::SomeFilterInstance)`

* `RadiationDetectorDSP.flt_output_length(fi::SomeFilterInstance)`

* `RadiationDetectorDSP.flt_input_length(fi::SomeFilterInstance)`

* `RadiationDetectorDSP.flt_output_time_axis(fi::SomeFilterInstance, time::AbstractVector{<:RealQuantity})`

Invertible filter instances should implement

* `InverseFunctions.inverse(fi::SomeFilterInstance)`

Default methods are implemented for

* `RadiationDetectorDSP.rdfilt(fi::AbstractRadSigFilterInstance, x::AbstractSamples)`
* `RadiationDetectorDSP.rdfilt(fi::AbstractRadSigFilterInstance, wf::RDWaveform)`
* `RadiationDetectorDSP.bc_rdfilt(fi::AbstractRadSigFilterInstance, inputs)`

The default methods that operate on `RadiationDetectorSignals.RDWaveform`s require
[`RadiationDetectorDSP.flt_output_time_axis`](@ref).
"""
abstract type AbstractRadSigFilterInstance{FT<:FilteringType} end
export AbstractRadSigFilterInstance

# # ToDo: Should filter instances be callable?
#(fi::AbstractRadSigFilterInstance)(input) = rdfilt(fltinstance(flt, smplinfo(input)), input)
## ToDo: Support mutating broadcasts (needs a bc_rdfilt! that takes chained filters into account)
#Base.Broadcast.broadcasted(fi::AbstractRadSigFilterInstance, inputs) = bc_rdfilt(fi, Base.materialize(inputs))


"""
    fltinstance(flt::AbstractRadSigFilter, si::SamplingInfo)::AbstractRadSigFilterInstance

Create a filter instance of the filter `flt`, specialized for the given
input, resp. input characteristics.
"""
function fltinstance end
export fltinstance

fltinstance(flt::AbstractRadSigFilter, si::SamplingInfo{T}) where T = throw(ArgumentError("fltinstance not defined for type $(nameof(typeof(flt)))"))

#=
ToDo: Do we want these convenience methods?
fltinstance(flt::AbstractRadSigFilter, input::AbstractSamples) = fltinstance(flt, smplinfo(input))
fltinstance(flt::AbstractRadSigFilter, input::RDWaveform) = fltinstance(flt, smplinfo(input))
=#


"""
    rdfilt!(output, fi::AbstractRadSigFilterInstance, input)

Apply filter `flt` or filter instance `fi` to signal `input` and store the
filtered signal in `output`. Return `output`.
"""
function rdfilt! end
export rdfilt!

# ToDo: Do we want this?
# rdfilt!(output, flt::AbstractRadSigFilter, input) = rdfilt!(output, fltinstance(flt, smplinfo(input)), input)


"""
    rdfilt(fi::AbstractRadSigFilterInstance, input)

Apply filter instance `fi` to signal `input`, return the filtered signal.

Returns `output`.
"""
function rdfilt end
export rdfilt

# ToDo: Do we want this?
# rdfilt(flt::AbstractRadSigFilter, input) = rdfilt(fltinstance(flt, smplinfo(input)), input)

function rdfilt(fi::AbstractRadSigFilterInstance, input::AbstractSamples)
    T_out = flt_output_smpltype(fi)
    n_out = flt_output_length(fi)
    output = similar(input, T_out, n_out)
    rdfilt!(output, fi, input)
end


function rdfilt(fi::AbstractRadSigFilterInstance, input::RDWaveform)
    y = rdfilt(fi, input.signal)
    t = flt_output_time_axis(fi, input.time)
    RDWaveform(t, y)
end


#=
function rdfilt!(output::RDWaveform, fi::AbstractRadSigFilterInstance, input::RDWaveform)
    @argcheck output.time == flt_output_time_axis(fi, input.time)
    y = rdfilt!(output.signal, fi, input.signal)
    output
end
=#


"""
    RadiationDetectorDSP::adapt_memlayout(
        fi::AbstractRadSigFilterInstance,
        backend::KernelAbstractions.Backend,
        A::AbstractArray{<:Number}
    )

Adapts the memory layout of `A` in a suitable fashion for `fi` on computing
device `backend`.

Returns a column-major version of `A` on all backends by default, filter instance
types may specialize this behavior.
"""
function adapt_memlayout end

adapt_memlayout(::AbstractRadSigFilterInstance, ::_KA_Backend, A::AbstractArray{<:Number}) = _col_major(A)

# Row-major on GPU backends, column-major on the CPU. KernelAbstractions v0.10 turned
# `GPU` into an alias of `Backend`, so GPU backends can't be dispatched on directly.
_row_major_on_gpu(::KernelAbstractions.CPU, A::AbstractArray{<:Number}) = _col_major(A)
_row_major_on_gpu(::_KA_Backend, A::AbstractArray{<:Number}) = _row_major(A)


"""
    bc_rdfilt(flt::AbstractRadSigFilter, input)
    bc_rdfilt(fi::AbstractRadSigFilterInstance, input)

Broadcast filter instance `fi` over signals `input`, return the filtered
signals.
"""
function bc_rdfilt end

# bc_rdfilt(flt::AbstractRadSigFilter, inputs) = bc_rdfilt(fltinstance(flt, elsmplinfo(inputs)), inputs)

function bc_rdfilt(fi::AbstractRadSigFilterInstance, inputs)
    broadcast((input -> rdfilt(fi, input)), inputs)
end

function bc_rdfilt(
    fi::AbstractRadSigFilterInstance,
    inputs::ArrayOfSimilarArrays{T_in,M,N}
) where {T_in<:Number,M,N}
    T_out = flt_output_smpltype(fi)
    n_out = flt_output_length(fi)
    flat_inputs = flatview(inputs)
    backend = _ka_get_backend(flat_inputs)
    new_flat_inputs = adapt_memlayout(fi, backend, flat_inputs)
    flat_output = _similar_memlayout(new_flat_inputs, T_out, (n_out, size(inputs)...))
    new_inputs = ArrayOfSimilarArrays{T_in,M,N}(new_flat_inputs)
    outputs = ArrayOfSimilarArrays{eltype(flat_output),M,N}(flat_output)
    bc_rdfilt!(outputs, fi, new_inputs)
end

function bc_rdfilt(
    fi::AbstractRadSigFilterInstance,
    inputs::ArrayOfRDWaveforms{<:RealQuantity}
)
    output_signal = bc_rdfilt(fi, inputs.signal)
    output_time = broadcast(t -> flt_output_time_axis(fi, t), inputs.time)
    return ArrayOfRDWaveforms((output_time, output_signal))
end


"""
    bc_rdfilt!(outputs, fi::AbstractRadSigFilterInstance, inputs)

Broadcast filter `flt` or filter instance `fi` over signals `inputs`, storing
the results in `outputs`.

`inputs` and `outputs` must be of type `AbstractVector{<:AbstractSamples}`.

Returns `outputs`.
"""
function bc_rdfilt! end

# bc_rdfilt!(outputs, flt::AbstractRadSigFilter, inputs) = bc_rdfilt!(outputs, fltinstance(flt, elsmplinfo(inputs)), inputs)

function bc_rdfilt!(
    outputs::AbstractVector{<:AbstractSamples},
    fi::AbstractRadSigFilterInstance,
    inputs::AbstractVector{<:AbstractSamples}
)
    Base.Threads.@threads for i in eachindex(outputs)
        rdfilt!(outputs[i], fi, inputs[i])
    end
    return outputs
end


@kernel function _ka_filter_kernel!(
    outputs::ArrayOfSimilarArrays{<:RealQuantity,M,N}, @Const(inputs::ArrayOfSimilarArrays{<:RealQuantity,M,N}),
    @Const(fi::AbstractRadSigFilterInstance)
) where {N,M}
    idxs = @index(Global, NTuple)
    rdfilt!(outputs[idxs...], fi, inputs[idxs...])
end

# Work-group size for kernels that write to `A`:
_ka_threads(::KernelAbstractions.CPU, A::AbstractArray) = (Base.Threads.nthreads(),)
_ka_threads(::_KA_Backend, A::AbstractArray) = ()

function _ka_bc_rdfilt!(
    outputs::ArrayOfSimilarVectors{<:Number},
    fi::AbstractRadSigFilterInstance,
    inputs::ArrayOfSimilarVectors{<:Number}
)
    @argcheck axes(inputs) == axes(outputs)

    adaptor = get_compute_unit(inputs)
    adapted_fi = adapt(adaptor, fi)

    backend = _ka_get_backend(flatview(outputs))
    if _use_staged_kernel(backend, fi, flatview(inputs))
        return _ka_staged_bc_rdfilt!(backend, outputs, adapted_fi, inputs)
    end
    kernel! = _ka_filter_kernel!(backend, _ka_threads(backend, flatview(outputs))...)
    kernel_ret = kernel!(outputs, inputs, adapted_fi, ndrange=size(inputs))
    _ka_synchronize(kernel!, kernel_ret)
    return outputs
end


# Whether `rdfilt!(x, fi, x)` filters `x` in place, i.e. the output has the length of the input
# and sample `i` of the input is not read after output sample `i` is written:
_rdfilt_in_place(::AbstractRadSigFilterInstance) = false

# On GPUs, filters that run in place are applied by a kernel that holds whole signals in shared
# memory (_ka_staged_filter_kernel!): loads and stores are contiguous in either memory layout,
# while each work item filters one signal as on the CPU.
_use_staged_kernel(::KernelAbstractions.CPU, fi::AbstractRadSigFilterInstance, X::AbstractArray) = false
_use_staged_kernel(::_KA_Backend, fi::AbstractRadSigFilterInstance, X::AbstractArray) =
    _rdfilt_in_place(fi) && X isa AbstractMatrix && size(X, 1) <= _ka_smem_smpls

# Memory layout for filters that use the staged kernel on GPUs, row-major otherwise:
_staged_memlayout(fi::AbstractRadSigFilterInstance, backend::_KA_Backend, A::AbstractArray{<:Number}) =
    _use_staged_kernel(backend, fi, A) ? A : _row_major_on_gpu(backend, A)

function _ka_staged_bc_rdfilt!(backend::_KA_Backend, outputs::ArrayOfSimilarVectors, fi::AbstractRadSigFilterInstance, inputs::ArrayOfSimilarVectors)
    X, Y = flatview(inputs), flatview(outputs)
    @argcheck size(X) == size(Y)
    isempty(X) && return outputs
    col_major_in, col_major_out = !(X isa Transpose), !(Y isa Transpose)
    # Filtering is one signal per work item, so as many signals per work group as fit:
    n_grp = _ka_smem_group(size(X, 1), col_major_in && col_major_out)
    kernel! = _ka_staged_filter_kernel!(backend, 256)
    kernel_ret = kernel!(Y, X, fi, n_grp, col_major_in, col_major_out, ndrange = 256 * cld(size(X, 2), n_grp))
    _ka_synchronize(kernel!, kernel_ret)
    return outputs
end

@kernel function _ka_staged_filter_kernel!(
    Y::AbstractMatrix{<:Real}, @Const(X::AbstractMatrix{<:Real}), @Const(fi::AbstractRadSigFilterInstance),
    n_grp::Int, col_major_in::Bool, col_major_out::Bool
)
    buf = @localmem eltype(Y) (8192,)
    n = size(X, 1)
    s0 = (@index(Group, Linear) - 1) * n_grp
    n_sig = min(n_grp, size(X, 2) - s0)

    _ka_load_signals!(buf, X, @index(Local, Linear) - 1, s0, n_sig, n_grp, col_major_in)
    @synchronize

    s = @index(Local, Linear) - 1
    if s < n_sig
        signal = view(buf, (s * n + 1):(s * n + n))
        rdfilt!(signal, fi, signal)
    end
    @synchronize

    _ka_store_signals!(Y, buf, @index(Local, Linear) - 1, s0, n_sig, n_grp, col_major_out)
end


# GPU kernels that hold whole signals in shared memory use work groups of 256 work items and
# n_grp signals of n <= 8192 / n_grp samples in a buffer of 8192 elements, sample i of signal s
# of the group (both 0-based) at element s * n + i + 1. With row-major data, n_grp must be a
# power of two.
const _ka_smem_smpls = 8192

_ka_smem_group(n::Integer, col_major::Bool) = col_major ? min(256, _ka_smem_smpls ÷ n) : min(256, prevpow(2, _ka_smem_smpls ÷ n))

# Loads signals s0+1:s0+n_sig of X into buf, by work item tt (0-based). All loads of a work item
# are issued before its first store, to keep enough memory accesses in flight.
@inline function _ka_load_signals!(buf, X, tt, s0, n_sig, n_grp, col_major)
    T = eltype(buf)
    n = size(X, 1)
    if col_major
        # The signals are consecutive in linear indexing order:
        n_total = n_sig * n
        offs = s0 * n
        Base.Cartesian.@nexprs 32 k -> (v_k = tt + 256 * (k - 1) < n_total ? T(@inbounds X[offs + tt + 256 * (k - 1) + 1]) : zero(T))
        Base.Cartesian.@nexprs 32 k -> (tt + 256 * (k - 1) < n_total && (@inbounds buf[tt + 256 * (k - 1) + 1] = v_k; true))
    else
        # n_grp consecutive work items read the signals at the same sample:
        s = tt & (n_grp - 1)
        i_1 = tt ÷ n_grp
        i_step = 256 ÷ n_grp
        Base.Cartesian.@nexprs 32 k -> (v_k = (s < n_sig && i_1 + i_step * (k - 1) < n) ? T(@inbounds X[i_1 + i_step * (k - 1) + 1, s0 + s + 1]) : zero(T))
        Base.Cartesian.@nexprs 32 k -> ((s < n_sig && i_1 + i_step * (k - 1) < n) && (@inbounds buf[s * n + i_1 + i_step * (k - 1) + 1] = v_k; true))
    end
    return nothing
end

# Stores buf to signals s0+1:s0+n_sig of Y, the inverse of _ka_load_signals!.
@inline function _ka_store_signals!(Y, buf, tt, s0, n_sig, n_grp, col_major)
    n = size(Y, 1)
    if col_major
        for s in 0:(n_sig - 1), i in tt:256:(n - 1)
            @inbounds Y[i + 1, s0 + s + 1] = buf[s * n + i + 1]
        end
    else
        s = tt & (n_grp - 1)
        if s < n_sig
            for i in (tt ÷ n_grp):(256 ÷ n_grp):(n - 1)
                @inbounds Y[i + 1, s0 + s + 1] = buf[s * n + i + 1]
            end
        end
    end
    return nothing
end


"""
    RadiationDetectorDSP.flt_output_smpltype(fi::AbstractRadSigFilterInstance)

Get the output sample type for
    
* a filter `flt` given an input sample type `input_smpltype`
* a filter instance `fi`

Must be implemented for all subtypes of [`AbstractRadSigFilter`](@ref).
"""
function flt_output_smpltype end


"""
    RadiationDetectorDSP.flt_input_smpltype(fi::AbstractRadSigFilterInstance)

Get the input sample type of a filter instance `fi`.
    
Must be implemented for all subtypes of [`AbstractRadSigFilterInstance`](@ref).
"""
function flt_input_smpltype end


"""
    RadiationDetectorDSP.flt_output_length(fi::SomeFilterInstance)::Integer

Get the output signal length of filter instance `fi`.

Must be implemented for all subtypes of [`AbstractRadSigFilterInstance`](@ref).
"""
function flt_output_length end


"""
    RadiationDetectorDSP.flt_input_length(fi::AbstractRadSigFilterInstance)::Integer

Get the output signal length of a filter instance `fi`.

Must be implemented for all subtypes of [`AbstractRadSigFilterInstance`](@ref).
"""
function flt_input_length end


"""
    RadiationDetectorDSP.flt_output_time_axis(fi::SomeFilterInstance, time::AbstractVector{<:RealQuantity})::AbstractVector{<:RealQuantity}

Get the output time axis of a filter instance `fi`, given an input time axis
`time`.

Must be implemented for subtypes of [`AbstractRadSigFilter`](@ref) and
[`AbstractRadSigFilterInstance`](@ref) only if the filter's output time axis
can be computed directly from the input time axis.
"""
function flt_output_time_axis end



#function check_input_compat(fi::AbstractRadSigFilterInstance, input::RDWaveform)
#    @argcheck flt_input_smpltype(fi) == eltype(input)
#    @argcheck flt_input_length(fi) == length(eachindex(input))
#end


#function create_output(ArrayType::Type{<:AbstractArray}, fi::AbstractRadSigFilterInstance)
#    T_out = flt_output_smpltype(fi)
#    n_out = flt_output_length(fi)
#    ArrayType(Fill(zero(T_out), n_out))
#end



"""
    abstract type AbstractRadIIRFilter <: AbstractRadSigFilter{LinearFiltering}

Abstract type for IIR filters.
"""
abstract type AbstractRadIIRFilter <: AbstractRadSigFilter{LinearFiltering} end
export AbstractRadIIRFilter



"""
    abstract type AbstractRadFIRFilter <: AbstractRadSigFilter{LinearFiltering}

Abstract type for FIR filters.
"""
abstract type AbstractRadFIRFilter <: AbstractRadSigFilter{LinearFiltering} end
export AbstractRadFIRFilter



# ToDo:
#=
"""
    RadiationDetectorDSP.getoperator(fi::AbstractRadSigFilterInstance{LinearFiltering})::Matrix

Return a matrix/operator representation of filter instance `fi`.
"""
function getoperator end
export getoperator

RadiationDetectorDSP.getoperator(fi::AbstractRadSigFilterInstance{LinearFiltering}) = ...
=#
