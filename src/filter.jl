# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    abstract type FilteringType

Indended as a type parameter to designate the behavior of a filter
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

Not that while a filter may have an inverse, it may (often depending on the
filter paramters) be very unstable in the presence of additional noise. So
the inverse may be not be useful to deconvolve signals in practical
applications.
"""
abstract type AbstractRadSigFilter{FT<:FilteringType} <: Function end
export AbstractRadSigFilter

(flt::AbstractRadSigFilter)(input) = rdfilt(fltinstance(flt, smplinfo(input)), input)

# ToDo: Support mutating broadcasts (needs a bc_rdfilt! that takes chained filters into account)
function Base.Broadcast.broadcasted(flt::AbstractRadSigFilter, inputs)
    X = Base.materialize(inputs)
    fi = fltinstance(flt, elsmplinfo(inputs))
    bc_rdfilt(fi, inputs)
end



"""
    abstract type AbstractRadSigFilterInstance{FT<:FilteringType}

Abstract type for signal filter instances. Filter instances are specilized to
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
function dfilt! end
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
    inputs::ArrayOfSimilarArrays{<:RealQuantity,M,N}
) where {M,N}
    T_out = flt_output_smpltype(fi)
    n_out = flt_output_length(fi)
    flat_output = similar(flatview(inputs), T_out, n_out, size(inputs)...)
    outputs = ArrayOfSimilarArrays{T_out,M,N}(flat_output)
    bc_rdfilt!(outputs, fi, inputs)
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
    Y::AbstractArray{<:RealQuantity,N}, @Const(X::AbstractArray{<:RealQuantity,N}),
    fi::AbstractRadSigFilterInstance
) where N
    idxs = @index(Global, NTuple)
    rdfilt!(view(Y, :, idxs...), fi, view(X, :, idxs...))
end

function _ka_bc_rdfilt!(
    outputs::ArrayOfSimilarVectors{<:RealQuantity},
    fi::AbstractRadSigFilterInstance,
    inputs::ArrayOfSimilarVectors{<:RealQuantity}
)
    X = flatview(inputs)
    Y = flatview(outputs)
    @argcheck Base.tail(axes(X)) == Base.tail(axes(Y))

    dev = KernelAbstractions.get_device(Y)
    kernel! = _ka_filter_kernel!(dev)
    evt = kernel!(Y, X, fi, ndrange=Base.tail(size(Y))) 
    wait(evt)
    return outputs
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
    abstract type AbstractRadFIRFilter <: AbstractRadLinearFilter

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
