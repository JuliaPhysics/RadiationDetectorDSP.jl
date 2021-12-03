# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    abstract type FilterLinearity <: FilteringType

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
    abstract type LinearFiltering <: FilteringType

When used as a type parameter value, marks non linear behavior of a filter.
"""
abstract type NonlinearFiltering <: FilteringType end
export NonlinearFiltering



"""
    abstract type AbstractRadSigFilter{<:FilteringType}

Abstract type for signal filters.

Subtypes of `AbstractRadSigFilter` must implement

```
fltinstance(
    flt::AbstractRadSigFilter,
    input_smpltype::Type{<:RealQuantity},
    input_length::Integer,
    input_dt::RealQuantity
)::[`AbstractRadSigFilterInstance`](@ref)
```

as well as:

* `RadiationDetectorDSP.flt_output_smpltype(flt::SomeFilter, input_smpltype::Type{<:RealQuantity})`

* `RadiationDetectorDSP.flt_output_length(flt::SomeFilter, input_length::Integer)`

* `RadiationDetectorDSP.flt_output_dt(fi::SomeFilter, input_dt::Real)`

Invertible filters also must implement

* `Base.inv(flt::SomeFilter)`

The default methods for

* `RadiationDetectorDSP.rdfilt(flt::AbstractRadSigFilter, input)`
* `RadiationDetectorDSP.rdfilt!(output, flt::AbstractRadSigFilter, input)`

should not be overloaded for `AbstractRadSigFilter`. Instead, overload these
methods for the filter instance type of `flt`.
"""
abstract type AbstractRadSigFilter{<:FilteringType} end
export AbstractRadSigFilter


"""
    abstract type AbstractRadSigFilterInstance{<:FilteringType}

Abstract type for signal filter instances. Filter instances are specilized to
a specific length and numerical type of input and output.

Subtypes of `AbstractRadSigFilterInstance` must implement

* `RadiationDetectorDSP.rdfilt!(output, fi::SomeFilterInstance, input)`

* `RadiationDetectorDSP.flt_input_smpltype(fi::SomeFilterInstance)`

* `RadiationDetectorDSP.flt_output_smpltype(fi::SomeFilterInstance)`

* `RadiationDetectorDSP.flt_input_length(fi::SomeFilterInstance)`

* `RadiationDetectorDSP.flt_output_length(fi::SomeFilterInstance)`

* `RadiationDetectorDSP.flt_input_dt(fi::SomeFilterInstance)`

* `RadiationDetectorDSP.flt_output_dt(fi::SomeFilterInstance)`

* `RadiationDetectorDSP.filterdef(fi::SomeFilterInstance)::AbstractRadSigFilter`

Linear filters (subtypes of `AbstractRadSigFilterInstance{LinearFiltering}`) must implement

* `Base.convert(LinearMaps.LinearMap, fi::SomeFilterInstance)`

Invertible filters may implement

* `Base.inv(fi::SomeFilterInstance)`

If a filter's output time axis can be computed from the input time axis, it
must implement

* `RadiationDetectorDSP.flt_output_time_axis(fi::SomeFilterInstance, time::Range)`

Default methods are implemented for

* `RadiationDetectorDSP.rdfilt(fi::AbstractRadSigFilterInstance, x::AbstractSamples)`
* `RadiationDetectorDSP.rdfilt(fi::AbstractRadSigFilterInstance, wf::RDWaveform)`
* `RadiationDetectorDSP.rdfilt!(output::RDWaveform, fi::AbstractRadSigFilterInstance, input::RDWaveform)`
* `RadiationDetectorDSP.rdfilt!(output::AbstractArray{<:RDWaveform,0}, fi::AbstractRadSigFilterInstance, input::RDWaveform)`

The default methods that operate on `RadiationDetectorSignals.RDWaveform`s require
[`RadiationDetectorDSP.flt_output_time_axis`](@ref).
"""
abstract type AbstractRadSigFilterInstance{<:FilteringType} end
export AbstractRadSigFilter



"""
    fltinstance(
        flt::AbstractRadSigFilter,
        input_smpltype::Type{<:RealQuantity},
        input_length::Integer,
        input_dt::RealQuantity
    )::AbstractRadSigFilterInstance

    fltinstance(flt::AbstractRadSigFilter, input::AbstractSamples)::AbstractRadSigFilterInstance
    fltinstance(flt::AbstractRadSigFilter, input::RDWaveform)::AbstractRadSigFilterInstance

Create a filter instance of the filter `flt`, specialized for the given
input, resp. input characteristics.
"""
function fltinstance end
export fltinstance

function fltinstance(flt::AbstractRadSigFilter, input::AbstractSamples)
    fltinstance(flt, eltype(input), length(eachindex(input)), one(Int32))
end

function fltinstance(flt::AbstractRadSigFilter, input::RDWaveform)
    fltinstance(flt, eltype(input.value), length(eachindex(input)), step(input.time))
end


"""
    rdfilt!(output, flt::AbstractRadSigFilter, input)
    rdfilt!(output, fi::AbstractRadSigFilterInstance, input)

Apply filter `flt` or filter instance `fi` to signal `input` and store the
filtered signal in `output`. Return `output`.
"""
function dfilt! end
export rdfilt!

rdfilt!(output, flt::AbstractRadSigFilter, input) = rdfilt!(output, fltinstance(flt, input), input)


"""
    rdfilt(flt::AbstractRadSigFilter, input)
    rdfilt(fi::AbstractRadSigFilterInstance, input)

Apply filter `flt` or filter instance `fi` to signal `input`, return the
filtered signal.

Returns `output`.
"""
function rdfilt end
export rdfilt

rdfilt(flt::AbstractRadSigFilter, input) = rdfilt(fltinstance(flt, input), input)


"""
    RadiationDetectorDSP.flt_output_smpltype(flt::AbstractRadSigFilter, input_smpltype::Type{<:RealQuantity})
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
    RadiationDetectorDSP.flt_output_length(flt::AbstractRadSigFilter, input_length::Integer)::Integer

Get the output signal length of

* a filter `flt`, given an input of length `input_length`
* a filter instance `fi`

Must be implemented for all subtypes of [`AbstractRadSigFilter`](@ref) and
[`AbstractRadSigFilterInstance`](@ref).
"""
function flt_output_length end


"""
    RadiationDetectorDSP.flt_input_length(fi::AbstractRadSigFilterInstance)::Integer

Get the output signal length of a filter instance `fi`.

Must be implemented for all subtypes of [`AbstractRadSigFilterInstance`](@ref).
"""
function flt_input_length end


"""
    RadiationDetectorDSP.flt_output_dt(flt::SomeFilter, input_dt::RealQuantity)::RealQuantity
    RadiationDetectorDSP.flt_output_dt(fi::SomeFilterInstance)::RealQuantity

Get the output sampling interval of

* a filter `flt`, given an input sampling interval `input_dt`
* a filter instance `fi`

Must be implemented for all subtypes of [`AbstractRadSigFilter`](@ref) and
[`AbstractRadSigFilterInstance`](@ref).        
"""
function flt_output_dt end


"""
    RadiationDetectorDSP.flt_input_dt(fi::AbstractRadSigFilterInstance)::RealQuantity

Get the input sampling interval of a filter instance `fi`.

Must be implemented for all subtypes of [`AbstractRadSigFilterInstance`](@ref).
"""
function flt_input_dt end


"""
    RadiationDetectorDSP.flt_output_time_axis(fi::AbstractRadSigFilterInstance, time::Range)::Range

Get the output time axis of a filter instance `fi`, given an input time axis
`time`.

Must be implemented for subtypes of [`AbstractRadSigFilter`](@ref) and
[`AbstractRadSigFilterInstance`](@ref) only if the filter's output time axis
can be computed directly from the input time axis.
"""
function flt_output_time_axis end



function check_input_compat(fi::AbstractRadSigFilterInstance, input::RDWaveform)
    @require flt_input_smpltype(fi) == eltype(input)
    @require flt_input_length(fi) == length(eachindex(input))
    @require flt_input_dt(fi) == step(input)
end

function create_output(ArrayType::Type{<:AbstractArray}, fi::AbstractRadSigFilterInstance)
    T_out = flt_output_smpltype(fi)
    n_out = flt_output_length(fi)
    ArrayType(Fill(zero(T_out), n_out))
end


function rdfilt(fi::AbstractRadSigFilterInstance, input::AbstractSamples)
    T_out = flt_output_smpltype(fi)
    n_out = flt_output_length(fi)
    output = similar(x, U, n_out)
    rdfilt!(output, fi, input)
end

# ToDo: Custom broadcasting over ArrayOfSimilarSamples


function rdfilt(fi::AbstractRadSigFilterInstance, input::RDWaveform)
    y = rdfilt(fi, input.value)
    t = flt_output_time_axis(fi, input.time)
    RDWaveform(t, y)
end


function rdfilt!(output::RDWaveform, fi::AbstractRadSigFilterInstance, input::RDWaveform)
    @require output.time == flt_output_time_axis(fi, input.time)
    y = rdfilt!(output.value, fi, input.value)
    output
end


function rdfilt!(output::AbstractArray{<:RDWaveform,0}, fi::AbstractRadSigFilterInstance, input::RDWaveform)
    t = flt_output_time_axis(fi, input.time)
    output_sc = output[]
    rdfilt!(output_sc.value, fi, input.value)
    if output_sc.time != t
        output[] = typeof(output_sc)(t, output_sc.value)
    end
    output
end



# ToDo: Custom broadcasting over ArrayOfRDWaveforms



struct GenericRadSigFilter{
    F<:Function,
    FOT<:Function,
    FOL<:Function,
    FOD<:Function,
    FOA<:Function
} <: AbstractRadSigFilter
    f_apply!::F
    f_output_type::FOT
    f_output_length::FOL
    f_output_dt::FOD
    f_output_time_axis::FOA
end


GenericRadSigFilter(
    f_apply!::Function;
    f_output_type::Function = identity,
    f_output_length::Function = identity,
    f_output_dt::Function = identity,
    f_output_time_axis::Function = identity  
) = GenericRadSigFilter{typeof(f_apply!), typeof(f_output_type), typeof(f_output_length), typeof(f_output_dt), typeof(f_output_time_axis)}(
    f_apply!, f_output_type, f_output_length, f_output_dt, f_output_time_axis
)


flt_output_smpltype(flt::GenericRadSigFilter, input_smpltype::Type{<:RealQuantity}) = flt.f_output_type(input_smpltype)

flt_output_length(flt::GenericRadSigFilter, input_length::Integer) = flt.f_output_length(input_length)

flt_output_time_axis(flt::GenericRadSigFilter, time::Range) = flt.f_output_time_axis(time)


struct GenericRadSigFilterInstance{
    F<:Function,
    TIN<:RealQuantity,
    TOUT<:RealQuantity,
    TIDT<:RealQuantity,
    TODT<:RealQuantity,
    FOA<:Function
} <: AbstractRadSigFilter
    f_apply!::F
    input_length::Integer
    output_length::Integer
    input_dt::TIDT
    output_dt::TODT
    f_output_time_axis::FOA
end

function fltinstance(
    flt::GenericRadSigFilter,
    input_smpltype::Type{<:RealQuantity},
    input_length::Integer,
    input_dt::RealQuantity
)
    output_smpltype = flt_output_smpltype(flt, input_smpltype)
    output_length = flt_output_length(flt, input_length)
    output_dt = flt_output_length(flt, input_dt)

    GenericRadSigFilterInstance{
        typeof(flt.f_apply!), input_smpltype, output_smpltype, typeof(output_length),
        typeof(output_dt), typeof(flt.f_output_time_axis)
    }(flt.f_apply!, input_length, output_length, input_dt, output_dt, flt.f_output_time_axis)
end

rdfilt!(output::AbstractSamples, flt::GenericRadSigFilter, input::AbstractSamples) = flt.f_apply!(output, input)
