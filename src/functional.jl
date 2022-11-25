# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    abstract type AbstractRadSigFunctional <: Function

Abstract type for signal functionals (feature extractors).

Functionals are callable as `(func::AbstractRadSigFunctional)(input)` and come
with specialized broadcasting.

Subtypes of `AbstractRadSigFunctional` must implement

```
funcinstance(func::AbstractRadSigFunctional, si::SamplingInfo)::[`AbstractRadSigFunctionalInstance`](@ref)
```
"""
abstract type AbstractRadSigFunctional <: Function end
export AbstractRadSigFunctional

(func::AbstractRadSigFunctional)(input) = rdfunc(funcinstance(func, smplinfo(input)), input)

# ToDo: Support mutating broadcasts (needs a bc_rdfunc! that takes chained filters into account)
function Base.Broadcast.broadcasted(func::AbstractRadSigFunctional, inputs)
    X = Base.materialize(inputs)
    fi = funcinstance(func, elsmplinfo(inputs))
    bc_rdfunc(fi, inputs)
end



"""
    abstract type AbstractRadSigFunctionalInstance

Abstract type for signal functionals. Functional instances are specilized to
a specific length and numerical type of signal input and value output.

Filter instances are callable as `(fi::SomeFunctionalInstance)(input)` and come
with specialized broadcasting.

Subtypes of `AbstractRadSigFunctionalInstance` must implement

* `RadiationDetectorDSP.rdfunc(fi::SomeFunctionalInstance, input)`

Default methods are implemented for

* `RadiationDetectorDSP.rdfunc(fi::AbstractRadSigFunctionalInstance, wf::RDWaveform)`
* `RadiationDetectorDSP.bc_rdfunc(fi::AbstractRadSigFunctionalInstance, inputs)`
"""
abstract type AbstractRadSigFunctionalInstance end
export AbstractRadSigFunctionalInstance


"""
    funcinstance(func::AbstractRadSigFunctional, si::SamplingInfo)::AbstractRadSigFunctionalInstance

Create a functional instance of the functional `func`, specialized for the given
input, resp. input characteristics.
"""
function funcinstance end
export funcinstance

funcinstance(func::AbstractRadSigFunctional, si::SamplingInfo{T}) where T = throw(ArgumentError("funcinstance not defined for type $(nameof(typeof(func)))"))


"""
    rdfunc::AbstractRadSigFunctionalInstance, input)

Apply functional instance `fi` to signal `input` and return the result.
"""
function rdfunc end
export rdfunc

rdfunc(fi::AbstractRadSigFunctionalInstance, input::RDWaveform) = rdfunc(fi, input.signal)


"""
    bc_rdfunc(outputs, fi::AbstractRadSigFunctionalInstance, inputs)

Broadcast filter `func` or filter instance `fi` over signals `inputs`, store
the results in `outputs`.

`inputs` and `outputs` must be of type `AbstractVector{<:AbstractSamples}`.

Returns `outputs`.
"""
function bc_rdfunc end

function bc_rdfunc(
    fi::AbstractRadSigFunctionalInstance,
    inputs::ArrayOfRDWaveforms{<:RealQuantity}
)
    bc_rdfunc(fi, inputs.signal)
end

function bc_rdfunc(
    fi::AbstractRadSigFunctionalInstance,
    inputs::AbstractVector{<:AbstractSamples}
)
    broadcast((input -> rdfunc(fi, input)), inputs)
end
