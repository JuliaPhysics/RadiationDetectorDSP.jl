# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    const MaybeWithUnits{T<:Number} = Union{T,Quantity{<:T}}

A numerical value with or without units
"""
const MaybeWithUnits{T<:Number} = Union{T,Quantity{<:T}}


"""
    const RealQuantity = MaybeWithUnits{<:Real}

A real value with or without units.
"""
const RealQuantity = MaybeWithUnits{<:Real}


"""
    const AbstractSamples{T<:RealQuantity} = AbstractVector{T}

A vector of signal samples.
"""
const AbstractSamples{T<:RealQuantity} = AbstractVector{T}


"""
    const ArrayOfSimilarSamples{T<:RealQuantity} = ArrayOfSimilarVectors{T}

An array of similar sample vectors.
"""
const ArrayOfSimilarSamples{T<:RealQuantity} = ArrayOfSimilarVectors{T}


dspfloattype(::Type{T}) where {T} = float(T)
dspfloattype(::Type{T}) where {T<:Union{Int8,UInt8,Int16,UInt16,Int32,UInt32}} = Float32



"""
    struct SamplingInfo{T<:RealQuantity,A<:AbstractVector{<:RealQuantity}}

Holds sampling information.

The numerical type of an individual sample is `T`, the (time) axis is
given by the `axis` field.

Constructors:

* ```$(FUNCTIONNAME){T,A}(axis)```
* ```$(FUNCTIONNAME){T}(axis)```

Fields:

$(TYPEDFIELDS)
"""
struct SamplingInfo{T<:RealQuantity,A<:AbstractVector{<:RealQuantity}}
    axis::A
end

export SamplingInfo

SamplingInfo{T}(axis::A) where {T<:RealQuantity,A<:AbstractVector{<:RealQuantity}} = SamplingInfo{T,A}(axis)

_smpltype(si::SamplingInfo{T}) where T = T
_smpllen(si::SamplingInfo) = size(si.axis, 1)



"""
    smplinfo(smpls::AbstractSamples)::SamplingInfo
    smplinfo(wf::RDWaveform{T,U})::RDWaveform

Get sampling information from a vector of samples, resp. a waveform.
"""
function smplinfo end
export smplinfo

smplinfo(smpls::AbstractSamples) = SamplingInfo{eltype(smpls)}(eachindex(smpls))
smplinfo(wf::RDWaveform{T,U}) where {T,U} = SamplingInfo{U}(wf.time)


"""
    smplinfo(smpls::AbstractSamples)::SamplingInfo
    smplinfo(wf::RDWaveform{T,U})::RDWaveform

Get sampling information an array of vectors of samples, resp. an array of
waveform. All elements must have equal samling information.
"""
function elsmplinfo end
export elsmplinfo

elsmplinfo(smpls::AbstractArray{<:AbstractSamples{T}}) where T = SamplingInfo{T}(_inneraxes(smpls, 1))

elsmplinfo(wfs::AbstractArray{<:RDWaveform{T,U}}) where {T,U} = SamplingInfo{T}(_uniqueelem(map(wf -> wf.time, wfs)))

elsmplinfo(wfs::ArrayOfRDWaveforms{T,U}) where {T,U} = SamplingInfo{T}(_uniqueelem(wfs.time))
