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
