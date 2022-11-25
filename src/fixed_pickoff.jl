# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    struct FixedPickoff <: AbstractRadSigFunctional

Retrieves the signal value at a given point in time.

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct FixedPickoff{
    T <: RealQuantity
} <: AbstractRadSigFunctional
    "pickoff time"
    time::T
end

export FixedPickoff

Adapt.adapt_structure(to, func::FixedPickoff) = func

function functinstance(func::FixedPickoff, si::SamplingInfo{T}) where T
    delta_t = step(si.axis)

    i = round(Int, uconvert(NoUnits, func.t / delta_t))

    FixedPickoffInstance(i)
end



struct FixedPickoffInstance <: AbstractRadSigFunctionalInstance
    i::Int
end


@inline function rdfunc(fi::FixedPickoffInstance, x::AbstractVector{<:RealQuantity})
    x[fi.i]
end

function bc_rdfunc(fi::FixedPickoffInstance, inputs::ArrayOfSimilarVectors{<:RealQuantity,N}) where N
    X = flatview(inputs)
    colons = ntuple(_ -> :, Val(N))
    X[fi.i, colons...]
end
