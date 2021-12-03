# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    const AbstractRadLinearFilter = AbstractRadSigFilter{LinearFiltering}

Convenience type alias, abstract linear filter.
"""
const AbstractRadLinearFilter = AbstractRadSigFilter{LinearFiltering}
export AbstractRadLinearFilter


"""
    const AbstractRadLinearFilterInstance = AbstractRadSigFilterInstance{LinearFiltering}

Convenience type alias, abstract linear filter instance.
"""
const AbstractRadLinearFilterInstance = AbstractRadSigFilterInstance{LinearFiltering}
export AbstractRadLinearFilterInstance



# ToDo: RadLinearFilterChain and RadLinearFilterChainInstance



"""
    struct FilterLinearMap{T} <: LinearMaps.LinearMap{T}

Representation of a linear filter as a `LinearMaps.LinearMap`.

Constructors:

* ```$(FUNCTIONNAME){T<:RealQuantity}(flt::AbstractRadLinearFilterInstance)```
"""
struct FilterLinearMap{
    T <: RealQuantity,
    F <: AbstractRadLinearFilterInstance,
    R <: Range
} <: LinearMaps.LinearMap{T}
    flt::F
end



function FilterLinearMap{T}(flt::AbstractRadLinearFilterInstance) where {T<:RealQuantity}
    FilterLinearMap{T,typeof(flt),typeof(time)}(flt, time)
end

LinearMaps.LinearMap{T}(flt::AbstractRadLinearFilter, time::Range) where {T<:RealQuantity} =
    FilterLinearMap{T}(flt, time)


# ToDo: Allow non-square filter matrices:
Base.size(A::FilterLinearMap) = (length(flt_output_time_axis(A.flt, A.time)), length(A.time))

# LinearAlgebra.issymmetric(A::FilterLinearMap) = false
# LinearAlgebra.ishermitian(A::FilterLinearMap) = false
# LinearAlgebra.isposdef(A::FilterLinearMap) = true

# LinearAlgebra.adjoint(A::FilterLinearMap) = ...

# LinearAlgebra.transpose(A::FilterLinearMap{<:Real}) = ...
# LinearAlgebra.transpose(A::FilterLinearMap) = ...

LinearMaps.MulStyle(A::FilterLinearMap) = LinearMaps.ThreeArg()

@inline function LinearAlgebra.mul!(y::AbstractSamples, A::FilterLinearMap, x::AbstractSamples)
    # ToDo: Use more appropriate exception
    @boundscheck @require (A.time,) == axes(x)
    mul!(y, A.flt, x)
end


# ToDo:
# LinearAlgebra.mul!(y::RDWaveform, A::FilterLinearMap, x::RDWaveform)


# ToDo: Decision, do we want this?
#
# function Base.:(*)(A::FilterLinearMap{T}, B::FilterLinearMap{U}) where {T,U}
#     @require A.time == B.time  # ToDo: Allow non-square filter matrices
#     FilterLinearMap{promote_type(T, U)}(A.flt * B.flt)
# end
