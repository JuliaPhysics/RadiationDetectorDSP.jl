# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

module RadiationDetectorDSPReactantExt

using Reactant: Reactant, TracedRNumber

import RadiationDetectorDSP
using RadiationDetectorDSP: SamplingInfo


const _TracedMatrix = AbstractMatrix{<:TracedRNumber}
const _ConcreteMatrix = Union{Reactant.ConcretePJRTArray{<:Number,2}, Reactant.ConcreteIFRTArray{<:Number,2}}

# The generic method transposes element by element:
RadiationDetectorDSP._nonlazy_transpose(X::_TracedMatrix) = permutedims(X)
RadiationDetectorDSP._nonlazy_transpose(X::_ConcreteMatrix) = permutedims(X)


# Filter instances are parameterized by the numerical type of the samples,
# not by the type of the symbolic numbers that stand in for them:
function RadiationDetectorDSP.elsmplinfo(smpls::AbstractArray{<:AbstractVector{TracedRNumber{T}}}) where {T}
    SamplingInfo{T}(RadiationDetectorDSP._inneraxes(smpls, 1))
end


RadiationDetectorDSP._convert_smpls(::Type{T}, X::AbstractArray{TracedRNumber{T}}) where {T} = X
RadiationDetectorDSP._convert_smpls(::Type{T}, X::AbstractArray{<:TracedRNumber}) where {T} = TracedRNumber{T}.(X)


end # module RadiationDetectorDSPReactantExt
