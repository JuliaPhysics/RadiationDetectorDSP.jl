# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


_getdummyelem(A::AbstractArray) = first(A)
_getdummyelem(A::FillArrays.Fill) = A.value


_bcgetindex(A::AbstractArray{<:AbstractArray}, idxs...) = broadcast(getindex, A, map(Ref, idxs)...)

function _bcgetindex(A::ArrayOfSimilarArrays{T,M,N}, idxs...) where {T,M,N}
    flat_A = flatview(A)
    flat_B = getindex(flat_A, ntuple(_ -> :, Val(M))..., idxs...)
    #ArrayOfSimilarArrays{T,M,N + length(size(flat_B)) - length(size(flat_A))}(flat_B)
    _bcgetindex_renest(flat_B, Val(M), Val{N + length(size(flat_B)) - length(size(flat_A))}())
end

_bcgetindex_renest(flat_B::AbstractArray{T}, ::Val{M}, ::Val{N}) where {T,M,N} = ArrayOfSimilarArrays{T,M,N}(flat_B)
_bcgetindex_renest(flat_B::AbstractArray{T}, ::Val{M}, ::Val{0}) where {T,M} = flat_B
