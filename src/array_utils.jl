# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


#_getdummyelem(A::AbstractArray) = first(A)
#_getdummyelem(A::FillArrays.Fill) = A.value

function _uniqueelem(A::AbstractArray)
    x = first(A)
    @argcheck all(isequal(x), A)
    return x
end

_uniqueelem(A::FillArrays.Fill) = A.value


_bcgetindex(A::AbstractArray{<:AbstractArray}, idxs...) = broadcast(getindex, A, map(Ref, idxs)...)

function _bcgetindex(A::ArrayOfSimilarArrays{T,M,N}, idxs...) where {T,M,N}
    flat_A = flatview(A)
    flat_B = getindex(flat_A, ntuple(_ -> :, Val(M))..., idxs...)
    #ArrayOfSimilarArrays{T,M,N + length(size(flat_B)) - length(size(flat_A))}(flat_B)
    _bcgetindex_renest(flat_B, Val(M), Val{N + length(size(flat_B)) - length(size(flat_A))}())
end

_bcgetindex_renest(flat_B::AbstractArray{T}, ::Val{M}, ::Val{N}) where {T,M,N} = ArrayOfSimilarArrays{T,M,N}(flat_B)
_bcgetindex_renest(flat_B::AbstractArray{T}, ::Val{M}, ::Val{0}) where {T,M} = flat_B


function _inneraxes(A::AbstractArray{<:AbstractArray{T,M},N}) where {T,M,N}
    ax = if !isempty(A)
        axes(first(A))
    else
        ntuple(_ ->  Base.OneTo(0), Val(M))
    end

    let ax = ax
        if any(X -> axes(X) != ax, A)
            throw(DimensionMismatch("Axes of element arrays of A are not equal, can't determine common shape"))
        end
    end

    ax
end

@inline _inneraxes(A::AbstractArray{<:AbstractArray}, dim::Integer) =
    _inneraxes(A)[dim]

@inline function _inneraxes(A::ArrayOfSimilarArrays{T,M,N}) where {T,M,N}
    ArraysOfArrays.front_tuple(axes(A.data), Val{M}())
end

@inline function _inneraxes(A::FillArrays.Fill{<:AbstractArray})
    axes(A.value)
end


# ToDo: Replace _to_same_device_as workaround as soon as ArrayInferface
# and Adapt have proper support for computing devices:

_to_same_device_as(::T, Y::T) where T = Y
_to_same_device_as(::Array, Y::Array) = Y
_to_same_device_as(::Array, Y::SubArray{Float32,1,<:Array}) = Y

function _to_same_device_as(X, Y)
    new_Y = similar(X, eltype(Y), size(Y))
    copy!(new_Y, Y)
    return new_Y
end
