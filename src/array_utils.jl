# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


#_getdummyelem(A::AbstractArray) = first(A)
#_getdummyelem(A::FillArrays.Fill) = A.value

function _uniqueelem(A::AbstractArray)
    x = first(A)
    @argcheck all(isequal(x), A)
    return x
end

_uniqueelem(A::FillArrays.Fill) = A.value


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


_similar_memlayout(A::AbstractArray{<:T}, ::Type{U}, sz::Dims) where {T<:Real,U<:Real} = similar(A, U, sz)
_similar_memlayout(A::LinearAlgebra.Transpose{<:T}, ::Type{U}, sz::Dims) where {T<:Real,U<:Real} = transpose(similar(A, U, reverse(sz)))

_row_major(A::AbstractMatrix{T}) where {T<:Number} = _lazy_transpose(_nonlazy_transpose(A))
_row_major(A::LinearAlgebra.Transpose{T}) where {T<:Number} = A

_col_major(A::AbstractMatrix{T}) where {T<:Number} = A
_col_major(A::LinearAlgebra.Transpose{T}) where {T<:Number} = _nonlazy_transpose(_lazy_transpose(A))


@inline _kbc_view(A::AbstractArray{<:Any,N}, idxs...) where N =  view(A, _firstn_of(idxs, Val(N))...)

@inline _kbc_getindex(A::AbstractArray{<:Any,N}, idxs...) where N =  getindex(A, _firstn_of(idxs, Val(N))...)
@inline _kbc_getindex(x::Number, idxs...) =  x
@inline _kbc_getindex(x::Tuple{<:Number}, idxs...) =  x[1]
@inline _kbc_getindex(x::Ref{<:Number}, idxs...) =  x[]

@inline _kbc_setindex!(A::AbstractArray{<:Any,N}, x, idxs...) where N =  setindex!(A, x, _firstn_of(idxs, Val(N))...)



const CollectionLike{T} = Union{AbstractArray{T},(NTuple{N,T} where N),Ref{T}}
const FlattableCollection{T} = Union{CollectionLike{T},CollectionLike{<:CollectionLike{T}}}

@inline _bc_flat_getindex(A::AbstractArray{<:Number,N}, idxs...) where N =  getindex(A, _front_tuple(idxs, Val(N))...)
@inline _bc_flat_getindex(x::Number, idxs::Integer...) =  x
@inline _bc_flat_getindex(x::Tuple{<:Number}, ::Integer, idxs::Integer...) =  x[1]
@inline _bc_flat_getindex(x::Ref{<:Number}, idxs::Integer...) =  x[]
@inline _bc_flat_getindex(x::NTuple{N,<:Number}, idxs::Integer...) where N =  x[_firstof(idxs)]

@inline function _bc_flat_getindex(x::AbstractArray{<:CollectionLike{<:Number},N}, idxs...) where N
    idxs_inner, idxs_outer = _split_lastn_of(idxs, Val{N}())
    _bc_flat_getindex(_bc_flat_getindex(x, idxs_outer...), idxs_inner...)
end

@inline function _bc_flat_getindex(x::Union{(NTuple{N,<:CollectionLike{<:Number}} where N),Ref{<:CollectionLike{<:Number}}}, idxs...)
    idxs_inner, i_outer = _split_lastof(idxs)
    _bc_flat_getindex(_bc_flat_getindex(x, i_outer), idxs_inner...)
end



struct _Reslice{M,N} end
(f::_Reslice{M,N})(A::AbstractArray{T}) where {M,N,T<:RealQuantity} = ArrayOfSimilarArrays{T,M,N}(A)

_kbc_flatview(x::RealQuantity) = x, identity
_kbc_flatview(ref::Ref{<:RealQuantity}) = ref[], identity
_kbc_flatview(ref::Tuple{<:RealQuantity}) = ref[], identity
_kbc_flatview(ref::Ref{<:AbstractArray{<:RealQuantity}}) = ref[], identity
_kbc_flatview(ref::Tuple{<:AbstractArray{<:RealQuantity}}) = ref[], identity
_kbc_flatview(As::ArrayOfSimilarArrays{<:RealQuantity,M,N}) where {M,N} = flatview(As), _Reslice{M,N}()

_kbc_size(::Tuple{}) = (1,)
_kbc_size(sz::Tuple) = sz

_kbc_result(x) = x
_kbc_result(x::AbstractArray{<:Any,0}) = x[]


const _BC_RQs = Union{AbstractArray{<:RealQuantity}, Ref{<:RealQuantity}, Tuple{<:RealQuantity}, RealQuantity}
const _BC_RQ_Arrays = Union{AbstractArray{<:AbstractArray{<:RealQuantity}}, Ref{<:AbstractArray{<:RealQuantity}}, Tuple{<:AbstractArray{<:RealQuantity}}}
const _BC_RQ_AosAs = Union{ArrayOfSimilarArrays{<:RealQuantity}, Ref{<:AbstractArray{<:RealQuantity}}, Tuple{<:AbstractArray{<:RealQuantity}}}


_similar_maybe_structarray(A, ::Type{T}, sz::Dims) where {T<:Real} =  similar(A, T, sz)
_similar_maybe_structarray(A, ::Type{T}, sz::Dims) where T =  similar(StructArray((A,)), T, sz)


const GPULikeArray{T,N} = Union{AbstractGPUArray{T,N},SubArray{T,N,<:AbstractGPUArray{T,N}}}


"""
    RadiationDetectorDSP.CPUNormAdaptor

To be used with `Adapt.adapt`.

`Adapt.adapt(RadiationDetectorDSP.CPUNormAdaptor, x)` adapts `x` to reside on
the CPU and tries to ensure that arrays are stored in column-major order.
"""
struct CPUNormAdaptor end

Adapt.adapt_storage(::CPUNormAdaptor, A::AbstractArray) = adapt(Array, A)
Adapt.adapt_structure(to::CPUNormAdaptor, A::LinearAlgebra.Transpose{<:Number}) = adapt(to, _nonlazy_transpose(transpose(A)))
