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


"""
    RadiationDetectorDSP.CPUNormAdaptor

To be used with `Adapt.adapt`.

`Adapt.adapt(RadiationDetectorDSP.CPUNormAdaptor, x)` adapts `x` to reside on
the CPU and tries to ensure that arrays are stored in column-major order.
"""
struct CPUNormAdaptor end

Adapt.adapt_storage(::CPUNormAdaptor, A::AbstractArray) = adapt(Array, A)
Adapt.adapt_structure(to::CPUNormAdaptor, A::LinearAlgebra.Transpose{<:Number}) = adapt(to, _nonlazy_transpose(transpose(A)))
