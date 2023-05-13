# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT)

Base.@propagate_inbounds _get_or_view(x, idxs...) = getindex(x, idxs...)

Base.@propagate_inbounds _get_or_view(x::AbstractArray, idxs::Vararg{N,T}) where {N,T<:Integer} = getindex(x, idxs...)
Base.@propagate_inbounds _get_or_view(x::AbstractArray, idxs...) = _get_or_view_array(x, Base.to_indices(x, idxs))
Base.@propagate_inbounds _get_or_view_array(x::AbstractArray, idxs::NTuple{N,T}) where {N,T<:Integer} = getindex(x, idxs...)
Base.@propagate_inbounds _get_or_view_array(x::AbstractArray, idxs::Tuple) = view(x, idxs...)

Base.@propagate_inbounds _get_or_view(x::AbstractRange, idxs::Vararg{N,T}) where {N,T<:Integer} = getindex(x, idxs...)
Base.@propagate_inbounds _get_or_view(x::AbstractRange, idxs...) = getindex(x, idxs...)


Base.@propagate_inbounds _firstof(x) = first(x)
Base.@propagate_inbounds _firstof(x::Union{Tuple,NamedTuple}) = x[1]
Base.@propagate_inbounds _firstof(x::Ref) = x[]


Base.@propagate_inbounds _secondof(x) = x[firstindex(x) + 1]
Base.@propagate_inbounds _secondof(x::Union{Tuple,NamedTuple}) = x[2]
Base.@propagate_inbounds _secondof(x::Ref) = throw(BoundsError(x, 2))


Base.@propagate_inbounds _thirdof(x) = x[firstindex(x) + 2]
Base.@propagate_inbounds _thirdof(x::Union{Tuple,NamedTuple}) = x[3]
Base.@propagate_inbounds _thirdof(x::Ref) = throw(BoundsError(x, 3))


Base.@propagate_inbounds _lastof(x) = last(x)
Base.@propagate_inbounds _lastof(x::NTuple{N,<:Any}) where N = x[N]
Base.@propagate_inbounds _lastof(x::NamedTuple) = _lastof(values(x))
Base.@propagate_inbounds _lastof(x::Ref) = x[]


Base.@propagate_inbounds _all_except_firstof(x) = _get_or_view(x, firstindex(x)+1:lastindex(x))
Base.@propagate_inbounds _all_except_firstof(::Tuple{}) = () # Base.tail doesn't allow empty tuples
Base.@propagate_inbounds _all_except_firstof(x::Tuple) = Base.tail(x)
Base.@propagate_inbounds _all_except_firstof(x::NamedTuple{names}) where names = NamedTuple{_all_except_firstof(names)}(_all_except_firstof(values(x)))
Base.@propagate_inbounds _all_except_firstof(::Ref) = ()


Base.@propagate_inbounds _all_except_lastof(x) = _get_or_view(x, firstindex(x):lastindex(x)-1)
Base.@propagate_inbounds _all_except_lastof(x::NTuple{N,<:Any}) where N = _firstn_of(x, Val(N-1))
Base.@propagate_inbounds _all_except_lastof(x::NamedTuple{names}) where names = NamedTuple{_all_except_lastof(names)}(_all_except_lastof(values(x)))
Base.@propagate_inbounds _all_except_lastof(::Ref) = ()


Base.@propagate_inbounds _firstn_of(x, n::Integer) = _get_or_view(x, firstindex(x):firstindex(x)+n-1)
Base.@propagate_inbounds _firstn_of(x::Union{Tuple,NamedTuple,Ref}, n::Integer) = Base.ntuple(i -> x[i], n)
Base.@propagate_inbounds _firstn_of(x::NamedTuple{names}, n::Integer) where names = NamedTuple{_firstn_of(names, n)}(_firstn_of(values(x), n))
Base.@propagate_inbounds _firstn_of(x::Ref, n::Integer) = _firstn_of((x[],), n)

Base.@propagate_inbounds _firstn_of(x, ::Val{N}) where N = Base.ntuple(i -> x[i], Val{N}())
Base.@propagate_inbounds _firstn_of(x::NamedTuple{names}, ::Val{N}) where {names,N} = NamedTuple{_firstn_of(names, Val{N}())}(_firstn_of(values(x), Val{N}()))
Base.@propagate_inbounds _firstn_of(x::Ref, ::Val{N}) where N = _firstn_of((x[],), Val{N}())


Base.@propagate_inbounds _lastn_of(x, n::Integer) = _get_or_view(x, lastindex(x)+1-n:lastindex(x))
Base.@propagate_inbounds _lastn_of(x::NTuple{M,<:Any}, n::Integer) where M = Base.ntuple(i -> x[i + (M - n)], n)
Base.@propagate_inbounds _lastn_of(x::NamedTuple{names}, n::Integer) where names = NamedTuple{_lastn_of(names, n)}(_lastn_of(values(x), n))
Base.@propagate_inbounds _lastn_of(x::Ref, n::Integer) = _lastn_of((x[],), n)

Base.@propagate_inbounds _lastn_of(x, ::Val{N}) where N = Base.ntuple(i -> x[lastindex(x) - N + i], Val{N}())
Base.@propagate_inbounds _lastn_of(x::NTuple{M,<:Any}, ::Val{N}) where {M,N} = Base.ntuple(i -> x[i + (M - N)], Val{N}())
Base.@propagate_inbounds _lastn_of(x::NamedTuple{names}, ::Val{N}) where {names,N} = NamedTuple{_lastn_of(names, Val{N}())}(_lastn_of(values(x), Val{N}()))
Base.@propagate_inbounds _lastn_of(x::Ref, ::Val{N}) where N = _lastn_of((x[],), Val{N}())


Base.@propagate_inbounds _split_firstof(x) = _firstof(x), _lastn_of(x, length(eachindex(x)) - 1)
Base.@propagate_inbounds _split_firstof(x::NTuple{M,<:Any}) where M = _firstof(x), _lastn_of(x, Val{M - 1}())
Base.@propagate_inbounds _split_firstof(x::NamedTuple{names,<:NTuple{M,<:Any}}) where {names,M} = _firstof(x), _lastn_of(x, Val{M - 1}())
Base.@propagate_inbounds _split_firstof(x::Ref) = _split_firstof((x[],))

Base.@propagate_inbounds _split_lastof(x) = _firstn_of(x, length(eachindex(x)) - 1), _lastof(x)
Base.@propagate_inbounds _split_lastof(x::NTuple{M,<:Any}) where M = _firstn_of(x, Val{M - 1}()), _lastof(x)
Base.@propagate_inbounds _split_lastof(x::NamedTuple{names,<:NTuple{M,<:Any}}) where {names,M} = _firstn_of(x, Val{M - 1}()), _lastof(x)
Base.@propagate_inbounds _split_lastof(x::Ref) = _split_lastof((x[],))


Base.@propagate_inbounds _split_firstn_of(x, n::Integer) = _firstn_of(x, n), _lastn_of(x, length(eachindex(x)) - n)
Base.@propagate_inbounds _split_firstn_of(x::NTuple{M,<:Any}, n::Integer) where M = firstn_of(x, n), _lastn_of(x, M - n)
Base.@propagate_inbounds _split_firstn_of(x::NamedTuple{names,<:NTuple{M,<:Any}}, n::Integer) where {names,M} = _firstn_of(x, n), _lastn_of(x, M - n)
Base.@propagate_inbounds _split_firstn_of(x::Ref, n::Integer) = _split_firstn_of((x[],), n)

Base.@propagate_inbounds _split_firstn_of(x, ::Val{N}) where N = _firstn_of(x, Val{N}()), _lastn_of(x, length(eachindex(x)) - N)
Base.@propagate_inbounds _split_firstn_of(x::NTuple{M,<:Any}, ::Val{N}) where {M,N} = _firstn_of(x, Val{N}()), _lastn_of(x, Val{M - N}())
Base.@propagate_inbounds _split_firstn_of(x::NamedTuple{names,<:NTuple{M,<:Any}}, ::Val{N}) where {names,M,N} = _firstn_of(x, Val{N}()), _lastn_of(x, Val{M - N}())
Base.@propagate_inbounds _split_firstn_of(x::Ref, ::Val{N}) where N = _split_firstn_of((x[],), Val{N}())


Base.@propagate_inbounds _split_lastn_of(x, n::Integer) = _firstn_of(x, length(eachindex(x)) - n), _lastn_of(x, n)
Base.@propagate_inbounds _split_lastn_of(x::NTuple{M,<:Any}, n::Integer) where M = _firstn_of(x, M - n), _lastn_of(x, n)
Base.@propagate_inbounds _split_lastn_of(x::NamedTuple{names,<:NTuple{M,<:Any}}, n::Integer) where {names,M} = _firstn_of(x, M - n), _lastn_of(x, n)
Base.@propagate_inbounds _split_lastn_of(x::Ref, n::Integer) = _split_lastn_of((x[],), n)

Base.@propagate_inbounds _split_lastn_of(x, ::Val{N}) where N = _firstn_of(x, length(eachindex(x)) - N), _lastn_of(x, Val{N}())
Base.@propagate_inbounds _split_lastn_of(x::NTuple{M,<:Any}, ::Val{N}) where {M,N} = _firstn_of(x, Val{M - N}()), _lastn_of(x, Val{N}())
Base.@propagate_inbounds _split_lastn_of(x::NamedTuple{names,<:NTuple{M,<:Any}}, ::Val{N}) where {names,M,N} = _firstn_of(x, Val{M - N}()), _lastn_of(x, Val{N}())
Base.@propagate_inbounds _split_lastn_of(x::Ref, ::Val{N}) where N = _split_lastn_of((x[],), Val{N}())


Base.@pure _ncolons(::Val{N}) where N = ntuple(_ -> Colon(), Val{N}())
