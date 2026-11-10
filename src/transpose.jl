# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


_lazy_transpose(A::AbstractMatrix{T}) where {T<:Number} = transpose(A)
_nonlazy_transpose(A::AbstractMatrix{T}) where {T<:Number} = copy(_lazy_transpose(A))


# Each work group transposes a tile of 32 x 32 elements, with 32 x 8 work items. Consecutive
# work groups cover the tiles along the shorter dimension first, so that the short columns of
# the input or output are read or written as a whole while they are in the cache.
@kernel function transpose_kernel_impl(
    Y::AbstractMatrix{T_out},
    @Const(X::AbstractMatrix{T_in}),
    n_tiles_fast::Int, i_fast::Bool
) where {T_out<:Number,T_in<:Number}
    group = @index(Group, Linear)
    local_i, local_j = @index(Local, NTuple)

    # Padded, so that reading columns of the tile causes no shared memory bank conflicts:
    buf = @localmem T_in (33, 32)

    tile_fast = (group - 1) % n_tiles_fast + 1
    tile_slow = (group - 1) ÷ n_tiles_fast + 1
    tile_i, tile_j = i_fast ? (tile_fast, tile_slow) : (tile_slow, tile_fast)

    offs_i = (tile_i-1) * 32 + first(axes(X,1)) - 1
    offs_j = (tile_j-1) * 32 + first(axes(X,2)) - 1
    for k in 0:3
        j = local_j + 8 * k
        @inbounds if offs_i + local_i <= last(axes(X,1)) && offs_j + j <= last(axes(X,2))
            buf[local_i, j] = X[offs_i + local_i, offs_j + j]
        end
    end

    @synchronize

    offs_i = (tile_j-1) * 32 + first(axes(Y,1)) - 1
    offs_j = (tile_i-1) * 32 + first(axes(Y,2)) - 1
    for k in 0:3
        j = local_j + 8 * k
        @inbounds if offs_i + local_i <= last(axes(Y,1)) && offs_j + j <= last(axes(Y,2))
            Y[offs_i + local_i, offs_j + j] = buf[j, local_i]
        end
    end
end

function _ka_nonlazy_transpose!(
    Y::AbstractMatrix{<:Real},
    X::AbstractMatrix{<:Real},
)
    @argcheck !Base.mightalias(X, Y)
    @argcheck axes(Y) == reverse(axes(X))

    backend = _ka_get_backend(X)

    n_tiles_i = div(size(X, 1) + 31, 32)
    n_tiles_j = div(size(X, 2) + 31, 32)
    i_fast = size(X, 1) <= size(X, 2)
    n_tiles_fast = i_fast ? n_tiles_i : n_tiles_j

    kernel! = transpose_kernel_impl(backend, (32, 8))
    kernel_ret = kernel!(Y, X, n_tiles_fast, i_fast, ndrange = (32 * n_tiles_i * n_tiles_j, 8))
    _ka_synchronize(kernel!, kernel_ret)

    return Y
end

function _nonlazy_transpose(X::GPUArraysCore.AbstractGPUArray{<:Number,2})
    Y = similar(X, reverse(size(X)));
    _ka_nonlazy_transpose!(Y, X)
end
