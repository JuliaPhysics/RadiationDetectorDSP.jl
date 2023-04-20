# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


_lazy_transpose(A::AbstractMatrix{T}) where {T<:Number} = transpose(A)
_nonlazy_transpose(A::AbstractMatrix{T}) where {T<:Number} = copy(_lazy_transpose(A))


@kernel function transpose_kernel_impl(
    Y::AbstractMatrix{T_out},
    X::AbstractMatrix{T_in},
) where {T_out<:Number,T_in<:Number}
    tile_size = @uniform @groupsize()[1] # Must be static
    # Sanity check:
    #@assert tile_size == @groupsize()[2]

    tile_i, tile_j = @index(Group, NTuple)
    local_i, local_j = @index(Local, NTuple)

    buf = @localmem T_in (tile_size, tile_size)

    global_i = (tile_i-1) * tile_size + first(axes(X,1)) - 1 + local_i
    global_j = (tile_j-1) * tile_size + first(axes(X,2)) - 1 + local_j

    # Sanity check:
    #in_global_i, in_global_j = @index(Global, NTuple)
    #@assert global_i == in_global_i
    #@assert global_j == in_global_j

    @inbounds if global_i in axes(X,1) && global_j in axes(X,2)
        buf[local_j, local_i] = X[global_i, global_j]
    end

    @synchronize

    global_i = (tile_j-1) * tile_size + first(axes(Y,1)) - 1 + local_i
    global_j = (tile_i-1) * tile_size + first(axes(Y,2)) - 1 + local_j

    @inbounds if global_i in axes(Y,1) && global_j in axes(Y,2)
        Y[global_i, global_j] = buf[local_i, local_j]
    end
end

function _ka_nonlazy_transpose!(
    Y::AbstractMatrix{<:Real},
    X::AbstractMatrix{<:Real},
)
    @argcheck !Base.mightalias(X, Y)
    @argcheck axes(Y) == reverse(axes(X))

    backend = _ka_get_backend(X)

    tile_size = 16
    n_tiles_i = div(size(X, 1) + tile_size - 1, tile_size)
    n_tiles_j = div(size(X, 2) + tile_size - 1, tile_size)
    ndrange = n_tiles_i * tile_size, n_tiles_j * tile_size

    kernel! = transpose_kernel_impl(backend, (tile_size, tile_size))
    kernel_ret = kernel!(Y, X, ndrange = ndrange) 
    _ka_synchronize(kernel!, kernel_ret)

    return Y
end

function _nonlazy_transpose(X::AbstractGPUArray{<:Number,2})
    Y = similar(X, reverse(size(X)));
    _ka_nonlazy_transpose!(X, Y)
end
