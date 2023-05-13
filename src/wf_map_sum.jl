# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


_acc_add(acc, x) = acc + x
_acc_add(acc::NTuple{N}, x::NTuple{N}) where N = map(_acc_add, acc, x)
_acc_add(acc::NamedTuple{names}, x::NamedTuple{names}) where names = NamedTuple{names}(map(_acc_add, values(acc), values(x)))

_acc_weighted_add(acc, x, weight::Number) = muladd(weight, x, acc)
_acc_weighted_add(acc::NTuple{N}, x::NTuple{N}, weight::Number) where N = ntuple(i -> _acc_weighted_add(acc[i], x[i], weight), Val(N))
_acc_weighted_add(acc::NamedTuple{names}, x::NamedTuple{names}, weight::Number) where names = NamedTuple{names}(_acc_weighted_add(values(acc), values(x), weight))

_deep_mul(x, weight::Number) = weight * x
_deep_mul(x::NTuple{N}, weight::Number) where N = ntuple(i -> _deep_mul(x[i], weight), Val(N))
_deep_mul(x::NamedTuple{names}, weight::Number) where names = NamedTuple{names}(_deep_mul(values(x), weight))


@inline function _wf_map_sum_single(
    f_presum, acc, x::AbstractRange{<:RealQuantity}, Y_buf::AbstractMatrix{<:Real}
    weight_last::Real, is_first_tile::Bool, is_last_tile::Bool,
    n_i::Integer, local_j::Integer,
)
    weight_first = one(weight_last) - weight_last
    
    tmp_acc::typeof(acc) = acc
    if n_i > 1
        tmp_acc = let local_i = 1
            #@debug "wf_map_sum_kernel_impl summing: $((;local_j, local_i))"
            f_x = f_presum(Y_buf[local_j, local_i])
            weight = ifelse(is_first_tile, weight_first, one(weight_first))
            _acc_weighted_add(tmp_acc, f_x, weight)
        end
    end
    @fastmath @inbounds @simd for local_i in 2:(n_i - 1)
        #@debug "wf_map_sum_kernel_impl summing: $((;local_j, local_i))"
        f_x = f_presum(Y_buf[local_j, local_i])
        tmp_acc = _acc_add(tmp_acc, f_x)
    end
    tmp_acc = let local_i = n_i
        #@debug "wf_map_sum_kernel_impl summing: $((;local_j, local_i))"
        f_x = f_presum(Y_buf[local_j, local_i])
        weight = ifelse(is_last_tile, weight_last, one(weight_last))
        _acc_weighted_add(tmp_acc, f_x, weight)
    end
    return tmp_acc
end


@kernel function wf_map_sum_kernel_impl(
    f_presum,
    f_postsum,
    Z::AbstractVector,
    X::Union{AbstractRange{<:RealQuantity},AbstractVector{<:AbstractRange{<:RealQuantity}}},
    Y::AbstractArray{<:Real,2},
    I_start::Union{Real,AbstractVector{<:Real}},
    n::Integer
)
    @uniform T_in = _floattype(eltype(Y))
    @uniform n_eff = n + 1
    @uniform weight_one = one(_floattype(eltype(Y)))

    # Assumes that f_presum(one(T_in) is not NaN or Inf:
    @uniform acc_initval = _deep_mul(f_presum(one(T_in)), zero(weight_one))

    onebased_j, = @index(Global, NTuple)
    global_j = onebased_j + first(axes(Y,2)) - 1
    if global_j in axes(Y,2)
        i_start_real = _kbc_getindex(I_start, global_j)
        weight_last = i_start_real - floor(i_start_real)
        i_start = floor(Int, i_start_real)
    
        Y_buf = transpose(view(Y, i_start:i_start+n_eff-1, global_j))
        tmp_sum = _wf_map_sum_single(
            f_presum, acc_initval, Y_buf,
            weight_last, true, true, n_eff, 1
        )
        @inbounds Z[global_j] = f_postsum(tmp_sum)
    end
end


@kernel function wf_map_sum_kernel_impl(
    f_presum,
    f_postsum,
    Z::AbstractVector,
    X::Union{AbstractRange{<:RealQuantity},AbstractVector{<:AbstractRange{<:RealQuantity}}},
    Y::GPULikeArray{<:Real, 2},
    I_start::Union{Real,AbstractVector{<:Real}},
    n::Integer
)
    # i indexes samples, j indexes waveforms

    @uniform T_in = _floattype(eltype(Y))
    weight_one = one(_floattype(eltype(Y)))
    n_workers = @uniform @groupsize()[1] # Must be static
    worker_group, = @index(Group, NTuple)
    worker, = @index(Local, NTuple)

    @uniform tile_size = n_workers
    @uniform n_eff = n + 1
    @uniform n_tiles = div(n_eff + n_workers - 1, n_workers)
    @uniform weight_one = one(_floattype(eltype(Y)))

    # Assumes that f_presum(one(T_in) is not NaN or Inf:
    @uniform acc_initval = _deep_mul(f_presum(one(T_in)), zero(weight_one))

    #worker == 1 && worker_group ==1 && @debug "wf_map_sum_kernel_impl:" tile_size n_tiles

    Y_buf = @localmem T_in (tile_size, tile_size)

    acc = @private typeof(acc_initval) 1
    @inbounds acc[1] = acc_initval # necessary?

    for tile in 1:n_tiles
        let
            # Note: Indexing is Y[i, j] but Y_buf[j, i]

            nan_value = convert(T_in, NaN)
            global_i_offset = (tile-1) * tile_size + first(axes(Y,1)) - 1
            global_j_offset = (worker_group-1) * tile_size + first(axes(Y,2)) - 1
            local_i = worker
            global_i_base = global_i_offset + local_i

            @fastmath @inbounds @simd for local_j in Base.OneTo(tile_size)
                Y_buf[local_j, local_i] = nan_value
            end

            @fastmath @inbounds @simd for local_j in Base.OneTo(tile_size)
                global_j = global_j_offset + local_j
                if global_j in axes(Y,2)
                    i_start_real = _kbc_getindex(I_start, global_j)
                    i_start = floor(Int, i_start_real)
                    i_start_offset = i_start - first(axes(Y,1))
                    global_i = global_i_base + i_start_offset
                    #@debug "wf_map_sum_kernel_impl buffering: $((;local_j, global_j, local_i, global_i))"
                    if global_i in axes(Y,1)
                        Y_buf[local_j, local_i] = convert(T_in, Y[global_i, global_j])
                    end
                end
            end
        end

        @synchronize

        let
            #worker == 1 && @debug "wf_map_sum_kernel_impl" Y_buf

            global_j_offset = (worker_group-1) * tile_size + first(axes(Y,2)) - 1
            local_j = worker
            global_j = global_j_offset + local_j
            n_i = min(n_eff - (tile-1) * tile_size, tile_size)
            if global_j in axes(Y,2)
                i_start_real = _kbc_getindex(I_start, global_j)
                weight_last = i_start_real - floor(i_start_real)
                #@debug "wf_map_sum_kernel_impl: $((;tile, tile_size))"
                tmp_sum = _wf_map_sum_single(
                    f_presum, acc_initval, Y_buf,
                    weight_last, tile == 1, tile == n_tiles, n_i, local_j
                )
                acc[1] = _acc_add(acc[1], tmp_sum)
            end
        end
    end

    # Sanity check:
    global_j, = @index(Global, NTuple)
    @assert global_j == (worker_group-1) * tile_size + first(axes(Y,2)) - 1 + worker

    if global_j in axes(Y,2)
        @inbounds Z[global_j] = f_postsum(acc[1])
    end
end


function run_wf_map_sum_kernel!(
    f_presum,
    f_postsum,
    Z::AbstractVector,
    Y::AbstractMatrix{<:Real},
    I_start::Union{Real,AbstractVector{<:Real}},
    n::Integer
)
    backend = _ka_get_backend(Y)
    n_workers = 32
    kernel! = wf_map_sum_kernel_impl(backend, (n_workers,))
    n_waveforms = length(axes(Y, 2))
    n_tiles = div(n_waveforms + n_workers - 1, n_workers)
    kernel_ret = kernel!(f_presum, f_postsum, Z, Y, I_start, n, ndrange=(n_tiles * n_workers,)) 
    _ka_synchronize(kernel!, kernel_ret)
    return Z
end


function _construct_wf_map_sum_output(f_presum, f_postsum, Y::AbstractArray{<:Real})
    # T_presum = Core.Compiler.return_type(f_presum, Tuple{eltype(Y)})
    # T_postsum = Core.Compiler.return_type(f_postsum, Tuple{T_presum})
    T_postsum = typeof(f_postsum(f_presum(one(eltype(Y)))))
    _similar_maybe_structarray(Y, T_postsum, (size(Y, 2),))
end
