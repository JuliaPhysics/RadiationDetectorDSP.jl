# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


_acc_add(acc, x) = acc + x
_acc_add(acc::NTuple{N}, x::NTuple{N}) where N = map(_acc_add, acc, x)
_acc_add(acc::NamedTuple{names}, x::NamedTuple{names}) where names = NamedTuple{names}(map(_acc_add, values(acc), values(x)))

_acc_weighted_add(acc, x, weight::Number) = muladd(weight, x, acc)
_acc_weighted_add(acc::NTuple{N}, x::NTuple{N}, weight::Number) where N = ntuple(i -> _acc_weighted_add(acc[i], x[i], weight), Val(N))
_acc_weighted_add(acc::NamedTuple{names}, x::NamedTuple{names}, weight::Number) where names = NamedTuple{names}(_acc_weighted_add(values(acc), values(x), weight))


@kernel function wf_map_sum_kernel_impl(
    f_presum,
    f_postsum,
    Y::AbstractVector{<:Real},
    X::AbstractMatrix{<:Real},
    I_start::Union{Integer,AbstractVector{<:Integer}},
    n::Integer,
    w_first::Number,
)
    # i indexes samples, j indexes waveforms

    @uniform T_in = _floattype(eltype(X))

    n_workers = @uniform @groupsize()[1] # Must be static
    worker_group, = @index(Group, NTuple)
    worker, = @index(Local, NTuple)

    @uniform tile_size = n_workers
    @uniform n_tiles = div(n + n_workers - 1, n_workers)
    @uniform w_one = one(typeof(w_first))
    @uniform w_last = w_one - w_first
    @uniform acc_initval = f_presum(-zero(T_in)) * w_one

    worker == 1 && worker_group ==1 && @info "DEBUG" tile_size n_tiles

    buf = @localmem T_in (tile_size, tile_size)

    acc = @private typeof(acc_initval) 1
    #!!!@inbounds 
    acc[1] = acc_initval # necessary?

    for tile in 1:n_tiles
        # Note: Indexing is X[i, j] but buf[j, i]

        nan_value = convert(T_in, NaN)
        global_i_offset = (tile-1) * tile_size + first(axes(X,1)) - 1
        global_j_offset = (worker_group-1) * tile_size + first(axes(X,2)) - 1
        local_i = worker
        global_i_base = global_i_offset + local_i

        #!!!@fastmath @inbounds @simd
        for local_j in Base.OneTo(tile_size)
            buf[local_j, local_i] = nan_value
        end

        #!!!@fastmath @inbounds @simd
        for local_j in Base.OneTo(tile_size)
            global_j = global_j_offset + local_j
            if global_j in axes(X,2)
                i_start_offset = _kbc_getindex(I_start, global_j) - first(axes(X,1))
                global_i = global_i_base + i_start_offset
                #@info "DEBUG $((;local_j, global_j, local_i, global_i))"
                if global_i in axes(X,1)
                    buf[local_j, local_i] = convert(T_in, X[global_i, global_j])
                end
            end
        end

        @synchronize

        #worker == 1 && @info "DEBUG" buf

        global_j_offset = (worker_group-1) * tile_size + first(axes(X,2)) - 1
        local_j = worker
        global_j = global_j_offset + local_j
        n_i = min(n - (tile-1) * tile_size, tile_size)
        tmp_acc::typeof(acc_initval) = acc_initval
        if global_j in axes(X,2)
            #@info "DEBUG $((;tile, tile_size, n, n_i))"
            if n_i > 1
                tmp_acc = let local_i = 1
                    #@info "DEBUG $((;local_j, local_i))"
                    f_x = f_presum(buf[local_j, local_i])
                    weight = ifelse(tile == 1, w_first, w_one)
                    _acc_weighted_add(tmp_acc, f_x, weight)
                end
            end
            #!!!@fastmath @inbounds @simd
            for local_i in 2:(n_i - 1)
                #@info "DEBUG $((;local_j, local_i))"
                f_x = f_presum(buf[local_j, local_i])
                tmp_acc = _acc_add(tmp_acc, f_x)
            end
            tmp_acc = let local_i = n_i
                #@info "DEBUG $((;local_j, local_i))"
                f_x = f_presum(buf[local_j, local_i])
                weight = ifelse(tile == n_tiles, w_last, w_one)
                _acc_weighted_add(tmp_acc, f_x, weight)
            end
            acc[1] += tmp_acc
        end
    end

    # Sanity check:
    global_j, = @index(Global, NTuple)
    @assert global_j == (worker_group-1) * tile_size + first(axes(X,2)) - 1 + worker

    if global_j in axes(X,2)
    #!!!@inbounds
    Y[global_j] = f_postsum(acc[1])
    end
end
