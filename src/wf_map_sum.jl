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
    I_start::AbstractVector{<:Integer},
    n::Integer,
    w_first::Number,
)
    # i indexes samples, j indexes waveforms

    @uniform T_in = _floattype(eltype(X))

    n_workers = @uniform @groupsize()[1] # Must be static
    worker_group = @index(Group, NTuple)[1]
    worker = @index(Local, NTuple)[1]

    @uniform tile_size = n_workers
    @uniform n_tiles = div(n + n_workers - 1, n_workers)
    @uniform w_one = one(typeof(w_first))
    @uniform w_last = w_one - w_first
    @uniform acc_initval = f_presum(-zero(T_in)) * w_one

    buf = @localmem T_in (tile_size, tile_size)

    # Sanity check:
    global_j, = @index(Global, NTuple)
    @assert global_j == (worker_group-1) * tile_size + first(axes(X,2)) - 1 + worker

    acc = @private typeof(acc_initval) acc_initval
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
            global_j = global_j_offset + local_j
            i_start_offset = I_start[global_j] - first(axes(X,1))
            global_i = global_i_base + i_start_offset
            if global_i in axes(X,1) && global_j in axes(X,2)
                buf[local_j, local_i] = convert(T_in, X[global_i, global_j])
            else
                buf[local_j, local_i] = nan_value
            end
        end

        @synchronize

        tmp_acc::eltype(acc_initval) = acc_initval
        tmp_acc = let i = 1
            f_x = f_presum(buf[worker, i])
            weight = ifelse(tile == 1, w_first, w_one)
            _acc_weighted_add(tmp_acc, f_x, weight)
        end
        #!!!@fastmath @inbounds @simd
        for i in 2:(tile_size - 1)
            f_x = f_presum(buf[worker, i])
            tmp_acc = _acc_add(tmp_acc, f_x)
        end
        tmp_acc = let i = tile_size
            f_x = f_presum(buf[worker, i])
            weight = ifelse(tile == n_tiles, w_last, w_one)
            _acc_weighted_add(tmp_acc, f_x, weight)
        end
        acc[1] += tmp_acc
    end

    #!!!@inbounds
    Y[global_i] = f_postsum(acc[1])
end
