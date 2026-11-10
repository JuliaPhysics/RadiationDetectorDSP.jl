# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    struct TrapezoidalChargeFilter <: AbstractRadFIRFilter

Filter that responds to a step signal with a trapezoidal pulse.

The filter is equivalent to two moving averages separated by a gap.

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)

A sharp step on the input will result in a trapezoid with rise time and fall
time `avgtime` and a flat top of length `gaptime`.
"""
Base.@kwdef struct TrapezoidalChargeFilter{
    T <: RealQuantity
} <: AbstractRadFIRFilter
    "pre-rise averaging time"
    avgtime::T

    "gap time"
    gaptime::T

    "post-rise averaging time"
    avgtime2::T = avgtime
end

TrapezoidalChargeFilter(avgtime::RealQuantity, gaptime::RealQuantity) = TrapezoidalChargeFilter(avgtime, gaptime, avgtime)

export TrapezoidalChargeFilter

Adapt.adapt_structure(to, flt::TrapezoidalChargeFilter) = flt

function ConvolutionFilter(flt::TrapezoidalChargeFilter)
    navg = Int(flt.avgtime)
    navg2 = Int(flt.avgtime2)
    ngap = Int(flt.gaptime)
    norm_factor = inv(navg)
    norm_factor2 = inv(navg2)
    T = typeof(norm_factor)
    params = zeros(T, navg + ngap + navg2)
    fill!(view(params, firstindex(params):(firstindex(params)+navg-1)), +norm_factor)
    fill!(view(params, firstindex(params)+navg+ngap:lastindex(params)), -norm_factor2)
    ConvolutionFilter(FFTConvolution(), params)
end


function fltinstance(flt::TrapezoidalChargeFilter, si::SamplingInfo{T}) where {T <: RealQuantity} 
    delta_t = step(si.axis)

    navg  = round(Int, uconvert(NoUnits, flt.avgtime  / delta_t))
    ngap  = round(Int, uconvert(NoUnits, flt.gaptime  / delta_t))
    navg2 = round(Int, uconvert(NoUnits, flt.avgtime2 / delta_t))

    navg >= 1 && ngap >= 0 && navg2 >= 1 || throw(ArgumentError("Require navg/navg2 >= 1 and ngap >= 0"))
    (navg + ngap + navg2) <= _smpllen(si) || throw(ArgumentError("filter must not be longer than input"))

    TrapezoidalChargeFilterInstance{T}(navg, ngap, navg2, _smpllen(si))
end


struct TrapezoidalChargeFilterInstance{T<:RealQuantity} <: AbstractRadSigFilterInstance{LinearFiltering}
    navg::Int
    ngap::Int
    navg2::Int
    n_input::Int
end


_filterlen(fi::TrapezoidalChargeFilterInstance) = fi.navg + fi.ngap + fi.navg2

function flt_output_time_axis(fi::TrapezoidalChargeFilterInstance, time::AbstractVector{<:RealQuantity})
    valid_range = (firstindex(time) + _filterlen(fi) - 1):lastindex(time)
    time[valid_range]
end


@inline function rdfilt!(y::AbstractVector{T}, fi::TrapezoidalChargeFilterInstance{U}, x::AbstractVector{U}) where {T<:Real, U<:Real}
    @fastmath begin
        navg  = fi.navg
        navg2 = fi.navg2
        ngap  = fi.ngap

        @assert firstindex(y) == firstindex(x)
        @assert lastindex(y) == lastindex(x) - _filterlen(fi) + 1

        norm_factor = inv(T(navg))
        norm_factor2 = inv(T(navg2))
        acc::T = zero(T)

        offs1::Int = 0
        offs2::Int = offs1 + navg
        offs3::Int = offs2 + ngap
        offs4::Int = offs3 + navg2

        @assert lastindex(y) + offs4 - 1 == lastindex(x)
        
        @inbounds @simd for i in firstindex(x) + offs1:(firstindex(x) + offs2 - 1)
            acc = acc + x[i - offs1] * norm_factor
        end
        @inbounds @simd for i in firstindex(x) + offs2:(firstindex(x) + offs3 - 1)
            acc = acc + (x[i - offs1] - x[i - offs2]) * norm_factor
        end
        @inbounds @simd for i in firstindex(x) + offs3:(firstindex(x) + offs4 - 1)
            acc = acc + (x[i - offs1] - x[i - offs2]) * norm_factor - x[i - offs3] * norm_factor2
        end

        y[firstindex(y)] = acc

        @inbounds @simd for i in firstindex(x) + offs4:(lastindex(x))
            acc = acc + (x[i - offs1] - x[i - offs2]) * norm_factor - (x[i - offs3] - x[i - offs4]) * norm_factor2
            y[i + 1 - offs4] = acc
        end
    end

    return y
end


# The GPU kernel below works with either memory layout. Longer signals are filtered one per work
# item, which needs row-major data on GPUs.
adapt_memlayout(::TrapezoidalChargeFilterInstance, ::KernelAbstractions.CPU, A::AbstractArray{<:Number}) = _col_major(A)
function adapt_memlayout(fi::TrapezoidalChargeFilterInstance, backend::_KA_Backend, A::AbstractArray{<:Number})
    _trapezoid_kernel_fits(fi) && A isa AbstractMatrix ? A : _row_major_on_gpu(backend, A)
end

function bc_rdfilt!(
    outputs::ArrayOfSimilarVectors{<:Number},
    fi::TrapezoidalChargeFilterInstance,
    inputs::ArrayOfSimilarVectors{<:Number}
)
    _trapezoid_bc_rdfilt!(_ka_get_backend(flatview(outputs)), outputs, fi, inputs)
end

_trapezoid_bc_rdfilt!(::KernelAbstractions.CPU, outputs, fi, inputs) = _ka_bc_rdfilt!(outputs, fi, inputs)

function _trapezoid_bc_rdfilt!(backend::_KA_Backend, outputs, fi, inputs)
    X, Y = flatview(inputs), flatview(outputs)
    if _trapezoid_kernel_fits(fi) && X isa AbstractMatrix
        @argcheck axes(X, 2) == axes(Y, 2)
        @argcheck size(X, 1) == fi.n_input && size(Y, 1) == flt_output_length(fi)
        isempty(X) && return outputs
        n_grp = _trapezoid_kernel_group(fi)
        kernel! = _trapezoid_kernel!(backend, 256)
        kernel_ret = kernel!(Y, X, fi.navg, fi.ngap, fi.navg2, n_grp, !(X isa Transpose), !(Y isa Transpose),
                             ndrange = 256 * cld(size(X, 2), n_grp))
        _ka_synchronize(kernel!, kernel_ret)
        return outputs
    else
        _ka_bc_rdfilt!(outputs, fi, inputs)
    end
end

# Signals per work group of the GPU kernel: a power of two, at most eight.
_trapezoid_kernel_group(fi::TrapezoidalChargeFilterInstance) = min(8, prevpow(2, 8192 ÷ fi.n_input))
_trapezoid_kernel_fits(fi::TrapezoidalChargeFilterInstance) = fi.n_input <= 8192

# Lane of the work item in its warp, and the samples of the warp's signal it sums up (positions
# in buf). The chunk length is odd, which avoids shared memory bank conflicts.
@inline function _trapezoid_kernel_chunk(tt::Integer, n_in::Integer)
    w, lane = divrem(tt, 32)
    chunk = cld(n_in, 32)
    chunk += iseven(chunk)
    first_n = w * n_in + lane * chunk + 1
    last_n = w * n_in + min(lane * chunk + chunk, n_in)
    return w, lane, first_n, last_n
end

# Each work group of 256 work items loads n_grp signals into shared memory, computes their
# prefix sums S with one warp per signal, and then all outputs in parallel as
# y[j] = (S[i] - S[i - navg]) / navg - (S[i - navg - ngap] - S[i - navg - ngap - navg2]) / navg2,
# with i = j + navg + ngap + navg2 - 1 and S[0] = 0. The sums are taken over the signal minus its
# first sample, which does not change the result and keeps the sums small. Sample n of signal s
# of the group (both 0-based) is element s * n_in + n + 1 of buf. The memory accesses of the
# work items follow the memory layout of X and Y. Values that differ between work items are
# computed anew after each `@synchronize`, as required by the KernelAbstractions CPU backend.
@kernel function _trapezoid_kernel!(
    Y::AbstractMatrix{<:Real}, @Const(X::AbstractMatrix{<:Real}),
    navg::Int, ngap::Int, navg2::Int, n_grp::Int, col_major_in::Bool, col_major_out::Bool
)
    T = eltype(Y)
    group = @index(Group, Linear)
    tt = @index(Local, Linear) - 1
    buf = @localmem eltype(Y) (8192,)
    totals = @localmem eltype(Y) (32, 8)

    n_in = size(X, 1)
    n_out = size(Y, 1)
    s0 = (group - 1) * n_grp
    n_sig = min(n_grp, size(X, 2) - s0)

    # All loads of a work item are issued before its first store, to keep enough memory
    # accesses in flight.
    if col_major_in
        # The signals of the group are consecutive in linear indexing order:
        n_total = n_sig * n_in
        offs = s0 * n_in
        Base.Cartesian.@nexprs 32 k -> (v_k = tt + 256 * (k - 1) < n_total ? T(@inbounds X[offs + tt + 256 * (k - 1) + 1]) : zero(T))
        Base.Cartesian.@nexprs 32 k -> (tt + 256 * (k - 1) < n_total && (@inbounds buf[tt + 256 * (k - 1) + 1] = v_k; true))
    else
        # n_grp consecutive work items read the signals at the same sample:
        s = tt & (n_grp - 1)
        n_1 = tt ÷ n_grp
        n_step = 256 ÷ n_grp
        Base.Cartesian.@nexprs 32 k -> (v_k = (s < n_sig && n_1 + n_step * (k - 1) < n_in) ? T(@inbounds X[n_1 + n_step * (k - 1) + 1, s0 + s + 1]) : zero(T))
        Base.Cartesian.@nexprs 32 k -> ((s < n_sig && n_1 + n_step * (k - 1) < n_in) && (@inbounds buf[s * n_in + n_1 + n_step * (k - 1) + 1] = v_k; true))
    end
    @synchronize

    # Prefix sums, one chunk of samples per lane:
    w, lane, first_n, last_n = _trapezoid_kernel_chunk(@index(Local, Linear) - 1, n_in)
    acc = zero(T)
    if w < n_sig
        # Not taken from buf, which lane 0 overwrites:
        @inbounds x_1 = T(X[1, s0 + w + 1])
        for n in first_n:last_n
            @inbounds acc += buf[n] - x_1
            @inbounds buf[n] = acc
        end
    end
    @inbounds totals[lane + 1, w + 1] = acc
    @synchronize
    w, lane, first_n, last_n = _trapezoid_kernel_chunk(@index(Local, Linear) - 1, n_in)
    acc_before = zero(T)
    for l in 1:lane
        @inbounds acc_before += totals[l, w + 1]
    end
    if w < n_sig
        for n in first_n:last_n
            @inbounds buf[n] += acc_before
        end
    end
    @synchronize

    tt = @index(Local, Linear) - 1
    norm_factor, norm_factor2 = inv(T(navg)), inv(T(navg2))
    offs4 = navg + ngap + navg2
    if col_major_out
        for s in 0:(n_sig - 1), j in tt:256:(n_out - 1)
            b = s * n_in + j
            @inbounds S_first = j > 0 ? buf[b] : zero(T)
            @inbounds Y[j + 1, s0 + s + 1] = (buf[b + offs4] - buf[b + offs4 - navg]) * norm_factor - (buf[b + navg2] - S_first) * norm_factor2
        end
    else
        s = tt & (n_grp - 1)
        if s < n_sig
            for j in (tt ÷ n_grp):(256 ÷ n_grp):(n_out - 1)
                b = s * n_in + j
                @inbounds S_first = j > 0 ? buf[b] : zero(T)
                @inbounds Y[j + 1, s0 + s + 1] = (buf[b + offs4] - buf[b + offs4 - navg]) * norm_factor - (buf[b + navg2] - S_first) * norm_factor2
            end
        end
    end
end



flt_output_smpltype(fi::TrapezoidalChargeFilterInstance) = float(flt_input_smpltype(fi))
flt_input_smpltype(fi::TrapezoidalChargeFilterInstance{T}) where T = T
flt_output_length(fi::TrapezoidalChargeFilterInstance) = flt_input_length(fi) - _filterlen(fi) + 1
flt_input_length(fi::TrapezoidalChargeFilterInstance) = fi.n_input
