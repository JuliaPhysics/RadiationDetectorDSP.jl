# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    abstract type ConvolutionMethod

Indended as a type parameter to designate the behavior of a filter
as linear or nonlinear.

Subtypes are [`DirectConvolution`](@ref) and [`FFTConvolution`](@ref).
"""
abstract type ConvolutionMethod end


"""
    DirectConvolution() isa ConvolutionMethod

Compute filter convolutions directly, without FFT.
"""
struct DirectConvolution <: ConvolutionMethod end
export DirectConvolution


"""
    FFTConvolution() isa ConvolutionMethod

Compute filter convolutions via FFT.
"""
struct FFTConvolution <: ConvolutionMethod end
export FFTConvolution



"""
    struct ConvolutionFilter{T<:RealQuantity} <: AbstractRadFIRFilter

A [FIR filter](https://en.wikipedia.org/wiki/Finite_impulse_response) defined
by it's filter taps, applied via convolution with the input signal.

Constructors:

* ```$(FUNCTIONNAME)(fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct ConvolutionFilter{C<:ConvolutionMethod,T<:RealQuantity,TV<:AbstractVector{T}} <: AbstractRadFIRFilter
    "Convolution method"
    method::C

    "Filter taps"
    coeffs::TV

    "Time axis offset"
    offset::Int = 0 
end

ConvolutionFilter(method::ConvolutionMethod, coeffs::AbstractVector{<:RealQuantity}) = ConvolutionFilter(method, coeffs, 0)

export ConvolutionFilter

Adapt.adapt_structure(to, flt::ConvolutionFilter) = ConvolutionFilter(flt.method, Adapt.adapt_structure(to, flt.coeffs))

function fltinstance(flt::ConvolutionFilter{DirectConvolution}, si::SamplingInfo{T}) where T
    reverse_h = reverse(T.(flt.coeffs)) # ToDo: Optimize
    DirectConvFilterInstance(reverse_h, _smpllen(si), flt.offset)
end

function fltinstance(flt::ConvolutionFilter{FFTConvolution}, si::SamplingInfo{T}) where T
    n_filter = size(flt.coeffs, 1)
    n_input = length(si.axis)
    @argcheck 1 <= n_filter <= n_input
    # FFTs are fast for lengths with only small prime factors. Zero-padding the input to n_fft
    # leaves the output samples n_filter:n_input free of circular wrap-around.
    n_fft = nextprod((2, 3, 5), n_input)
    h_ext = similar(flt.coeffs, T, n_fft)
    fill!(h_ext, zero(eltype(h_ext)))
    copy!(view(h_ext, firstindex(h_ext):(firstindex(h_ext) + n_filter - 1)), flt.coeffs)
    # includes the 1/n_fft normalization of the inverse transform (applied with brfft)
    rfft_h = rfft(h_ext) ./ n_fft
    FFTConvFilterInstance(rfft_h, n_filter, n_input, n_fft, flt.offset)
end


DSP.FIRFilter(flt::ConvolutionFilter) = DSP.FIRFilter(flt.coeffs)



abstract type AbstractConvFilterInstance{T<:Real} <: AbstractRadSigFilterInstance{LinearFiltering} end

flt_output_smpltype(fi::AbstractConvFilterInstance) = flt_input_smpltype(fi)
flt_input_smpltype(fi::AbstractConvFilterInstance{T}) where T = T
flt_output_length(fi::AbstractConvFilterInstance) = flt_input_length(fi) - _filterlen(fi) + 1
flt_input_length(fi::AbstractConvFilterInstance) = fi.n_input

function flt_output_time_axis(fi::AbstractConvFilterInstance, time::AbstractVector{<:RealQuantity})
    valid_range = (firstindex(time) + _filterlen(fi) - 1 + fi.offset):(lastindex(time) + fi.offset)
    time[valid_range]
end



struct DirectConvFilterInstance{T<:Real,TV<:AbstractVector{T}} <: AbstractConvFilterInstance{T}
    reverse_h::TV
    n_input::Int
    offset::Int
end

Adapt.adapt_structure(to, fi::DirectConvFilterInstance) = DirectConvFilterInstance(adapt(to, fi.reverse_h), fi.n_input, fi.offset)

_filterlen(fi::DirectConvFilterInstance) = size(fi.reverse_h, 1)

@inline function rdfilt!(y::AbstractVector{T}, fi::DirectConvFilterInstance{T}, x::AbstractVector{T}) where {T<:Real}
    rh = fi.reverse_h

    @assert firstindex(y) == firstindex(x) == firstindex(rh) && lastindex(x) >= lastindex(rh)
    @assert lastindex(y) == lastindex(x) - (lastindex(rh) - firstindex(rh))

    @inbounds for i in eachindex(y)
        acc = zero(T)
        @simd for j in eachindex(rh)
            acc = fma(rh[j], x[i+j-1], acc)
        end
        y[i] = acc
    end
    return y
end

function bc_rdfilt!(
    outputs::ArrayOfSimilarVectors{<:RealQuantity},
    fi::DirectConvFilterInstance,
    inputs::ArrayOfSimilarVectors{<:RealQuantity}
)
    _ka_bc_rdfilt!(outputs, fi, inputs)
end


"""
    FFTConvFilterInstance

Filter instance of a `ConvolutionFilter{FFTConvolution}` for inputs of length `n_input`.

`rfft_h` is the real FFT of the filter taps zero-padded to length `n_fft >= n_input`, divided
by `n_fft`. The input is zero-padded to `n_fft` as well, so `brfft(rfft(x_ext) .* rfft_h, n_fft)`
is the circular convolution of length `n_fft`, and its samples `n_filter:n_input` equal the
linear convolution of `x` with the taps.
"""
struct FFTConvFilterInstance{T<:Real,TV<:AbstractVector{Complex{T}}} <: AbstractConvFilterInstance{T}
    rfft_h::TV
    n_filter::Int
    n_input::Int
    n_fft::Int
    offset::Int
end

Adapt.adapt_structure(to, fi::FFTConvFilterInstance) = FFTConvFilterInstance(adapt(to, fi.rfft_h), fi.n_filter, fi.n_input, fi.n_fft, fi.offset)

_filterlen(fi::FFTConvFilterInstance) = fi.n_filter


function rdfilt(fi::FFTConvFilterInstance{T}, x::AbstractVector{<:Real}) where T
    @argcheck length(x) == fi.n_input
    x_ext = fill!(similar(x, T, fi.n_fft), zero(T))
    copyto!(x_ext, firstindex(x_ext), x, firstindex(x), length(x))
    y_ext = brfft(rfft(x_ext) .* fi.rfft_h, fi.n_fft)
    y_ext[(firstindex(y_ext) + fi.n_filter - 1):(firstindex(y_ext) + fi.n_input - 1)]
end

function rdfilt!(y::AbstractVector{T}, fi::FFTConvFilterInstance, x::AbstractVector{T}) where {T<:Real}
    y .= rdfilt(fi, x)
end


function bc_rdfilt(
    fi::FFTConvFilterInstance,
    inputs::ArrayOfSimilarArrays{<:Real,1,N}
) where N
    T_out = flt_output_smpltype(fi)
    X = flatview(inputs)
    Y = similar(X, T_out, (flt_output_length(fi), Base.tail(size(X))...))
    _fft_conv_columns!(reshape(Y, size(Y, 1), :), fi, reshape(X, size(X, 1), :))
    ArrayOfSimilarArrays{T_out,1,N}(Y)
end

function bc_rdfilt!(
    outputs::ArrayOfSimilarArrays{<:Real,1},
    fi::FFTConvFilterInstance,
    inputs::ArrayOfSimilarArrays{<:Real,1}
)
    X, Y = flatview(inputs), flatview(outputs)
    @argcheck Base.tail(size(Y)) == Base.tail(size(X))
    _fft_conv_columns!(reshape(Y, size(Y, 1), :), fi, reshape(X, size(X, 1), :))
    return outputs
end


# Convolution of each column of X, kept samples written to the columns of Y.
function _fft_conv_columns!(Y::AbstractMatrix, fi::FFTConvFilterInstance{T}, X::AbstractMatrix) where T
    @argcheck size(X, 1) == fi.n_input
    @argcheck size(Y, 1) == flt_output_length(fi)
    @argcheck axes(Y, 2) == axes(X, 2)
    if get_compute_unit(X) isa CPUnit && get_compute_unit(Y) isa CPUnit
        _fft_conv_columns_blocked!(Y, fi, X)
    else
        adapted_fi = adapt(get_compute_unit(X), fi)
        X_ext = fill!(similar(X, T, (fi.n_fft, size(X, 2))), zero(T))
        X_ext[firstindex(X_ext, 1):(firstindex(X_ext, 1) + fi.n_input - 1), :] .= X
        Y_ext = brfft(rfft(X_ext, 1) .* adapted_fi.rfft_h, fi.n_fft, 1)
        Y .= view(Y_ext, (firstindex(Y_ext, 1) + fi.n_filter - 1):(firstindex(Y_ext, 1) + fi.n_input - 1), :)
    end
    return Y
end

# On the CPU, the columns are transformed in blocks of at most `blocksize` columns, one task per
# group of blocks with its own buffers and FFTW plans, so the temporary arrays stay small and all
# threads are used.
function _fft_conv_columns_blocked!(Y::AbstractMatrix, fi::FFTConvFilterInstance{T}, X::AbstractMatrix; blocksize::Int = 256) where T
    n_cols = size(X, 2)
    n_cols == 0 && return Y
    bs = min(blocksize, n_cols)
    blocks = collect(Iterators.partition(axes(X, 2), bs))
    n_tasks = min(Threads.nthreads(), length(blocks))
    groups = collect(Iterators.partition(eachindex(blocks), cld(length(blocks), n_tasks)))
    rfft_h = convert(Vector{Complex{T}}, fi.rfft_h)
    in_rows = firstindex(X, 1):lastindex(X, 1)
    kept_rows = fi.n_filter:fi.n_input
    # Rows n_input+1:n_fft of X_ext stay zero.
    bufs = map(groups) do _
        X_ext = zeros(T, fi.n_fft, bs)
        F = Matrix{Complex{T}}(undef, fi.n_fft ÷ 2 + 1, bs)
        Y_ext = Matrix{T}(undef, fi.n_fft, bs)
        (; X_ext, F, Y_ext, p_fwd = plan_rfft(X_ext, 1), p_bwd = plan_brfft(F, fi.n_fft, 1))
    end
    tasks = map(groups, bufs) do group, buf
        Threads.@spawn for b in group
            cols = blocks[b]
            n = length(cols)
            view(buf.X_ext, 1:fi.n_input, 1:n) .= view(X, in_rows, cols)
            mul!(buf.F, buf.p_fwd, buf.X_ext)
            buf.F .*= rfft_h
            mul!(buf.Y_ext, buf.p_bwd, buf.F)
            view(Y, :, cols) .= view(buf.Y_ext, kept_rows, 1:n)
        end
    end
    foreach(wait, tasks)
    return Y
end
