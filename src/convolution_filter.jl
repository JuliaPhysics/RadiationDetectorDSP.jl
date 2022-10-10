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
end

export ConvolutionFilter


function fltinstance(flt::ConvolutionFilter{DirectConvolution}, si::SamplingInfo{T}) where T
    reverse_h = reverse(T.(flt.coeffs)) # ToDo: Optimize
    DirectConvFilterInstance(reverse_h, _smpllen(si))
end

function fltinstance(flt::ConvolutionFilter{FFTConvolution}, si::SamplingInfo{T}) where T
    rfft_h = rfft(T.(flt.coeffs)) # ToDo: Optimize
    FFTConvFilterInstance(rfft_h)
end


DSP.FIRFilter(flt::ConvolutionFilter) = DSP.FIRFilter(flt.coeffs)



abstract type AbstractConvFilterInstance{T<:Real} <: AbstractRadSigFilterInstance{LinearFiltering} end

flt_output_smpltype(fi::AbstractConvFilterInstance) = flt_input_smpltype(fi)
flt_input_smpltype(fi::AbstractConvFilterInstance{T}) where T = T
flt_output_length(fi::AbstractConvFilterInstance) = flt_input_length(fi) - filterlen(fi) + 1
flt_input_length(fi::AbstractConvFilterInstance) = fi.n_input

function flt_output_time_axis(fi::AbstractConvFilterInstance, time::AbstractVector{<:RealQuantity})
    valid_range = (firstindex(time) + fi.n_filter - 1):lastindex(time)
    time[valid_range]
end



struct DirectConvFilterInstance{T<:Real,TV<:AbstractVector{T}} <: AbstractConvFilterInstance{T}
    reverse_h::TV
    n_input::Int
end

_filterlen(fi::DirectConvFilterInstance) = fi.n_input


@inline function rdfilt!(y::AbstractVector{T}, fi::DirectConvFilterInstance{T}, x::AbstractVector{T}) where {T<:Real}
    rh = fi.reverse_h

    @assert firstindex(y) == firstindex(x) == firstindex(rh) && lastindex(x) >= lastindex(rh)
    @assert lastindex(y) == lastindex(x) - (lastindex(rh) - firstindex(rh))

    for i in eachindex(y)
        y[i] = 0
        @simd for j in eachindex(rh)
            y[i] = fma(rh[j], x[i+j-1], y[i])
        end
    end
    return y
end



struct FFTConvFilterInstance{T<:Real,TV<:AbstractVector{Complex{T}}} <: AbstractConvFilterInstance{T}
    rfft_h::TV
end

_filterlen(fi::FFTConvFilterInstance) = size(fi.rfft_h, 1)


# # true starting from i = non-zero-length of h:
# filt(FIRFilter(h), x) .≈ irfft(rfft(h) .* rfft(x), length(h))

# # Multi-waveform:
# irfft(rfft(H, 1) .* rfft(X, 1), size(H, 1), 1)


function rdfilt(fi::FFTConvFilterInstance, x::AbstractVector{<:Real})
    y_ext = irfft(rfft(x) .* fi.rfft_h, size(x, 1))
    valid_range = (firstindex(y) + fi.n_filter - 1):lastindex(y)
    y_ext[valid_range]
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
    Y = irfft(rfft(X, 1) .* fi.rfft_h, size(X, 1), 1)
    valid_range = (firstindex(y) + fi.n_filter - 1):lastindex(y)
    flat_output = Y[valid_range, :]
    ArrayOfSimilarArrays{T_out,1,N}(flat_output)
end

function bc_rdfilt!(
    outputs::ArrayOfSimilarArrays{<:Real,1},
    fi::FFTConvFilterInstance,
    inputs::ArrayOfSimilarArrays{<:Real,1}
)
    flatview(outputs) .= flatview(bc_rdfilt(fi, inputs))
    return outputs
end
