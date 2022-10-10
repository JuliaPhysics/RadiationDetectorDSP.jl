# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    abstract type ConvolutionMethod

Indended as a type parameter to designate the behavior of a filter
as linear or nonlinear.

Subtypes are [`DirectConvolution`](@ref) and [`FFTConvolution`](@ref).
"""
abstract type ConvolutionMethod end


"""
    abstract type DirectConvolution <: ConvolutionMethod

Compute filter convolutions directly, without FFT.
"""
abstract type DirectConvolution <: ConvolutionMethod end
export DirectConvolution


"""
    abstract type FFTConvolution <: ConvolutionMethod

Compute filter convolutions via FFT.
"""
abstract type FFTConvolution <: ConvolutionMethod end
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

    "Offset of coefficients on time axis"
    offset::Int
end

export ConvolutionFilter

## For testing:
#fltparameters(f::ConvolutionFilter) = (b_012 = f.b_012, a_12 = f.a_12)
#fltparameters(f::DSP.Biquad) = (b_012 = SVec(f.b0, f.b1, f.b2), a_12 = SVec(one(f.a1), f.a1, f.a2))


function fltinstance(flt::ConvolutionFilter, si::SamplingInfo{T}) where T
    AbstractConvFilterInstance{T}(flt.b_012, flt.a_12, _smpllen(si))
end


# # true starting from i = non-zero-length of h:
# filt(FIRFilter(h), x) .≈ irfft(rfft(h) .* rfft(x), length(h))

# # Multi-waveform:
# irfft(rfft(H, 1) .* rfft(X, 1), size(H, 1), 1)





function InverseFunctions.inverse(flt::ConvolutionFilter)
    # In direct form 1:
    # y_i[i] = b0 * x_i[i] + b1 * x_i[i-1] + b2 * x_i[i-2] - a1 * y_i[i-1] - a2 * y_i[i-2]
    # x_i[i] = 1/b0 * y_i[i] + a1/b0 * y_i[i-1] + a2/b0 * y_i[i-2] - b1/b0 * x_i[i-1] - b2/b0 * x_i[i-2]

    b0, b1, b2 = flt.b_012
    a1, a2 = flt.a_12
    inv_b0 = inv(b0)

    ConvolutionFilter((inv_b0, inv_b0 * a1, inv_b0 * a2), (inv_b0 * b1, inv_b0 * b2))
end


function DSP.Biquad(flt::ConvolutionFilter{T}) where {T<:RealQuantity}
    DSP.Biquad(map(T, flt.b_012)..., map(T, flt.a_12)...)
end



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
    @inbounds @simd for i in eachindex(X)
        rh = fi.reverse_h

        @assert firstindex(y) == firstindex(x) == firstindex(rh) && lastindex(x) > lastindex(x)
        @assert lastindex(y) == lastindex(x) - (lastindex(rh) - firstindex(rh))

        for i in eachindex(y)
            y[i] = 0
            for j in eachindex(rh)
                y[i] = fma(rh[j], x[i+j-1], y[i])
            end
        end

        @assert axes(x) == firstindex(h):fir
        assert(firstindex(h) == firstindex(x))
    end
    Y
end



struct FFTConvFilterInstance{T<:Real,TV<:AbstractVector{Complex{T}}} <: AbstractConvFilterInstance{T}
    rfft_h::TV
end

_filterlen(fi::FFTConvFilterInstance) = size(fi.rfft_h, 1)


function rdfilt(fi::FFTConvFilterInstance, x::AbstractVector{T})
    y_ext = irfft(rfft(x) .* fi.rfft_h, size(x, 1))
    valid_range = (firstindex(y) + fi.n_filter - 1):lastindex(y)
    y_ext[valid_range]
end

function rdfilt!(y::AbstractVector{T}, fi::FFTConvFilterInstance, x::AbstractVector{T})
    y .= rdfilt(fi, x)
end


function bc_rdfilt(
    fi::FFTConvFilterInstance,
    inputs::ArrayOfSimilarArrays{<:RealQuantity,1}
) where {M,N}
    T_out = flt_output_smpltype(fi)
    X = flatview(inputs)
    Y = irfft(rfft(X, 1) .* fi.rfft_h, size(X, 1), 1)
    valid_range = (firstindex(y) + fi.n_filter - 1):lastindex(y)
    flat_output = Y[valid_range, :]
    ArrayOfSimilarArrays{T_out,M,N}(flat_output)
end

function bc_rdfilt!(
    outputs::ArrayOfSimilarArrays{<:RealQuantity,1},
    fi::FFTConvFilterInstance,
    inputs::ArrayOfSimilarArrays{<:RealQuantity,1}
)
    flatview(outputs) .= flatview(bc_rdfilt(fi, inputs))
    return outputs
end
