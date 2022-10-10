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

When used as a type parameter value, marks linear behavior of a filter.
"""
abstract type DirectConvolution <: ConvolutionMethod end
export DirectConvolution


"""
    abstract type FFTConvolution <: ConvolutionMethod

When used as a type parameter value, marks non linear behavior of a filter.
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
    ConvolutionFilterInstance{T}(flt.b_012, flt.a_12, _smpllen(si))
end


# # true starting from i = non-zero-length of h:
# filt(FIRFilter(h), x) .≈ irfft(rfft(h) .* rfft(x), length(h))

# # Multi-waveform:
# irfft(rfft(H, 1) .* rfft(X, 1), size(H, 1), 1)


fcv.(Ref(x2), Ref(h2), 1:10) ≈ irfft(rfft(x2) .* rfft(h2), length(h2))

heatmap(deconv_img) 




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



struct ConvolutionFilterInstance{T} <: AbstractRadSigFilterInstance{DirectConvolution}
    b_012::NTuple{3,T}
    a_12::NTuple{2,T}
    n::Int
end

@inline function rdfilt!(Y::AbstractVector{T}, fi::ConvolutionFilterInstance{T}, X::AbstractVector{T}) where {T<:Real}
    a1, a2 = fi.a_12
    neg_a1, neg_a2 = -a1, -a2
    b0, b1, b2 = fi.b_012
    s1::T = zero(T) # s_init[1]
    s2::T = zero(T) # s_init[2]
    s3::T = zero(T)

    #!!! @assert eachindex(X) == eachindex(Y)

    @inbounds @simd for i in eachindex(X)
        x_i = T(X[i])

        z1 = fma(b0, x_i, s1)
        z2 = fma(b1, x_i, s2)
        z3 = fma(b2, x_i, s3)

        y_i = z1
        Y[i] = y_i

        s1 = fma(neg_a1, y_i, z2)
        s2 = fma(neg_a2, y_i, z3)
    end
    Y
end



flt_output_smpltype(fi::ConvolutionFilterInstance) = flt_input_smpltype(fi)
flt_input_smpltype(fi::ConvolutionFilterInstance{T}) where T = T
flt_output_length(fi::ConvolutionFilterInstance) = flt_input_length(fi) - length(eachindex(fi.coeffs))
flt_input_length(fi::ConvolutionFilterInstance) = fi.n_input
flt_output_time_axis(fi::ConvolutionFilterInstance, time::AbstractVector{<:RealQuantity}) = time
