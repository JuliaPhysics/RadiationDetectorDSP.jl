# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

module RadiationDetectorDSPReactantExt

using Reactant: Reactant, TracedRNumber

import RadiationDetectorDSP
using RadiationDetectorDSP: SamplingInfo, FFTConvFilterInstance, _conv_fft_length, _conv_spectrum, _convert_smpls, _filterlen

using FFTW: fft


const _TracedMatrix = AbstractMatrix{<:TracedRNumber}
const _ConcreteMatrix = Union{Reactant.ConcretePJRTArray{<:Number,2}, Reactant.ConcreteIFRTArray{<:Number,2}}

# The generic method transposes element by element:
RadiationDetectorDSP._nonlazy_transpose(X::_TracedMatrix) = permutedims(X)
RadiationDetectorDSP._nonlazy_transpose(X::_ConcreteMatrix) = permutedims(X)


# Filter instances are parameterized by the numerical type of the samples,
# not by the type of the symbolic numbers that stand in for them:
function RadiationDetectorDSP.elsmplinfo(smpls::AbstractArray{<:AbstractVector{TracedRNumber{T}}}) where {T}
    SamplingInfo{T}(RadiationDetectorDSP._inneraxes(smpls, 1))
end


RadiationDetectorDSP._convert_smpls(::Type{T}, X::AbstractArray{TracedRNumber{T}}) where {T} = X
RadiationDetectorDSP._convert_smpls(::Type{T}, X::AbstractArray{<:TracedRNumber}) where {T} = TracedRNumber{T}.(X)


# Without a fixed work-group size, kernels may be launched as a single block over the whole
# range, which fails silently for large batches of signals:
RadiationDetectorDSP._ka_threads(::RadiationDetectorDSP._KA_Backend, ::AbstractArray{<:TracedRNumber}) = (256,)


# XLA copies the input of inverse real FFTs and scales their output in a separate pass. So pairs
# of real signals are convolved as the real and imaginary part of one complex signal, using
# forward complex FFTs only: ifft(A) == conj(fft(conj(A))) / n_fft, and the spectrum of the
# taps includes the factor 1 / n_fft.
function RadiationDetectorDSP._fft_conv_ext(fi::FFTConvFilterInstance{T}, X::AbstractMatrix{<:TracedRNumber}) where T
    n_input, n_signals = size(X)
    n_fft = _conv_fft_length(fi, X)
    h = _conv_spectrum(fi, n_fft)
    # Conjugate of the full spectrum, from the half spectrum of the real taps:
    g = vcat(conj.(h), h[(end - iseven(n_fft)):-1:2])
    Xt = _convert_smpls(T, X)
    X_even = isodd(n_signals) ? hcat(Xt, zero(similar(Xt, n_input, 1))) : Xt
    n_pairs = size(X_even, 2) ÷ 2
    X_pairs = reshape(X_even, n_input, 2, n_pairs)
    Z = complex.(X_pairs[:, 1, :], X_pairs[:, 2, :])
    Z_ext = vcat(Z, zero(similar(Z, n_fft - n_input, n_pairs)))
    V = fft(conj.(fft(Z_ext, 1)) .* g, 1)[_filterlen(fi):n_input, :]
    n_out = size(V, 1)
    Y_pairs = cat(reshape(real.(V), n_out, 1, n_pairs), reshape(-imag.(V), n_out, 1, n_pairs); dims = 2)
    reshape(Y_pairs, n_out, 2 * n_pairs)[:, 1:n_signals]
end


end # module RadiationDetectorDSPReactantExt
