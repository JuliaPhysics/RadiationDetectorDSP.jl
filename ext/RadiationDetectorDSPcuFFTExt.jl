# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

module RadiationDetectorDSPcuFFTExt

using CUDACore: CUDACore, CuContext, CuMatrix, DenseCuMatrix, HandleCache, context, stream
using cuFFT: cuFFT, cufftHandle

import RadiationDetectorDSP
using RadiationDetectorDSP: FFTConvFilterInstance, _conv_fft_length, _conv_spectrum, _filterlen

using Adapt: adapt
using HeterogeneousComputing: get_compute_unit


# Plans for batches of real-to-complex (forward) or complex-to-real transforms of length n_fft,
# with real signals rdist samples apart and complex spectra stored contiguously. Key: context,
# sample type, forward, n_fft, rdist, batch size.
const _PlanKey = Tuple{CuContext,DataType,Bool,Int,Int,Int}

function _make_plan((ctx, T, forward, n_fft, rdist, batch)::_PlanKey)
    CUDACore.context!(ctx) do
        handle_ref = Ref{cufftHandle}()
        cuFFT.cufftCreate(handle_ref)
        handle = handle_ref[]
        rtype, ctype = convert(CUDACore.cudaDataType, T), convert(CUDACore.cudaDataType, Complex{T})
        cdist = n_fft ÷ 2 + 1
        itype, idist, otype, odist = forward ? (rtype, rdist, ctype, cdist) : (ctype, cdist, rtype, rdist)
        n = Clonglong[n_fft]
        cuFFT.cufftXtMakePlanMany(handle, 1, n, n, 1, idist, itype, n, 1, odist, otype, batch, Ref{Csize_t}(), ctype)
        handle
    end
end

_destroy_plan(key::_PlanKey, handle::cufftHandle) = CUDACore.context!(() -> cuFFT.cufftDestroy(handle), key[1])

const _plans = HandleCache{_PlanKey,cufftHandle}(_make_plan, _destroy_plan)

function _exec_fft!(output::DenseCuMatrix, input::DenseCuMatrix, ::Type{T}, forward::Bool, n_fft::Integer,
                    rdist::Integer, batch::Integer, offset_in::Integer = 0, offset_out::Integer = 0) where {T}
    key = (context(), T, forward, Int(n_fft), Int(rdist), Int(batch))
    handle = pop!(_plans, key)
    try
        cuFFT.cufftSetStream(handle, stream())
        direction = forward ? cuFFT.CUFFT_FORWARD : cuFFT.CUFFT_INVERSE
        cuFFT.cufftXtExec(handle, pointer(input, offset_in + 1), pointer(output, offset_out + 1), direction)
    finally
        push!(_plans, key, handle)
    end
    return output
end


# The kept output samples n_filter:n_input of the circular convolution only depend on the first
# n_input samples of each padded signal. So the forward FFT can read the columns of X in place,
# with transform length n_fft and n_input samples between signals: the padding then holds samples
# of the next signal. Only the last signals, whose transforms would read past the end of X, are
# padded explicitly (all signals, if X is not aligned as cuFFT requires). A non-finite sample may
# therefore also spoil the output of the previous signal. The complex-to-real transform overwrites its input, which is a temporary array here.
function RadiationDetectorDSP._fft_conv_ext(fi::FFTConvFilterInstance{T}, X::DenseCuMatrix{T}) where {T<:Union{Float32,Float64}}
    n_input, n_signals = size(X)
    n_fft = _conv_fft_length(fi, X)
    n_freq = n_fft ÷ 2 + 1
    rfft_h = adapt(get_compute_unit(X), _conv_spectrum(fi, n_fft))
    F = CuMatrix{Complex{T}}(undef, n_freq, n_signals)
    # cuFFT requires real input to be aligned like complex numbers:
    aligned = UInt(pointer(X)) % sizeof(Complex{T}) == 0
    n_direct = aligned && n_signals * n_input >= n_fft ? (n_signals * n_input - n_fft) ÷ n_input + 1 : 0
    n_direct > 0 && _exec_fft!(F, X, T, true, n_fft, n_input, n_direct)
    if n_direct < n_signals
        n_tail = n_signals - n_direct
        X_tail = CuMatrix{T}(undef, n_fft, n_tail)
        fill!(view(X_tail, (n_input + 1):n_fft, :), zero(T))
        copyto!(view(X_tail, 1:n_input, :), view(X, :, (n_direct + 1):n_signals))
        _exec_fft!(F, X_tail, T, true, n_fft, n_fft, n_tail, 0, n_direct * n_freq)
    end
    F .*= rfft_h
    # In place, the real signals are 2 * n_freq samples apart:
    Y_ext = reshape(reinterpret(T, F), 2 * n_freq, n_signals)
    _exec_fft!(Y_ext, F, T, false, n_fft, 2 * n_freq, n_signals)
    view(Y_ext, _filterlen(fi):n_input, :)
end


end # module RadiationDetectorDSPcuFFTExt
