# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    shift_waveform(signal::AbstractSamples, a::RealQuantity)
    shift_waveform(wf::RDWaveform, a::RealQuantity)

Shifts each sample of a waveform up by `a`.
"""
function shift_waveform end
export shift_waveform

shift_waveform(signal::AbstractSamples, a::RealQuantity) = signal .+ a


"""
    multiply_waveform(signal::AbstractSamples, a::RealQuantity)
    multiply_waveform(wf::RDWaveform, a::RealQuantity)

Multiplies each sample of a waveform by `a`.
"""
function multiply_waveform end
export multiply_waveform

multiply_waveform(signal::AbstractSamples, a::RealQuantity) = signal .* a


"""
    reverse_waveform(signal::AbstractSamples)
    reverse_waveform(wf::RDWaveform)

Reverses the order of samples in a waveform.
"""
function reverse_waveform end
export reverse_waveform

reverse_waveform(signal::AbstractSamples) = reverse(signal)


const _WaveformOp = Union{typeof(shift_waveform), typeof(multiply_waveform), typeof(reverse_waveform)}

(f::_WaveformOp)(wf::RDWaveform, args::RealQuantity...) = RDWaveform(wf.time, f(wf.signal, args...))

function Base.Broadcast.broadcasted(f::_WaveformOp, bc_inputs, bc_args...)
    _bc_wfop(f, Base.materialize(bc_inputs), map(Base.materialize, bc_args)...)
end

_bc_wfop(f::_WaveformOp, inputs, args...) = broadcast((input, as...) -> f(input, as...), inputs, args...)

function _bc_wfop(f::_WaveformOp, inputs::ArrayOfRDWaveforms, args...)
    ArrayOfRDWaveforms((inputs.time, _bc_wfop(f, inputs.signal, args...)))
end

function _bc_wfop(f::_WaveformOp, inputs::ArrayOfSimilarVectors{<:RealQuantity}, args...)
    ArrayOfSimilarVectors(_bc_wfop_flat(f, flatview(inputs), args...))
end

_bc_wfop_flat(::typeof(shift_waveform), X::AbstractMatrix{<:RealQuantity}, a) = X .+ a'
_bc_wfop_flat(::typeof(multiply_waveform), X::AbstractMatrix{<:RealQuantity}, a) = X .* a'
_bc_wfop_flat(::typeof(reverse_waveform), X::AbstractMatrix{<:RealQuantity}) = reverse(X, dims = 1)
