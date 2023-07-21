# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    shift_waveform(signal::AbstractSamples, a::RealQuantity)
    shift_waveform(wf::RDWaveform, a::RealQuantity)

Shifts each sample of a waveform up by `a`.
"""
function shift_waveform end
export shift_waveform

shift_waveform(signal::AbstractSamples, a::RealQuantity) = signal .+ a

shift_waveform(wf::RDWaveform, a::RealQuantity) = RDWaveform(wf.time, shift_waveform(wf.signal, a))

function Base.Broadcast.broadcasted(::typeof(shift_waveform), inputs, a)
    bc_shift_waveform(Base.materialize(inputs), Base.materialize(a))
end

_nobcspec_shift_waveform(inputs, a) = shift_waveform(inputs, a)
bc_shift_waveform(inputs, a) = _nobcspec_shift_waveform.(inputs, a)

function bc_shift_waveform(inputs::ArrayOfRDWaveforms, a)
    ArrayOfRDWaveforms((inputs.time, broadcast(shift_waveform, inputs.signal, a)))
end

function bc_shift_waveform(inputs::ArrayOfSimilarVectors{<:RealQuantity}, a)
    X = flatview(inputs)
    Y = X .+ a'
    ArrayOfSimilarVectors(Y)
end


"""
    multiply_waveform(signal::AbstractSamples, a::RealQuantity)
    multiply_waveform(wf::RDWaveform, a::RealQuantity)
    
Multiplies each sample of a waveform by `a`.
"""
function multiply_waveform end
export multiply_waveform

multiply_waveform(signal::AbstractSamples, a::RealQuantity) = signal .* a

multiply_waveform(wf::RDWaveform, a::RealQuantity) = RDWaveform(wf.time, multiply_waveform(wf.signal, a))

function Base.Broadcast.broadcasted(::typeof(multiply_waveform), inputs, a)
    bc_multiply_waveform(Base.materialize(inputs), Base.materialize(a))
end

_nobcspec_multiply_waveform(inputs, a) = multiply_waveform(inputs, a)
bc_multiply_waveform(inputs, a) = _nobcspec_multiply_waveform.(inputs, a)

function bc_multiply_waveform(inputs::ArrayOfRDWaveforms, a)
    ArrayOfRDWaveforms((inputs.time, broadcast(multiply_waveform, inputs.signal, a)))
end

function bc_multiply_waveform(inputs::ArrayOfSimilarVectors{<:RealQuantity}, a)
    X = flatview(inputs)
    Y = X .* a'
    ArrayOfSimilarVectors(Y)
end


"""
    reverse_waveform(signal::AbstractSamples)
    reverse_waveform(wf::RDWaveform)

Reverses the order of samples in a waveform.
"""
function reverse_waveform end
export reverse_waveform

reverse_waveform(signal::AbstractSamples) = reverse(signal)

reverse_waveform(wf::RDWaveform) = RDWaveform(wf.time, reverse_waveform(wf.signal))

function Base.Broadcast.broadcasted(::typeof(reverse_waveform), inputs)
    bc_reverse_waveform(Base.materialize(inputs))
end

_nobcspec_reverse_waveform(inputs) = reverse_waveform(inputs)
bc_reverse_waveform(inputs) = _nobcspec_reverse_waveform.(inputs)

function bc_reverse_waveform(inputs::ArrayOfRDWaveforms)
    ArrayOfRDWaveforms((inputs.time, broadcast(reverse_waveform, inputs.signal)))
end

function bc_reverse_waveform(inputs::ArrayOfSimilarVectors{<:RealQuantity})
    X = flatview(inputs)
    Y = reverse(X, dims = 1)
    ArrayOfSimilarVectors(Y)
end
