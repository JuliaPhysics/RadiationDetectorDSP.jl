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
