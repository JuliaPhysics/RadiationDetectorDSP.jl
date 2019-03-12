# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

function add_rect_pulse!(samples::AbstractSamples, start::Integer, pulselen::Integer, amplitude::Real = 1.0)
    T = eltype(samples)
    A = view(samples, start:(start + pulselen - 1))
    A .= A .+ T(amplitude)
    samples
end


gen_rect_pulse(tracelen::Integer, start::Integer, pulselen::Integer, amplitude::Real = 1.0) =
    add_rect_pulse!(fill(zero(amplitude), tracelen), start, pulselen, amplitude)


# add_gaussian_noise!(samples::AbstractSamples, sigma::AbstractFloat = 1.0) = ...
