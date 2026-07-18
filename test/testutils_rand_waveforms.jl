# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using Random
using ArraysOfArrays, FillArrays
using RadiationDetectorSignals


"""
    rand_detector_waveforms(
        rng::AbstractRNG, ::Type{T}, n_samples::Integer, n_waveforms::Integer;
        kwargs...
    )

Generate semi-realistic detector charge waveforms: rectangular current
pulses with random amplitudes, shaped by a CR-RC filter with random rise
and fall times, plus random baseline offsets and amplitude-dependent
Gaussian noise.

All times are in units of samples. Returns a named tuple with the field
`waveforms::ArrayOfRDWaveforms` and the true `amplitudes`, `tau_rises`,
`tau_falls`, `baseline_levels`, `noise_levels`, `i_trigger` and
`n_collect`.
"""
function rand_detector_waveforms(
    rng::AbstractRNG, ::Type{T}, n_samples::Integer, n_waveforms::Integer;
    amplitude_range::Tuple{Real,Real} = (100, 1000), amplitude_rand_map = (-) ∘ log,
    t_trigger::Real = 0.25 * n_samples, t_collect::Real = 0.025 * n_samples,
    tau_rise::Tuple{Real,Real} = n_samples .* (0.01, 0.05),
    tau_fall::Tuple{Real,Real} = n_samples .* (0.1, 0.2),
    baseline_sigma::Real = 0.025 * (amplitude_range[2] - amplitude_range[1]),
    noise_sigma::Real = 0.005 * (amplitude_range[2] - amplitude_range[1]),
    rel_fano_factor::Real = 4,
    time_axis::AbstractRange = Base.OneTo(n_samples),
) where {T<:Real}
    X = zeros(T, n_samples, n_waveforms)

    i_trigger = firstindex(X, 1) + round(Int, t_trigger)
    n_collect = round(Int, t_collect)

    min_amplitude = T(amplitude_range[1])
    amplitude_scale = T(amplitude_range[2] - amplitude_range[1])
    amplitudes = min_amplitude .+ amplitude_scale .* T.(amplitude_rand_map.(rand(rng, T, n_waveforms)))

    tau_rises = T(tau_rise[1]) .+ (T(tau_rise[2]) - T(tau_rise[1])) .* rand(rng, T, n_waveforms)
    tau_falls = T(tau_fall[1]) .+ (T(tau_fall[2]) - T(tau_fall[1])) .* rand(rng, T, n_waveforms)
    baseline_levels = T(baseline_sigma) .* randn(rng, T, n_waveforms)
    noise_levels = T(noise_sigma) .* (1 .+ sqrt.(T(rel_fano_factor) .* amplitudes ./ amplitude_scale))

    for j in axes(X, 2)
        current = amplitudes[j] / n_collect
        α_cr = tau_falls[j] / (tau_falls[j] + 1)
        β_rc = 1 / (tau_rises[j] + 1)
        icr = zero(T)
        rc = zero(T)
        for i in axes(X, 1)
            x = i_trigger <= i < i_trigger + n_collect ? current : zero(T)
            icr = α_cr * (icr + x)
            rc = β_rc * icr + (1 - β_rc) * rc
            X[i, j] = rc + baseline_levels[j] + noise_levels[j] * randn(rng, T)
        end
    end

    waveforms = ArrayOfRDWaveforms((Fill(time_axis, n_waveforms), ArrayOfSimilarVectors(X)))

    return (
        waveforms = waveforms,
        amplitudes = amplitudes,
        tau_rises = tau_rises,
        tau_falls = tau_falls,
        baseline_levels = baseline_levels,
        noise_levels = noise_levels,
        i_trigger = i_trigger,
        n_collect = n_collect,
    )
end
