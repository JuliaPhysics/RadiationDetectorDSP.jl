# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    rand_waveforms(
        rng::AbstractRNG, ::Type{T}, n_samples::Integer, n_waveforms::Integer;
        kwargs...
    )

Generate semi-realistic detector charge waveforms: rectangular current
pulses with random amplitudes, shaped by an integrating CR-RC filter with
random rise and fall times, plus random baseline levels and
amplitude-dependent Gaussian noise.

All times are in units of samples. The waveforms are generated on the
device of `rng`, as the columns of a matrix.

Returns a named tuple with the fields `waveforms`, `amplitudes`,
`tau_rises`, `tau_falls`, `baseline_levels`, `noise_levels`, `i_trigger`
and `n_collect`.
"""
function rand_waveforms(
    rng::AbstractRNG, ::Type{T}, n_samples::Integer, n_waveforms::Integer;
    amplitude_range::Tuple{Real,Real} = (100,1000), amplitude_rand_map::Function = (-) ∘ log,
    t_trigger::Real = 0.25 * n_samples, t_collect::Real = 0.025 * n_samples,
    tau_rise::Tuple{Real,Real} = n_samples .* (0.01, 0.05),
    tau_fall::Tuple{Real,Real} = n_samples .* (0.1, 0.2),
    baseline_sigma::Real = 0.025 * (amplitude_range[2] - amplitude_range[1]),
    noise_sigma::Real = 0.005 * (amplitude_range[2] - amplitude_range[1]),
    rel_fano_factor::Real = 4,
) where {T<:Real}
    # zero-values samples:
    X = rand(rng, T, n_samples, n_waveforms)
    fill!(X, zero(T))

    # Fill with rectangular pulses:
    i_trigger = first(axes(X, 1)) .+ round(Int, t_trigger)
    n_collect = round(Int, t_collect)
    view(X, i_trigger:(i_trigger+n_collect-1), :) .+= one(T)

    min_amplitude = T(amplitude_range[1])
    amplitude_scale = T(amplitude_range[2] - amplitude_range[1])

    # Scale with random amplitudes
    amplitudes = min_amplitude .+ amplitude_scale * amplitude_rand_map.(rand(rng, T, n_waveforms))
    # divide by n_collect due to integration in filter below:
    X .*= amplitudes' ./ T(n_collect)

    tau_rises = T(tau_rise[1]) .+ (T(tau_rise[2]) .- T(tau_rise[1])) .* rand(rng, T, n_waveforms)
    tau_falls = T(tau_fall[1]) .+ (T(tau_fall[2]) .- T(tau_fall[1])) .* rand(rng, T, n_waveforms)

    # Apply integrating CR-RC filter for rise and fall time:
    α_cr = tau_falls ./ (tau_falls .+ 1)
    β_rc = 1 ./ (tau_rises .+ 1)

    x_curr  = similar(X, n_waveforms)
    icr_curr = similar(X, n_waveforms)
    rc_curr = similar(X, n_waveforms)
    x_prev  = similar(X, n_waveforms)
    icr_prev = similar(X, n_waveforms)
    rc_prev = similar(X, n_waveforms)

    idxs_i = axes(X, 1)
    x_curr .= view(X, first(idxs_i), :)
    x_prev .= x_curr
    icr_prev .= x_prev
    rc_prev .= x_prev

    X_rows = eachrow(X)

    # ToDo: Doesn't trace with Reactant yet.
    @trace for i in eachindex(X_rows)[begin+1:end]
        x_curr .= X_rows[i]
        icr_curr .= α_cr .* (icr_prev .+ x_curr)
        rc_curr .= β_rc .* icr_curr .+ (1 .- β_rc) .* rc_prev
        X[i, :] .= rc_curr
        x_prev  .= x_curr
        icr_prev .= icr_curr
        rc_prev .= rc_curr
    end

    # Add random baseline levels:
    baseline_levels = T(baseline_sigma) .* randn(rng, T, n_waveforms)
    X .+= baseline_levels'

    # Add random Gaussian noise:
    noise_levels = T(noise_sigma) .* (1 .+ sqrt.(T(rel_fano_factor) .* amplitudes ./ amplitude_scale))

    noise = randn(rng, T, n_samples, n_waveforms) .* noise_levels'
    X .+= noise

    return (
        waveforms = X,
        amplitudes = amplitudes,
        tau_rises = tau_rises,
        tau_falls = tau_falls,
        baseline_levels = baseline_levels,
        noise_levels = noise_levels,
        i_trigger = i_trigger,
        n_collect = n_collect,
    )
end
