# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


function rc_filter(RC::Real)
    T = float(typeof(RC))
    α = 1 / (1 + RC)
    Biquad(T(α), T(0), T(0), T(α - 1), T(0))
end


function cr_filter(RC::Real)
    T = float(typeof(RC))
    α = RC / (RC + 1)
    Biquad(T(α), T(-α), T(0), T(-α), T(0))
end


function inv_cr_filter(RC::Real)
    T = float(typeof(RC))
    α = 1 / (1 + RC)
    k = 1 + 1/RC
    Biquad(T(k), T(k * (α - 1)), T(0), T(-1), T(0))
end


function crmod_filter(RC::Real)
    T = float(typeof(RC))
    α = RC / (RC + 1)
    Biquad(T(1), T(-1), T(0), T(-α), T(0))
end


function inv_crmod_filter(RC::Real)
    T = float(typeof(RC))
    α = 1 / (1 + RC)
    Biquad(T(1), T(α - 1), T(0), T(-1), T(0))
end


function integrator_filter(gain::Real)
    T = float(typeof(gain))
    Biquad(T(gain), T(0), T(0), T(-1), T(0))
end


function differentiator_filter(gain::Real)
    T = float(typeof(gain))
    Biquad(T(gain), T(-gain), T(0), T(0), T(0))
end


function integrator_cr_filter(gain::Real, RC::Real)
    T = float(promote_type(typeof(gain), typeof(RC)))
    α = 1 / (1 + RC)
    Biquad(T(gain), T(-α), T(0), T(α - 1), T(0))
end


function integrator_crmod_filter(gain::Real, RC::Real)
    T = float(promote_type(typeof(gain), typeof(RC)))
    α = 1 / (1 + RC)
    Biquad(T(gain), T(0), T(0), T(α - 1), T(0))
end


function simple_csa_response_filter(τ_rise::Real, τ_decay::Real, gain::Real = one(τ_rise))
    # TODO: Use a single biquad filter

    T = float(promote_type(promote_type(typeof(τ_rise), typeof(τ_decay)), typeof(gain)))
    rc_filter(T(τ_rise)) * integrator_cr_filter(T(gain), T(τ_decay))
end
