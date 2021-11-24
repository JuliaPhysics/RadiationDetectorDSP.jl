# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    rc_filter(RC::Real)

Return a DSP.jl-compatible RC-filter.
"""
function rc_filter(RC::Real)
    T = float(typeof(RC))
    α = 1 / (1 + RC)
    Biquad(T(α), T(0), T(0), T(α - 1), T(0))
end
export rc_filter


"""
    cr_filter(RC::Real)

Return a DSP.jl-compatible CR-filter.
"""
function cr_filter(RC::Real)
    T = float(typeof(RC))
    α = RC / (RC + 1)
    Biquad(T(α), T(-α), T(0), T(-α), T(0))
end
export cr_filter


"""
    inv_cr_filter(RC::Real)

Return a DSP.jl-compatible inverse CR-filter.
"""
function inv_cr_filter(RC::Real)
    T = float(typeof(RC))
    α = 1 / (1 + RC)
    k = 1 + 1/RC
    Biquad(T(k), T(k * (α - 1)), T(0), T(-1), T(0))
end
export inv_cr_filter


"""
    crmod_filter(RC::Real)

Return a DSP.jl-compatible modified CR-filter.
"""
function crmod_filter(RC::Real)
    T = float(typeof(RC))
    α = RC / (RC + 1)
    Biquad(T(1), T(-1), T(0), T(-α), T(0))
end
export crmod_filter


"""
    inv_crmod_filter(RC::Real)

Return a DSP.jl-compatible inverse modified CR-filter.
"""
function inv_crmod_filter(RC::Real)
    T = float(typeof(RC))
    α = 1 / (1 + RC)
    Biquad(T(1), T(α - 1), T(0), T(-1), T(0))
end
export inv_crmod_filter


"""
    integrator_filter(gain::Real)

Return a DSP.jl-compatible integrator filter.
"""
function integrator_filter(gain::Real)
    T = float(typeof(gain))
    Biquad(T(gain), T(0), T(0), T(-1), T(0))
end
export integrator_filter


"""
    differentiator_filter(gain::Real)

Return a DSP.jl-compatible differentiator filter.
"""
function differentiator_filter(gain::Real)
    T = float(typeof(gain))
    Biquad(T(gain), T(-gain), T(0), T(0), T(0))
end
export differentiator_filter


"""
    integrator_cr_filter(gain::Real, RC::Real)  

Return a DSP.jl-compatible integrator plus CR filter.
"""
function integrator_cr_filter(gain::Real, RC::Real)
    T = float(promote_type(typeof(gain), typeof(RC)))
    α = 1 / (1 + RC)
    Biquad(T(gain), T(-α), T(0), T(α - 1), T(0))
end
export integrator_cr_filter


"""
    integrator_crmod_filter(gain::Real, RC::Real)

Return a DSP.jl-compatible integrator plus modified CR filter.
"""
function integrator_crmod_filter(gain::Real, RC::Real)
    T = float(promote_type(typeof(gain), typeof(RC)))
    α = 1 / (1 + RC)
    Biquad(T(gain), T(0), T(0), T(α - 1), T(0))
end
export integrator_crmod_filter


"""
simple_csa_response_filter(τ_rise::Real, τ_decay::Real, gain::Real = one(τ_rise))

Return a DSP.jl-compatible filter that models the response of a typical
charge-sensitive amplifier (CSA).
"""
function simple_csa_response_filter(τ_rise::Real, τ_decay::Real, gain::Real = one(τ_rise))
    # TODO: Use a single biquad filter

    T = float(promote_type(promote_type(typeof(τ_rise), typeof(τ_decay)), typeof(gain)))
    rc_filter(T(τ_rise)) * integrator_cr_filter(T(gain), T(τ_decay))
end
export simple_csa_response_filter
