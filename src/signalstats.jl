# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    signalstats(signal::AbstractSamples, start::Real, stop::Real)
    signalstats(signal::RDWaveform, start::RealQuantity, stop::RealQuantity)

Get statistics on `signal` in the interval (`start`,`stop`).
"""

function signalstats end
export signalstats


function signalstats(input::SamplesOrWaveform, start::RealQuantity, stop::RealQuantity)
    X_axis, Y = _get_axis_and_signal(input)
    # ToDo: Lower numerical precision of x-axis to y-axis, if x-axis is a range
    first_x, step_x = first(X_axis), step(X_axis)
    from = round(Int, ustrip(NoUnits, (start - first_x) / step_x)) + firstindex(X_axis)
    until = round(Int, ustrip(NoUnits, (stop - first_x) / step_x)) + firstindex(X_axis)
    _signalstats_impl(X_axis, Y, from:until)
end

function _signalstats_impl(X::AbstractArray{<:RealQuantity}, Y::AbstractArray{<:RealQuantity}, idxs::AbstractUnitRange{<:Integer})
    @assert axes(X) == axes(Y)
    @assert firstindex(X) <= first(idxs) <= last(idxs) <= lastindex(X)
    @assert firstindex(Y) <= first(idxs) <= last(idxs) <= lastindex(Y)

    zx = zero(eltype(X))
    zy = zero(eltype(Y))

    sum_X::float(typeof(zx)) = zx
    sum_Y::float(typeof(zy)) = zy
    sum_X_sqr::float(typeof(zx * zx)) = zx * zx
    sum_Y_sqr::float(typeof(zy * zy)) = zy * zy
    sum_XY::float(typeof(zx * zy)) = zx * zy

    @inbounds @fastmath @simd for i in idxs
        x, y = X[i], Y[i]
        sum_X = x + sum_X
        sum_X_sqr = fma(x, x, sum_X_sqr)
        sum_Y = y + sum_Y
        sum_Y_sqr = fma(y, y, sum_Y_sqr)
        sum_XY = fma(x, y, sum_XY)
    end

    n = length(idxs)
    inv_n = inv(n)

    # mean_X = sum_X * inv_n
    mean_Y = sum_Y * inv_n
    # var_X = sum_X_sqr * inv_n - mean_X * mean_X
    var_Y = max(sum_Y_sqr * inv_n - mean_Y * mean_Y, zero(mean_Y))
    # cov_XY = sum_XY * inv_n - mean_X * mean_Y

    # mean_X_uncert = sqrt( (sum_X_sqr - sum_X * mean_X) / (n - 1) )
    # mean_Y_uncert = sqrt( (sum_Y_sqr - sum_Y * mean_Y) / (n - 1) )

    # corr_XY = cov_XY / sqrt(var_X * var_Y)

    slope = cov_XY / var_X
    offset = mean_Y - slope * mean_X
    # slope_uncert = sqrt( (var_X*var_Y - cov_XY*cov_XY) / (n - 2) ) / var_X
    # offset_uncert = slope_uncert * sqrt(sum_X_sqr * inv_n)

    (
        mean = mean_Y,
        sigma = sqrt(var_Y),
        slope = slope,
        offset = offset,
    )
end
