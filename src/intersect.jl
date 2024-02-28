# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    struct Intersect <: Function

Finds the intersects of a Y with a threshold

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct Intersect{T<:RealQuantity} <: Function
    "minimum time-over-threshold"
    mintot::T = 4
end

export Intersect

function (f::Intersect)(input::SamplesOrWaveform, threshold::RealQuantity)
    X_axis, Y = _get_axis_and_signal(input)
    min_n_over_thresh = max(1, round(Int, ustrip(NoUnits, f.mintot / step(X_axis))))
    _find_intersect_impl(X_axis, Y, threshold, min_n_over_thresh)
end


function _find_intersect_impl(X::AbstractVector{<:RealQuantity}, Y::AbstractVector{<:RealQuantity}, threshold::RealQuantity, min_n_over_thresh::Int)
    @assert axes(X) == axes(Y)

    # ToDo: What if eltype(Y) is based on ForwardDiff.Dual, but eltype(X) is not?
    R = float(eltype(X))

    if isempty(Y)
        return (
            x = R(NaN),
            multiplicity = -1
        )
    end

    # GPU-friendly branch-free code:

    cand_pos::Int = firstindex(Y) + 1
    intersect_pos::Int = firstindex(Y) + 1
    y_high_counter::Int = ifelse(first(Y) >= threshold, min_n_over_thresh + 1, 0)

    n_intersects::Int = 0
    @inbounds for i in eachindex(Y)
        y = Y[i]
        y_is_high = y >= threshold
        first_high_y = y_high_counter == 0
        cand_pos = ifelse(y_is_high && first_high_y, i, cand_pos)
        y_high_counter = ifelse(y_is_high, y_high_counter + 1, 0)
        new_intersect_found = y_high_counter == min_n_over_thresh
        n_intersects = ifelse(new_intersect_found, n_intersects + 1, n_intersects)
        no_previous_intersecs = n_intersects == 1
        first_intersect_found = new_intersect_found && no_previous_intersecs
        intersect_pos = ifelse(first_intersect_found, cand_pos, intersect_pos)
    end

    #TODO: return NaN if position found is unphysical but make sure it is compatible with the other routines
    # @assert intersect_pos > firstindex(Y)
    if intersect_pos <= firstindex(Y)
        intersect_pos = firstindex(Y) + 1
    end

    # Linear interpolation:
    x_l = X[intersect_pos - 1]
    x_r = X[intersect_pos]
    y_l = Y[intersect_pos - 1]
    y_r = Y[intersect_pos]
    intersect_cand_x = R(threshold - y_l) * R(x_r - x_l) / R(y_r - y_l) + R(x_l)
    @assert  y_l <= threshold <= y_r || n_intersects == 0
    # TODO: return NaN if no intersect found but make sure it is compatible with the other routines
    intersect_x  = ifelse(n_intersects > 0, intersect_cand_x, zero(intersect_cand_x))
    n_intersects = ifelse(n_intersects > 0, n_intersects, 0)

    return (
        x = intersect_x,
        multiplicity = n_intersects
    )
end
