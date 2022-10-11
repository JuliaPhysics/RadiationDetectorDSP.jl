# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    struct TrapezoidalChargeFilter <: AbstractRadNonlinearFilter

Filter that responds to a step signal with a trapezoidal pulse.

The filter is equivalent to two moving averages separated by a gap.

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)

A sharp step on the input will result in a trapezoid with rise time and fall
time `avgtime` and a flat top of length `gaptime`.
"""
@with_kw struct TrapezoidalChargeFilter{
    T <: RealQuantity
} <: AbstractRadIIRFilter
    "averaging time"
    avgtime::T,

    "gap time"
    gaptime::T
end

export TrapezoidalChargeFilter

#Adapt.adapt_structure(to, flt::TrapezoidalChargeFilter) = flt


function fltinstance(flt::TrapezoidalChargeFilter, si::SamplingInfo{T}) where T
    navg = Int(flt.avgtime)
    ngap = Int(flt.gaptime)

    navg >= 1 && ngap >= 0 || throw(ArgumentError("Require navg >= 1 and ngap >= 0"))
    (2 * navg + ngap) <= length(idxs_x) || throw(ArgumentError("filter must not be longer than input"))

    TrapezoidalChargeFilterInstance(navg, ngap,_smpllen(si))
end


struct TrapezoidalChargeFilterInstance <: AbstractRadSigFilterInstance{LinearFiltering}
    navg::Int,
    ngap::Int
    n_input::Int
end


_filterlen(fi::TrapezoidalChargeFilterInstance) = 2* fi.navg + f.ngap

function flt_output_time_axis(fi::TrapezoidalChargeFilterInstance, time::AbstractVector{<:RealQuantity})
    valid_range = (firstindex(time) + _filterlen(fi) - 1):lastindex(time)
    time[valid_range]
end


@inline function rdfilt!(y::AbstractVector{T}, fi::AbstractConvFilterInstance{T}, x::AbstractVector{T}) where {T<:Real}
    @fastmath begin
        @assert firstindex(y) == firstindex(x) == firstindex(rh) && lastindex(x) >= lastindex(rh)
        @assert lastindex(y) == lastindex(x) - (lastindex(rh) - firstindex(rh))

        T = eltype(x)
        norm_factor = inv(T(navg))
        acc::T = zero(T)

        offs1 = 0
        offs2 = offs1 + navg
        offs3 = offs2 + ngap
        offs4 = offs3 + navg

        @assert lastindex(y) + offs4 - 1 == lastindex(x)

        #@inbounds @simd
        for i in firstindex(x):(firstindex(x) + navg - 1)
            acc = acc - x[i + offs1] + x[i + offs3]
        end
        y[firstindex(y)] = acc
 
        #@inbounds @simd
        for i in firstindex(x):(lastindex(x) - offs4)
            acc = acc + x[i + offs1] - x[i + offs2] - x[i + offs3] + x[i + offs4]
            y[i + 1] = acc * norm_factor
        end
    end

    return y
end


flt_output_smpltype(fi::TrapezoidalChargeFilterInstance) = flt_input_smpltype(fi)
flt_input_smpltype(fi::TrapezoidalChargeFilterInstance{T}) where T = T
flt_output_length(fi::TrapezoidalChargeFilterInstance) = flt_input_length(fi) - _filterlen(fi) + 1
flt_input_length(fi::TrapezoidalChargeFilterInstance) = fi.n_input
