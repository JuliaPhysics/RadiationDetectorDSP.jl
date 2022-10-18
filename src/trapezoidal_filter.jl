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
Base.@kwdef struct TrapezoidalChargeFilter{
    T <: RealQuantity
} <: AbstractRadFIRFilter
    "averaging time"
    avgtime::T

    "gap time"
    gaptime::T
end

export TrapezoidalChargeFilter

Adapt.adapt_structure(to, flt::TrapezoidalChargeFilter) = flt

function ConvolutionFilter(flt::TrapezoidalChargeFilter)
    navg = Int(flt.avgtime)
    ngap = Int(flt.gaptime)
    norm_factor = inv(navg)
    T = typeof(norm_factor)
    params = zeros(T, 2*navg + ngap)
    fill!(view(params, firstindex(params):(firstindex(params)+navg-1)), +norm_factor)
    fill!(view(params, firstindex(params)+navg+ngap:lastindex(params)), -norm_factor)
    ConvolutionFilter(FFTConvolution(), params)
end


function fltinstance(flt::TrapezoidalChargeFilter, si::SamplingInfo{T}) where T
    delta_t = step(si.axis)

    navg = round(Int, uconvert(NoUnits, flt.avgtime / delta_t))
    ngap = round(Int, uconvert(NoUnits, flt.gaptime / delta_t))

    navg >= 1 && ngap >= 0 || throw(ArgumentError("Require navg >= 1 and ngap >= 0"))
    (2 * navg + ngap) <= _smpllen(si) || throw(ArgumentError("filter must not be longer than input"))

    TrapezoidalChargeFilterInstance{T}(navg, ngap, _smpllen(si))
end


struct TrapezoidalChargeFilterInstance{T<:RealQuantity} <: AbstractRadSigFilterInstance{LinearFiltering}
    navg::Int
    ngap::Int
    n_input::Int
end


_filterlen(fi::TrapezoidalChargeFilterInstance) = 2* fi.navg + fi.ngap

function flt_output_time_axis(fi::TrapezoidalChargeFilterInstance, time::AbstractVector{<:RealQuantity})
    valid_range = (firstindex(time) + _filterlen(fi) - 1):lastindex(time)
    time[valid_range]
end


@inline function rdfilt!(y::AbstractVector{T}, fi::TrapezoidalChargeFilterInstance{T}, x::AbstractVector{T}) where {T<:Real}
    @fastmath begin
        navg = fi.navg
        ngap = fi.ngap

        @assert firstindex(y) == firstindex(x)
        @assert lastindex(y) == lastindex(x) - _filterlen(fi) + 1

        norm_factor = inv(T(navg))
        acc::T = zero(T)

        offs1 = 0
        offs2 = offs1 + navg
        offs3 = offs2 + ngap
        offs4 = offs3 + navg

        @assert lastindex(y) + offs4 - 1 == lastindex(x)

        @inbounds @simd for i in firstindex(x):(firstindex(x) + navg - 1)
            #@info "YYYY" i + offs1 i + offs3
            acc = acc - x[i + offs1] + x[i + offs3]
        end
        y[firstindex(y)] = acc * norm_factor
 
        @inbounds @simd for i in firstindex(x):(lastindex(x) - offs4)
            acc = acc + x[i + offs1] - x[i + offs2] - x[i + offs3] + x[i + offs4]
            y[i + 1] = acc * norm_factor
        end
    end

    return y
end

function bc_rdfilt!(
    outputs::ArrayOfSimilarVectors{<:RealQuantity},
    fi::TrapezoidalChargeFilterInstance,
    inputs::ArrayOfSimilarVectors{<:RealQuantity}
)
    _ka_bc_rdfilt!(outputs, fi, inputs)
end



flt_output_smpltype(fi::TrapezoidalChargeFilterInstance) = flt_input_smpltype(fi)
flt_input_smpltype(fi::TrapezoidalChargeFilterInstance{T}) where T = T
flt_output_length(fi::TrapezoidalChargeFilterInstance) = flt_input_length(fi) - _filterlen(fi) + 1
flt_input_length(fi::TrapezoidalChargeFilterInstance) = fi.n_input
