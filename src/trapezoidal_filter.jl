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
    "pre-rise averaging time"
    avgtime::T

    "gap time"
    gaptime::T

    "post-rise averaging time"
    avgtime2::T
    # alternative constructor for symmetric trapezodial filters
    function TrapezoidalChargeFilter(avgtime::T, gaptime::T) where {T <: RealQuantity} 
        new{T}(avgtime, gaptime, avgtime)
    end
    function TrapezoidalChargeFilter(avgtime::T, gaptime::T, avgtime2::T) where {T <: RealQuantity} 
        new{T}(avgtime, gaptime, avgtime2)
    end
end

export TrapezoidalChargeFilter

Adapt.adapt_structure(to, flt::TrapezoidalChargeFilter) = flt

function ConvolutionFilter(flt::TrapezoidalChargeFilter)
    navg  = Int(flt.avgtime)
    navg2 = Int(flt.avgtime2)
    ngap  = Int(flt.gaptime)
    norm_factor = inv(navg)
    T = typeof(norm_factor)
    params = zeros(T, navg + ngap + navg2)
    fill!(view(params, firstindex(params):(firstindex(params)+navg-1)), +norm_factor)
    fill!(view(params, firstindex(params)+navg+ngap:lastindex(params)), -norm_factor)
    ConvolutionFilter(FFTConvolution(), params)
end


function fltinstance(flt::TrapezoidalChargeFilter, si::SamplingInfo{T}) where {T <: RealQuantity} 
    delta_t = step(si.axis)

    navg2  = round(Int, uconvert(NoUnits, flt.avgtime  / delta_t))
    ngap  = round(Int, uconvert(NoUnits, flt.gaptime  / delta_t))
    navg = round(Int, uconvert(NoUnits, flt.avgtime2 / delta_t))

    navg >= 1 && ngap >= 0 && navg2 >= 1 || throw(ArgumentError("Require navg/navg2 >= 1 and ngap >= 0"))
    (navg + ngap + navg2) <= _smpllen(si) || throw(ArgumentError("filter must not be longer than input"))

    TrapezoidalChargeFilterInstance{T}(navg, ngap, navg2, _smpllen(si))
end


struct TrapezoidalChargeFilterInstance{T<:RealQuantity} <: AbstractRadSigFilterInstance{LinearFiltering}
    navg::Int
    ngap::Int
    navg2::Int
    n_input::Int
end


_filterlen(fi::TrapezoidalChargeFilterInstance) = fi.navg + fi.ngap + fi.navg2

function flt_output_time_axis(fi::TrapezoidalChargeFilterInstance, time::AbstractVector{<:RealQuantity})
    valid_range = (firstindex(time) + _filterlen(fi) - 1):lastindex(time)
    time[valid_range]
end


@inline function rdfilt!(y::AbstractVector{T}, fi::TrapezoidalChargeFilterInstance{U}, x::AbstractVector{U}) where {T<:Real, U<:Real}
    @fastmath begin
        navg  = fi.navg
        navg2 = fi.navg2
        ngap  = fi.ngap

        @assert firstindex(y) == firstindex(x)
        @assert lastindex(y) == lastindex(x) - _filterlen(fi) + 1

        norm_factor = inv(T(navg))
        # println(norm_factor)
        norm_factor2 = inv(T(navg2))
        acc::T = zero(T)

        offs1 = 0
        offs2 = offs1 + navg
        offs3 = offs2 + ngap
        offs4 = offs3 + navg2

        @assert lastindex(y) + offs4 - 1 == lastindex(x)

        @inbounds @simd for i in firstindex(x):(firstindex(x) + navg - 1)
            #@info "YYYY" i + offs1 i + offs3
            # acc = acc + x[i + offs1] - x[i + offs3]
            acc = acc + x[i + offs1] * norm_factor - x[i + offs3] * norm_factor
        end
        # y[firstindex(y)] = acc * norm_factor
        y[firstindex(y)] = acc
        @inbounds @simd for i in firstindex(x):(lastindex(x) - offs4)
            # acc = acc + (x[i + offs1] - x[i + offs2]) * norm_factor - (x[i + offs3] - x[i + offs4]) * norm_factor2
            acc = acc + (x[i + offs1] - x[i + offs2]) * norm_factor - (x[i + offs3] - x[i + offs4]) * norm_factor2
            # acc = acc + x[i + offs1] - x[i + offs2] - x[i + offs3] + x[i + offs4]
            # y[i + 1] = acc * norm_factor2
            y[i + 1] = acc
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



flt_output_smpltype(fi::TrapezoidalChargeFilterInstance) = float(flt_input_smpltype(fi))
flt_input_smpltype(fi::TrapezoidalChargeFilterInstance{T}) where T = T
flt_output_length(fi::TrapezoidalChargeFilterInstance) = flt_input_length(fi) - _filterlen(fi) + 1
flt_input_length(fi::TrapezoidalChargeFilterInstance) = fi.n_input
