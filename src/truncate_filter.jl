# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    struct TruncateFilter <: AbstractRadSigFilter{LinearFiltering}

Filter that truncates the input signal.

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
Base.@kwdef struct TruncateFilter{FI<: ClosedInterval{<:RealQuantity}} <: AbstractRadSigFilter{LinearFiltering}
    "interval to keep"
    interval::FI
end

export TruncateFilter

Adapt.adapt_structure(to, flt::TruncateFilter) = flt


function fltinstance(flt::TruncateFilter, si::SamplingInfo{T}) where T
    delta_t = step(si.axis)
    i0 = firstindex(si.axis)
    t0 = first(si.axis)
    from = i0 + round(Int, uconvert(NoUnits, (minimum(flt.interval)-t0) / delta_t))
    until = i0 + round(Int, uconvert(NoUnits, (maximum(flt.interval)-t0) / delta_t))
    idxs = from:until
    TruncateFilterInstance{T,typeof(idxs)}(idxs, _smpllen(si))
end



struct TruncateFilterInstance{T<:Real,R<:AbstractUnitRange{<:Integer}} <: AbstractRadSigFilterInstance{LinearFiltering}
    idxs::R
    n_input::Int
end


flt_output_smpltype(fi::TruncateFilterInstance) = flt_input_smpltype(fi)
flt_input_smpltype(fi::TruncateFilterInstance{T}) where T = T
flt_output_length(fi::TruncateFilterInstance) = length(fi.idxs)
flt_input_length(fi::TruncateFilterInstance) = fi.n_input

function flt_output_time_axis(fi::TruncateFilterInstance, time::AbstractVector{<:RealQuantity})
    time[fi.idxs]
end


function rdfilt!(y::AbstractVector{T}, fi::TruncateFilterInstance, x::AbstractVector{T}) where {T<:Real}
    copyto!(y, view(x, fi.idxs))
    return y
end


function bc_rdfilt!(
    outputs::ArrayOfSimilarArrays{<:Real,N},
    fi::TruncateFilterInstance,
    inputs::ArrayOfSimilarArrays{<:Real,N}
) where N
    Y = flatview(outputs)
    X = flatview(inputs)
    colons = ntuple(_ -> :, Val(N))
    copyto!(Y, view(X, fi.idxs, colons...))
    return outputs
end
