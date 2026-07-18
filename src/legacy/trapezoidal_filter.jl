# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


const RealOrSIMD{T<:Real} = Union{T,<:SIMD.Vec{N,<:T} where N}


"""
    charge_trapflt!(samples::AbstractVector{<:RealOrSIMD{<:AbstractFloat}}, navg::Integer, ngap::Integer)

Apply a trapezoidal FIR filter to a charge signal in `samples`.
"""
function charge_trapflt!(samples::AbstractVector{<:RealOrSIMD{<:AbstractFloat}}, navg::Integer, ngap::Integer)
    @fastmath begin
        T = eltype(samples)
        idxs = eachindex(samples)
        n = length(idxs)

        norm_factor = inv(T(navg))

        navg >= 1 && ngap >= 0 || throw(ArgumentError("Require navg >= 1 and ngap >= 0"))
        (2 * navg + ngap) <= length(idxs) || throw(ArgumentError("filter must not be longer than input"))

        flt_shift = 0

        offs1 = flt_shift
        offs2 = flt_shift + navg
        offs3 = flt_shift + navg + ngap
        offs4 = flt_shift + navg + ngap + navg

        padding_value = samples[last(idxs)]

        acc::T = zero(T)

        @inbounds @simd for j in idxs[offs3+1:offs4] acc += samples[j] end
        @inbounds @simd for j in idxs[offs1+1:offs2] acc -= samples[j] end

        # Read window edges before overwriting samples[i], keeping the
        # in-place single pass valid:
        @inbounds for i in idxs[1:n-offs4]
            v1 = samples[i + offs1]; v2 = samples[i + offs2]
            v3 = samples[i + offs3]; v4 = samples[i + offs4]
            samples[i] = acc * norm_factor
            acc = acc + v1 - v2 - v3 + v4
        end

        @inbounds for i in idxs[n-offs4+1:n-offs3]
            v1 = samples[i + offs1]; v2 = samples[i + offs2]; v3 = samples[i + offs3]
            samples[i] = acc * norm_factor
            acc = acc + v1 - v2 - v3 + padding_value
        end

        @inbounds for i in idxs[n-offs3+1:n-offs2]
            v1 = samples[i + offs1]; v2 = samples[i + offs2]
            samples[i] = acc * norm_factor
            acc = acc + v1 - v2
        end

        @inbounds for i in idxs[n-offs2+1:n-offs1]
            v1 = samples[i + offs1]
            samples[i] = acc * norm_factor
            acc = acc + v1 - padding_value
        end
    end

    samples
end
export charge_trapflt!
