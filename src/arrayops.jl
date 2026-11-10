# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

# Filter implementations in terms of whole-array operations on the flat
# matrix of a batch of signals (one signal per column). They need no kernel
# compiler and no element-wise access, so they run on any array type that
# supports broadcasting, slicing and `cumsum`, and can be traced by Reactant.


"""
    RadiationDetectorDSP.bc_rdfilt_arrayops(flt::AbstractRadSigFilter, inputs::ArrayOfSimilarVectors)
    RadiationDetectorDSP.bc_rdfilt_arrayops(fi::AbstractRadSigFilterInstance, inputs::ArrayOfSimilarVectors)

Like `flt.(inputs)`, but computed with whole-array operations instead of one
kernel work item per signal.
"""
function bc_rdfilt_arrayops end

bc_rdfilt_arrayops(flt::AbstractRadSigFilter, inputs) = bc_rdfilt_arrayops(fltinstance(flt, elsmplinfo(inputs)), inputs)

function bc_rdfilt_arrayops(fi::AbstractRadSigFilterInstance, inputs::ArrayOfSimilarVectors)
    ArrayOfSimilarVectors(_arrayops_rdfilt(fi, _col_major(flatview(inputs))))
end


_convert_smpls(::Type{T}, X::AbstractArray{T}) where {T} = X
_convert_smpls(::Type{T}, X::AbstractArray) where {T} = T.(X)

# Y with Y[i, :] == X[i - k, :], zeros in the first k rows:
function _shift_down(X::AbstractMatrix, k::Integer)
    n = size(X, 1)
    Y = fill!(similar(X), zero(eltype(X)))
    k < n && (view(Y, (k + 1):n, :) .= view(X, 1:(n - k), :))
    return Y
end

# P with P[j + 1, :] == sum(X[1:j, :], dims = 1), so P[1, :] is zero:
function _prefix_sums(X::AbstractMatrix)
    P = fill!(similar(X, size(X, 1) + 1, size(X, 2)), zero(eltype(X)))
    view(P, 2:size(P, 1), :) .= cumsum(X, dims = 1)
    return P
end

# Y[i, :] = p * Y[i - 1, :] + U[i, :], with Y[0, :] = 0, by recursive doubling:
# after the step with shift s, Y[i, :] holds the terms p^k * U[i - k, :] for k < 2s.
function _linear_recurrence(p::Real, U::AbstractMatrix)
    isone(p) && return cumsum(U, dims = 1)
    n = size(U, 1)
    Y, q, s = U, p, 1
    while s < n && !iszero(q)
        Y = Y .+ q .* _shift_down(Y, s)
        q, s = q * q, 2 * s
    end
    return Y
end


function _arrayops_rdfilt(fi::FirstOrderIIRInstance{T}, X::AbstractMatrix) where {T}
    b0, b1 = fi.b_01
    U = b0 .* X .+ b1 .* _shift_down(X, 1)
    _linear_recurrence(-fi.a_1[1], U)
end


# The state (y[i], y[i-1]) follows a two-dimensional linear recurrence with
# transition matrix A = [-a1 -a2; 1 0], solved by recursive doubling.
function _arrayops_rdfilt(fi::BiquadFilterInstance{T}, X::AbstractMatrix) where {T}
    b0, b1, b2 = fi.b_012
    a1, a2 = fi.a_12
    n = size(X, 1)
    Y1 = b0 .* X .+ b1 .* _shift_down(X, 1) .+ b2 .* _shift_down(X, 2)
    Y2 = fill!(similar(Y1), zero(eltype(Y1)))
    A, s = Float64[-a1 -a2; 1 0], 1
    while s < n && !iszero(A)
        S1, S2 = _shift_down(Y1, s), _shift_down(Y2, s)
        a11, a12, a21, a22 = T(A[1, 1]), T(A[1, 2]), T(A[2, 1]), T(A[2, 2])
        Y1, Y2 = Y1 .+ a11 .* S1 .+ a12 .* S2, Y2 .+ a21 .* S1 .+ a22 .* S2
        A, s = A * A, 2 * s
    end
    return Y1
end


function _arrayops_rdfilt(fi::TrapezoidalChargeFilterInstance, X::AbstractMatrix)
    T = flt_output_smpltype(fi)
    navg, navg2 = fi.navg, fi.navg2
    n, len = size(X, 1), _filterlen(fi)
    n_out = n - len + 1
    P = _prefix_sums(_convert_smpls(T, X))
    norm_factor, norm_factor2 = inv(T(navg)), inv(T(navg2))
    norm_factor .* (view(P, (len + 1):(n + 1), :) .- view(P, (len - navg + 1):(n - navg + 1), :)) .-
        norm_factor2 .* (view(P, (navg2 + 1):(n_out + navg2), :) .- view(P, 1:n_out, :))
end


function _arrayops_rdfilt(fi::DirectConvFilterInstance{T}, X::AbstractMatrix) where {T}
    rh = Array(fi.reverse_h)
    n_out = size(X, 1) - length(rh) + 1
    Y = rh[1] .* view(X, 1:n_out, :)
    for j in 2:length(rh)
        Y = Y .+ rh[j] .* view(X, j:(j + n_out - 1), :)
    end
    return Y
end


function _arrayops_rdfilt(fi::FFTConvFilterInstance{T}, X::AbstractMatrix) where {T}
    n_fft = _conv_fft_length(fi, X)
    rfft_h = adapt(get_compute_unit(X), _conv_spectrum(fi, n_fft))
    X_ext = fill!(similar(X, T, (n_fft, size(X, 2))), zero(T))
    view(X_ext, 1:fi.n_input, :) .= X
    Y_ext = brfft(rfft(X_ext, 1) .* rfft_h, n_fft, 1)
    Y_ext[_filterlen(fi):fi.n_input, :]
end
