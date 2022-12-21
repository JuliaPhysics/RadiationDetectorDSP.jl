# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


_vandermonde(X::AbstractVector{<:Real}, degree::Integer) = X .^ (0:degree)'

_onehot_vec(n::Integer, i::Integer, x::Number) = [ifelse(j == i, x, zero(x)) for j in 1:n]


function _lsq_fit_matrix(X, degree::Integer)
    V = _vandermonde(X, degree)
    V * inv(V'V)' # return transpose of inv(V'V) * V' for advantagous memory layout
end


function _lsqfitatpos_single_impl(A_slqfit::AbstractMatrix{<:Real}, Y::AbstractVector{<:Real}, x::Real)
    # Equivalent to dot(A_slqfit * Y, x.^(0:(size(A,1)-1)))

    @assert axes(A_slqfit, 1) == axes(Y, 1)
    R = promote_type(eltype(A_slqfit), eltype(Y), typeof(x))
    y::R = 0
    x_jm1::R = 1  # cache x^(j-1)
    @inbounds for j in axes(A_slqfit, 2)
        coeff_jm1::R = 0  # polynomial coefficient c_{j-1}
        @simd for i in axes(A_slqfit, 1)
            coeff_jm1 = muladd( A_slqfit[i,j], Y[i], coeff_jm1)
        end
        y += coeff_jm1 * x_jm1
        x_jm1 *= x
    end
    return y
end

function _lsqfitatpos_impl(A_slqfit::AbstractMatrix{<:Real}, X::AbstractRange{<:RealQuantity}, Y::AbstractVector{<:Real}, x::RealQuantity)
    n = size(A_slqfit, 1)
    n_m1_d2 = div(n-1, 2)
    first_x, step_x = first(X), step(X)
    cand_xi = ustrip(NoUnits, (x - first_x) / step_x) + firstindex(X)
    cand_from = floor(Integer, cand_xi) - n_m1_d2
    cand_until = cand_from + n - 1
    in_range = cand_from >= firstindex(Y) && cand_until + 1 <= lastindex(Y)
    dummy_from = oftype(cand_from, firstindex(Y))
    dummy_until = dummy_from + n - 1
    dummy_xi = oftype(cand_xi, dummy_from + n_m1_d2)
    xi, from, until = ifelse(
        in_range,
        (cand_xi, cand_from, cand_until),
        (dummy_xi, dummy_from, dummy_until)
    )

    range_r = (from+1):(until+1)
    shifted_x_r = xi - first(range_r)
    #weight_r = _smoothstep(mod(shifted_x_r, 1))
    weight_r = mod(shifted_x_r, 1)
    y_r = _lsqfitatpos_single_impl(A_slqfit, view(Y, range_r), shifted_x_r)

    range_l = (from):(until)
    shifted_x_l = xi - first(range_l)
    weight_l = 1 - weight_r
    y_l = _lsqfitatpos_single_impl(A_slqfit, view(Y, range_l), shifted_x_l)
  
    result = weight_r * y_r + weight_l * y_l
    ifelse(in_range, result, oftype(result, NaN))
end


function lsqfitatpos(degree::Integer, n::Integer, X::AbstractRange{<:RealQuantity}, Y::AbstractVector{<:Real}, x::RealQuantity)
    A_slqfit = _lsq_fit_matrix(0:(n - 1), degree)
    _lsqfitatpos_impl(A_slqfit, X, Y, x)
end


@kernel function _lsqfitatpos_kernel4!(
    Y_est::AbstractArray{<:Real},
    @Const(A_slqfit::AbstractMatrix{<:Real}),
    X_axis::AbstractRange{<:RealQuantity},
    @Const(flat_Ys::AbstractArray{<:Real}),
    @Const(X_pos::Union{RealQuantity, AbstractVector{<:RealQuantity}})
)
   idxs = @index(Global, NTuple)
   _kbc_setindex!(
        Y_est,
        _lsqfitatpos_impl(A_slqfit, X_axis, _kbc_view(flat_Ys, :, idxs...), _kbc_getindex(X_pos, idxs...)),
        idxs...
   )
end

function bc_lsqfitat!(
    Y_est::AbstractArray{<:Real},
    degree::Integer,
    n::Integer,
    X_axis::AbstractRange{<:RealQuantity},
    Ys::_BC_RQ_AosAs,
    X_pos::_BC_RQs
)
    flat_Ys, _ = _kbc_flatview(Ys)
    @argcheck isempty(Base.tail(axes(flat_Ys))) || axes(Y_est) == Base.tail(axes(flat_Ys))
    @argcheck axes(X_axis,1) == axes(flat_Ys)[1]
    @argcheck isempty(axes(X_pos)) || axes(Y_est) == axes(X_pos)

    # ToDo: Try to avoid additional copy:
    A_slqfit_tmp = _lsq_fit_matrix(0:(n - 1), degree)
    A_slqfit = _to_same_device_as(flat_Ys, A_slqfit_tmp)

    dev = KernelAbstractions.get_device(flat_Ys)
    kernel! = _lsqfitatpos_kernel4!(dev, _ka_threads(dev)...)
    evt = kernel!(Y_est, A_slqfit, X_axis, flat_Ys, X_pos, ndrange=_kbc_size(size(Y_est)))
    wait(evt)
    return _kbc_result(Y_est)
end
