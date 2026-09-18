# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using RadiationDetectorDSP: _lsq_fit_matrix, _lsqfitatpos_single_impl, _smoothstep, lsqfitatpos

@testset "lsqfitatpos" begin
    @testset "uses smoothstep window blending" begin
        degree = 1
        n = 3
        X = 1.0:1.0:10.0
        Y = [0.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        x = 5.25

        A_lsqfit = _lsq_fit_matrix(0:(n - 1), degree)
        xi = (x - first(X)) / step(X) + firstindex(X)
        from = floor(Integer, xi) - div(n - 1, 2)
        until = from + n - 1

        y_l = _lsqfitatpos_single_impl(A_lsqfit, view(Y, from:until), xi - from)
        y_r = _lsqfitatpos_single_impl(A_lsqfit, view(Y, (from + 1):(until + 1)), xi - (from + 1))
        w_linear = mod(xi - (from + 1), 1)
        w_smooth = _smoothstep(w_linear)

        @test lsqfitatpos(degree, n, X, Y, x) ≈ (1 - w_smooth) * y_l + w_smooth * y_r
        @test lsqfitatpos(degree, n, X, Y, x) != (1 - w_linear) * y_l + w_linear * y_r
    end
end
