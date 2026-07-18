# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using RadiationDetectorDSP: _get_or_view, _firstof, _secondof, _thirdof, _lastof,
    _all_except_firstof, _all_except_lastof, _firstn_of, _lastn_of,
    _split_firstn_of, _split_lastn_of, _ncolons


@testset "fast_indexing" begin
    v = [10, 20, 30, 40, 50]
    t = (10, 20, 30, 40, 50)
    nt = (a = 10, b = 20, c = 30)
    r = 10:10:50
    rf = Ref(42)

    @testset "_get_or_view" begin
        @test @inferred(_get_or_view(v, 2)) == 20
        @test @inferred(_get_or_view(v, 2:4)) == [20, 30, 40]
        @test _get_or_view(v, 2:4) isa SubArray
        @test @inferred(_get_or_view(r, 2)) == 20
        @test @inferred(_get_or_view(r, 2:4)) == 20:10:40
        @test _get_or_view(r, 2:4) isa AbstractRange
        A = reshape(1:12, 3, 4)
        @test @inferred(_get_or_view(A, 2, 3)) == A[2, 3]
        @test @inferred(_get_or_view(A, 1:2, 3)) == view(A, 1:2, 3)
    end

    @testset "element access" begin
        @test @inferred(_firstof(v)) == 10
        @test @inferred(_firstof(t)) == 10
        @test @inferred(_firstof(nt)) == 10
        @test @inferred(_firstof(rf)) == 42

        @test @inferred(_secondof(v)) == 20
        @test @inferred(_secondof(t)) == 20
        @test @inferred(_secondof(nt)) == 20
        @test_throws BoundsError _secondof(rf)

        @test @inferred(_thirdof(v)) == 30
        @test @inferred(_thirdof(t)) == 30
        @test @inferred(_thirdof(nt)) == 30
        @test_throws BoundsError _thirdof(rf)

        @test @inferred(_lastof(v)) == 50
        @test @inferred(_lastof(t)) == 50
        @test @inferred(_lastof(nt)) == 30
        @test @inferred(_lastof(rf)) == 42
    end

    @testset "_all_except_firstof / _all_except_lastof" begin
        @test _all_except_firstof(v) == [20, 30, 40, 50]
        @test @inferred(_all_except_firstof(t)) == (20, 30, 40, 50)
        @test @inferred(_all_except_firstof(nt)) == (b = 20, c = 30)
        @test @inferred(_all_except_firstof(())) == ()
        @test @inferred(_all_except_firstof(rf)) == ()

        @test _all_except_lastof(v) == [10, 20, 30, 40]
        @test @inferred(_all_except_lastof(t)) == (10, 20, 30, 40)
        @test @inferred(_all_except_lastof(nt)) == (a = 10, b = 20)
        @test @inferred(_all_except_lastof(rf)) == ()
    end

    @testset "_firstn_of / _lastn_of" begin
        @test _firstn_of(v, 3) == [10, 20, 30]
        @test _firstn_of(t, 3) == (10, 20, 30)
        @test _firstn_of(nt, 2) == (a = 10, b = 20)
        @test _firstn_of(rf, 1) == (42,)

        @test @inferred(_firstn_of(v, Val(3))) == (10, 20, 30)
        @test @inferred(_firstn_of(t, Val(3))) == (10, 20, 30)
        @test @inferred(_firstn_of(nt, Val(2))) == (a = 10, b = 20)
        @test @inferred(_firstn_of(rf, Val(1))) == (42,)

        @test _lastn_of(v, 3) == [30, 40, 50]
        @test _lastn_of(t, 3) == (30, 40, 50)
        @test _lastn_of(nt, 2) == (b = 20, c = 30)
        @test _lastn_of(rf, 1) == (42,)

        @test @inferred(_lastn_of(v, Val(3))) == (30, 40, 50)
        @test @inferred(_lastn_of(t, Val(3))) == (30, 40, 50)
        @test @inferred(_lastn_of(nt, Val(2))) == (b = 20, c = 30)
        @test @inferred(_lastn_of(rf, Val(1))) == (42,)
    end

    @testset "_split_firstn_of / _split_lastn_of" begin
        @test _split_firstn_of(v, 2) == ([10, 20], [30, 40, 50])
        @test _split_firstn_of(t, 2) == ((10, 20), (30, 40, 50))
        @test _split_firstn_of(nt, 2) == ((a = 10, b = 20), (c = 30,))
        @test _split_firstn_of(rf, 1) == ((42,), ())

        @test @inferred(_split_firstn_of(t, Val(2))) == ((10, 20), (30, 40, 50))
        @test @inferred(_split_firstn_of(nt, Val(2))) == ((a = 10, b = 20), (c = 30,))
        @test @inferred(_split_firstn_of(rf, Val(1))) == ((42,), ())

        @test _split_lastn_of(v, 2) == ([10, 20, 30], [40, 50])
        @test _split_lastn_of(t, 2) == ((10, 20, 30), (40, 50))
        @test _split_lastn_of(nt, 2) == ((a = 10,), (b = 20, c = 30))
        @test _split_lastn_of(rf, 1) == ((), (42,))

        @test @inferred(_split_lastn_of(t, Val(2))) == ((10, 20, 30), (40, 50))
        @test @inferred(_split_lastn_of(nt, Val(2))) == ((a = 10,), (b = 20, c = 30))
        @test @inferred(_split_lastn_of(rf, Val(1))) == ((), (42,))
    end

    @testset "_ncolons" begin
        @test @inferred(_ncolons(Val(0))) == ()
        @test @inferred(_ncolons(Val(3))) == (:, :, :)
    end
end
