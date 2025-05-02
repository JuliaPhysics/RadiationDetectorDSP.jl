# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using RadiationDetectorDSP: _uniqueelem
using RadiationDetectorDSP: _inneraxes

@testset "array_utils" begin
    @testset "uniqueleem_tests" begin
        @test _uniqueelem([1,1,1,1,1]) == 1
        @test _uniqueelem(["hello","hello", "hello"]) == "hello"
        @test _uniqueelem(["A", "A","A" ]) == "A"
        @test _uniqueelem([true, true, true ]) == true

        @test_throws ArgumentError _uniqueelem([5, 6 ,7 ,3 ,5]) == 5
        @test_throws ArgumentError _uniqueelem(["hello", "bye"]) == "hello"
        @test_throws ArgumentError _uniqueelem([true, true, false]) == true
    end



end


@testset "_inneraxes_tests" begin
    @test _inneraxes([[1 2; 3 4]]) == (Base.OneTo(2), Base.OneTo(2))
    @test _inneraxes([ [1, 2], [3, 4], [5, 6], [7, 8] ]) == (Base.OneTo(2),)
    @test _inneraxes([[1 2; 3 4], [5 6; 7 8]]) == (Base.OneTo(2), Base.OneTo(2))
    @test _inneraxes([[1, 2]]) == (Base.OneTo(2),)
    @test _inneraxes(AbstractVector{Int}[]) == (Base.OneTo(0),)
    @test _inneraxes([Int[], Int[]]) == (Base.OneTo(0),)

    @test_throws DimensionMismatch _inneraxes([ [1, 2], [3] ])
    @test_throws MethodError _inneraxes([1, 2])

end

