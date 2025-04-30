# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

using RadiationDetectorDSP
using Test

using RadiationDetectorDSP: _uniqueelem

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
