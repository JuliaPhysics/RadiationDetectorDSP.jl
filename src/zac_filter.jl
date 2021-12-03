# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    struct ZeroAreaCuspFilter <: AbstractRadFIRFilter

Zero area cusp (ZAC) filter.

For the definition the filter and a discussion of the filter properties, see
["Improvement of the energy resolution via an optimized digital signal processing in GERDA Phase I", Eur. Phys. J. C 75, 255 (2015)](https://doi.org/10.1140/epjc/s10052-015-3409-6).

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)

A sharp step on the input will result in a trapezoid with rise time and fall
time `avgtime` and a flat top of length `gaptime`.
"""
@with_kw struct ZeroAreaCuspFilter{
    TT <: RealQuantity
} <: AbstractRadNonlinearFilter
    "total length of the filter (L)"
    length::TT = 100

    "equivalent of shaping time (τₛ)"
    tau::TT = 20,

    "length of flat top (FT)"
    toplen::TT = 10

end
export ZeroAreaCuspFilter

#!!!!! ToDo: Implement ZeroAreaCuspFilter
