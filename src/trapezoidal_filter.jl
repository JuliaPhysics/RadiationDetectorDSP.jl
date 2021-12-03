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
@with_kw struct TrapezoidalChargeFilter{
    TT <: RealQuantity
} <: AbstractRadNonlinearFilter
    "averaging time"
    avgtime::TT = 20,

    "gap time"
    gaptime::TT = 10
end
export TrapezoidalChargeFilter

#!!!!! ToDo: Implement TrapezoidalChargeFilter
