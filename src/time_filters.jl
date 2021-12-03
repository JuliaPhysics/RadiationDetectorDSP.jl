# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    struct DropFilter <: AbstractRadLinearFilter

Drop a length of input, pass another length, then drop the rest (if any).

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
@with_kw struct DropFilter{
    TT <: RealQuantity
} <: AbstractRadNonlinearFilter
    "start delay of pulse"
    droplen::TT = 0,

    "rise time of pulse"
    takelen::TT = Inf
end
export DropFilter

#!!!!! ToDo: Implement DropFilter


"""
    struct ResamplingFilter <: AbstractRadLinearFilter

Resample input.

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
@with_kw struct ResamplingFilter <: AbstractRadNonlinearFilter
    "resampling ratio, output length = ratio * input length"
    ratio::Rational{Int} = 1//1
end
export ResamplingFilter

#!!!!! ToDo: Implement ResamplingFilter
