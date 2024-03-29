# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    struct PulseGenerator <: AbstractRadNonlinearFilter

Add a pulse with (optional) linear rise and fall to the input.

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
@with_kw struct PulseGenerator{
    TT <: RealQuantity,
    TA <: RealQuantity
} <: AbstractRadNonlinearFilter
    "start delay of pulse"
    delay::TT = 0,

    "rise time of pulse"
    risetime::TT = 0,

    "length of pulse"
    pulselen::TT = 1,

    "fall time of pulse"
    falltime::TT = risetime,

    "amplitude of pulse"
    amplitude::TA = 1

    # ToDo: Add rise/fall functional form, support linear, exponential, etc.?
end
export PulseGenerator

# ToDo: Implement PulseGenerator



"""
    struct StepGenerator <: AbstractRadNonlinearFilter

Add a pulse with (optional) linear rise and fall to the input.

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
@with_kw struct StepGenerator{
    TT <: RealQuantity,
    TA <: RealQuantity
} <: AbstractRadNonlinearFilter
    "start delay of step"
    delay::TT = 0,

    "rise time of step"
    risetime::TT = 0,

    "amplitude of step"
    amplitude::TA = 1

    # ToDo: Add rise functional form, support linear, exponential, etc.?
end
export StepGenerator

# ToDo: Implement StepGenerator



"""
    struct WhiteNoiseGenerator <: AbstractRadNonlinearFilter

Add IID random noise to the input.

Uses a deterministic RNG initialized from hash of the input.

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
@with_kw struct WhiteNoiseGenerator <: AbstractRadNonlinearFilter
    "noise level"
    level::Float64 = 1.0
end
export WhiteNoiseGenerator

# ToDo: Implement WhiteNoiseGenerator


"""
    struct ShapedNoiseGenerator <: AbstractRadNonlinearFilter

Add shaped noise to the input.

Uses a deterministic RNG initialized from hash of the input.

Constructors:

* ```$(FUNCTIONNAME)(; fields...)```

Fields:

$(TYPEDFIELDS)
"""
@with_kw struct ShapedNoiseGenerator{F} <: AbstractRadNonlinearFilter
    "spectrum transfer function"
    specfunc::F = inv
    "noise level"
    level::Float64 = 1.0
end
export ShapedNoiseGenerator

# ToDo: Implement ShapedNoiseGenerator


# ToDo: Add triangular noise generator for dithering.
