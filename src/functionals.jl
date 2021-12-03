# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

"""
    abstract type AbstractRadSigFunctional

Abstract type for functional that compute signal summary properties.
"""
abstract type AbstractRadSigFunctional end
export AbstractRadSigFunctional


# ToDo:
#
# SignalMean(...)
# SignalVariance(...)
# SignalMaximum(...)
# SignalLinearFit(...)
# SignalThresholdIntersect(...)
# ...
