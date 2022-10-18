# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    const AbstractRadLinearFilter = AbstractRadSigFilter{LinearFiltering}

Convenience type alias, abstract linear filter.
"""
const AbstractRadLinearFilter = AbstractRadSigFilter{LinearFiltering}
export AbstractRadLinearFilter


"""
    const AbstractRadLinearFilterInstance = AbstractRadSigFilterInstance{LinearFiltering}

Convenience type alias, abstract linear filter instance.
"""
const AbstractRadLinearFilterInstance = AbstractRadSigFilterInstance{LinearFiltering}
export AbstractRadLinearFilterInstance
