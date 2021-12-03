# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


"""
    const AbstractRadNonlinearFilter = AbstractRadSigFilter{NonlinearFiltering}

Convenience type alias, abstract linear filter.
"""
const AbstractRadNonlinearFilter = AbstractRadSigFilter{NonlinearFiltering}
export AbstractRadNonlinearFilter


# ToDo:
#
# Dither(...)
# FFTDenoise(...)
# WaveletDenoise(...)
# ...
