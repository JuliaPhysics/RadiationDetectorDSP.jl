var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/#Modules","page":"API","title":"Modules","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Order = [:module]","category":"page"},{"location":"api/#Types-and-constants","page":"API","title":"Types and constants","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Order = [:type, :constant]","category":"page"},{"location":"api/#Functions-and-macros","page":"API","title":"Functions and macros","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Order = [:macro, :function]","category":"page"},{"location":"api/#Documentation","page":"API","title":"Documentation","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [RadiationDetectorDSP]\nOrder = [:module, :type, :constant, :macro, :function]","category":"page"},{"location":"api/#RadiationDetectorDSP.AbstractRadFIRFilter","page":"API","title":"RadiationDetectorDSP.AbstractRadFIRFilter","text":"abstract type AbstractRadFIRFilter <: AbstractRadLinearFilter\n\nAbstract type for FIR filters.\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.AbstractRadIIRFilter","page":"API","title":"RadiationDetectorDSP.AbstractRadIIRFilter","text":"abstract type AbstractRadIIRFilter <: AbstractRadSigFilter{LinearFiltering}\n\nAbstract type for IIR filters.\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.AbstractRadSigFilter","page":"API","title":"RadiationDetectorDSP.AbstractRadSigFilter","text":"abstract type AbstractRadSigFilter{FT<:FilteringType} <: Function\n\nAbstract type for signal filters.\n\nFilters are callable as (flt::AbstractRadSigFilter)(input) and come with specialized broadcasting.\n\nSubtypes of AbstractRadSigFilter must implement\n\nfltinstance(flt::AbstractRadSigFilter, si::SamplingInfo)::[`AbstractRadSigFilterInstance`](@ref)\n\nInvertible filters should also implement\n\nInverseFunctions.inverse(flt::SomeFilter)\n\nNote that while a filter may have an inverse, it may, depending on the filter paramters, be very unstable in the presence of additional noise. Filters with a high-pass characteristic pass high-frequency noise, so their inverses pass such noise as well without amplifying it (substantially). Filters with a low-pass characteristic, on the other hand, attentuate high-frequency noise, so their inverses amplify such noise and are typically not useful to deconvolve signals in practical applications.\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.AbstractRadSigFilterInstance","page":"API","title":"RadiationDetectorDSP.AbstractRadSigFilterInstance","text":"abstract type AbstractRadSigFilterInstance{FT<:FilteringType}\n\nAbstract type for signal filter instances. Filter instances are specilized to a specific length and numerical type of input and output.\n\nFilter instances are callable as (fi::SomeFilterInstance)(input) and come with specialized broadcasting.\n\nSubtypes of AbstractRadSigFilterInstance must implement\n\nRadiationDetectorDSP.rdfilt!(output, fi::SomeFilterInstance, input)\nRadiationDetectorDSP.flt_output_smpltype(fi::SomeFilterInstance)\nRadiationDetectorDSP.flt_input_smpltype(fi::SomeFilterInstance)\nRadiationDetectorDSP.flt_output_length(fi::SomeFilterInstance)\nRadiationDetectorDSP.flt_input_length(fi::SomeFilterInstance)\nRadiationDetectorDSP.flt_output_time_axis(fi::SomeFilterInstance, time::AbstractVector{<:RealQuantity})\n\nInvertible filter instances should implement\n\nInverseFunctions.inverse(fi::SomeFilterInstance)\n\nDefault methods are implemented for\n\nRadiationDetectorDSP.rdfilt(fi::AbstractRadSigFilterInstance, x::AbstractSamples)\nRadiationDetectorDSP.rdfilt(fi::AbstractRadSigFilterInstance, wf::RDWaveform)\nRadiationDetectorDSP.bc_rdfilt(fi::AbstractRadSigFilterInstance, inputs)\n\nThe default methods that operate on RadiationDetectorSignals.RDWaveforms require RadiationDetectorDSP.flt_output_time_axis.\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.AbstractSamples","page":"API","title":"RadiationDetectorDSP.AbstractSamples","text":"const AbstractSamples{T<:RealQuantity} = AbstractVector{T}\n\nA vector of signal samples.\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.ArrayOfSimilarSamples","page":"API","title":"RadiationDetectorDSP.ArrayOfSimilarSamples","text":"const ArrayOfSimilarSamples{T<:RealQuantity} = ArrayOfSimilarVectors{T}\n\nAn array of similar sample vectors.\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.BiquadFilter","page":"API","title":"RadiationDetectorDSP.BiquadFilter","text":"struct BiquadFilter{T<:RealQuantity} <: AbstractRadIIRFilter\n\nA biquad filter.\n\nConstructors:\n\nBiquadFilter(fields...)\n\nFields:\n\nb_012::Tuple{T, T, T} where T<:(Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real): Coefficients b0 to b2\na_12::Tuple{T, T} where T<:(Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real): Coefficients a1 to a2, a_0 equals one implicitly\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.CPUNormAdaptor","page":"API","title":"RadiationDetectorDSP.CPUNormAdaptor","text":"RadiationDetectorDSP.CPUNormAdaptor\n\nTo be used with Adapt.adapt.\n\nAdapt.adapt(RadiationDetectorDSP.CPUNormAdaptor, x) adapts x to reside on the CPU and tries to ensure that arrays are stored in column-major order.\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.CRFilter","page":"API","title":"RadiationDetectorDSP.CRFilter","text":"struct CRFilter <: AbstractRadIIRFilter\n\nA first-order CR highpass filter.\n\nThe inverse filter is InvCRFilter, this is typically stable even in the presence of additional noise. This is because a CR filter passes high-frequency noise and so it's inverse passes such noise as well without amplifying it.\n\nConstructors:\n\nCRFilter(fields...)\n\nFields:\n\ncr::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: CR time constant\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.ConvolutionFilter","page":"API","title":"RadiationDetectorDSP.ConvolutionFilter","text":"struct ConvolutionFilter{T<:RealQuantity} <: AbstractRadFIRFilter\n\nA FIR filter defined by it's filter taps, applied via convolution with the input signal.\n\nConstructors:\n\nConvolutionFilter(fields...)\n\nFields:\n\nmethod::RadiationDetectorDSP.ConvolutionMethod: Convolution method\ncoeffs::AbstractVector{T} where T<:(Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real): Filter taps\noffset::Int64: Time axis offset\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.ConvolutionMethod","page":"API","title":"RadiationDetectorDSP.ConvolutionMethod","text":"abstract type ConvolutionMethod\n\nIndended as a type parameter to designate the behavior of a filter as linear or nonlinear.\n\nSubtypes are DirectConvolution and FFTConvolution.\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.DNIMethod","page":"API","title":"RadiationDetectorDSP.DNIMethod","text":"abstract type DNIMethod\n\nAbstract type for denoising and interpolation methods.\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.DifferentiatorFilter","page":"API","title":"RadiationDetectorDSP.DifferentiatorFilter","text":"struct DifferentiatorFilter <: AbstractRadIIRFilter\n\nAn integrator filter. It's inverse is IntegratorFilter.\n\nConstructors:\n\nDifferentiatorFilter(fields...)\n\nFields:\n\ngain::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: Filter gain\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.DirectConvolution","page":"API","title":"RadiationDetectorDSP.DirectConvolution","text":"DirectConvolution() isa ConvolutionMethod\n\nCompute filter convolutions directly, without FFT.\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.FFTConvolution","page":"API","title":"RadiationDetectorDSP.FFTConvolution","text":"FFTConvolution() isa ConvolutionMethod\n\nCompute filter convolutions via FFT.\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.FilteringType","page":"API","title":"RadiationDetectorDSP.FilteringType","text":"abstract type FilteringType\n\nIndended as a type parameter to designate the behavior of a filter as linear or nonlinear.\n\nSubtypes are LinearFiltering and NonlinearFiltering.\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.FirstOrderIIR","page":"API","title":"RadiationDetectorDSP.FirstOrderIIR","text":"struct FirstOrderIIR{T<:RealQuantity} <: AbstractRadIIRFilter\n\nA biquad filter.\n\nConstructors:\n\nFirstOrderIIR(fields...)\n\nFields:\n\nb_01::Tuple{T, T} where T<:(Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real): Coefficients b0 to b1\na_1::Tuple{T} where T<:(Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real): Coefficient a1, a0 equals one implicitly\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.IntegratorCRFilter","page":"API","title":"RadiationDetectorDSP.IntegratorCRFilter","text":"struct IntegratorCRFilter <: AbstractRadIIRFilter\n\nA modified CR-filter. The filter has an inverse.\n\nConstructors:\n\nIntegratorCRFilter(fields...)\n\nFields:\n\ngain::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: Filter gain\ncr::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: CR time constant\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.IntegratorFilter","page":"API","title":"RadiationDetectorDSP.IntegratorFilter","text":"struct IntegratorFilter <: AbstractRadIIRFilter\n\nAn integrator filter. It's inverse is DifferentiatorFilter.\n\nConstructors:\n\nIntegratorFilter(fields...)\n\nFields:\n\ngain::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: Filter gain\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.IntegratorModCRFilter","page":"API","title":"RadiationDetectorDSP.IntegratorModCRFilter","text":"struct IntegratorModCRFilter <: AbstractRadIIRFilter\n\nA modified CR-filter. The filter has an inverse.\n\nConstructors:\n\nIntegratorModCRFilter(fields...)\n\nFields:\n\ngain::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: Filter gain\ncr::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: CR time constant\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.Intersect","page":"API","title":"RadiationDetectorDSP.Intersect","text":"struct Intersect <: Function\n\nFinds the intersects of a Y with a threshold\n\nConstructors:\n\nIntersect(; fields...)\n\nFields:\n\nmintot::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: minimum time-over-threshold\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.InvCRFilter","page":"API","title":"RadiationDetectorDSP.InvCRFilter","text":"struct InvCRFilter <: AbstractRadIIRFilter\n\nInverse of CRFilter.\n\nConstructors:\n\nInvCRFilter(fields...)\n\nFields:\n\ncr::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: CR time constant\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.InvModCRFilter","page":"API","title":"RadiationDetectorDSP.InvModCRFilter","text":"struct InvModCRFilter <: AbstractRadIIRFilter\n\nInverse of ModCRFilter.\n\nConstructors:\n\nInvModCRFilter(fields...)\n\nFields:\n\ncr::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: CR time constant\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.InvRCFilter","page":"API","title":"RadiationDetectorDSP.InvRCFilter","text":"struct InvRCFilter <: AbstractRadIIRFilter\n\nInverse of RCFilter.\n\nConstructors:\n\nInvRCFilter(fields...)\n\nFields:\n\nrc::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: RC time constant\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.LinearFiltering","page":"API","title":"RadiationDetectorDSP.LinearFiltering","text":"abstract type LinearFiltering <: FilteringType\n\nWhen used as a type parameter value, marks linear behavior of a filter.\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.ModCRFilter","page":"API","title":"RadiationDetectorDSP.ModCRFilter","text":"struct ModCRFilter <: AbstractRadIIRFilter\n\nA first-order CR highpass filter, modified for full-amplitude step-signal response.\n\nThe resonse of the standard digital CRFilter will not recover the full amplitude of a digital step stignal since a step from one sample to the still has a finite rise time. This version of a CR filter compensates for this loss in amplitude, so it effectively treats a step as having\n\nThe inverse filter is InvModCRFilter, this is typically stable even in the presence of additional noise (see CRFilter).\n\nConstructors:\n\nModCRFilter(fields...)\n\nFields:\n\ncr::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: CR time constant\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.NonlinearFiltering","page":"API","title":"RadiationDetectorDSP.NonlinearFiltering","text":"abstract type NonlinearFiltering <: FilteringType\n\nWhen used as a type parameter value, marks non linear behavior of a filter.\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.PolynomialDNI","page":"API","title":"RadiationDetectorDSP.PolynomialDNI","text":"struct PolynomialDNI <: DNIMethod\n\nPolynomial denoising and interpolation method.\n\nOperates in a similar way as a Savitzky-Golay filter, but interpolates as well.\n\nConstructors:\n\nPolynomialDNI(; fields...)\n\nFields:\n\ndegree::Int64: polynomial degree Default: (2, \"length\")\nlength::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.RCFilter","page":"API","title":"RadiationDetectorDSP.RCFilter","text":"struct RCFilter <: AbstractRadII>RFilter\n\nA first-order RC lowpass filter.\n\nThe inverse filter is InvCRFilter, but note that this is unstable in the presence of additional noise. As an RC filter attenuates high-frequency noise, its inverse amplifies such noise and will typically not be useful to deconvolve signals in practical applications.\n\nConstructors:\n\nRCFilter(fields...)\n\nFields:\n\nrc::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: RC time constant\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.SamplingInfo","page":"API","title":"RadiationDetectorDSP.SamplingInfo","text":"struct SamplingInfo{T<:RealQuantity,A<:AbstractVector{<:RealQuantity}}\n\nHolds sampling information.\n\nThe numerical type of an individual sample is T, the (time) axis is given by the axis field.\n\nConstructors:\n\nSamplingInfo{T,A}(axis)\nSamplingInfo{T}(axis)\n\nFields:\n\naxis::AbstractVector{<:Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real}\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.SavitzkyGolayFilter","page":"API","title":"RadiationDetectorDSP.SavitzkyGolayFilter","text":"struct SavitzkyGolayFilter <: AbstractRadFIRFilter\n\nA Savitzky-Golay filter.\n\nConstructors:\n\nSavitzkyGolayFilter(; fields...)\n\nFields:\n\nlength::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: filter length\ndegree::Int64: Polynomial defgree\nderivative::Int64: n-th derivative (0 for no derivative)\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.SignalEstimator","page":"API","title":"RadiationDetectorDSP.SignalEstimator","text":"struct SignalEstimator <: Function\n\nEstimates a signal at a given position x.\n\nUsage:\n\n(f::SamplesOrWaveform)(input::RDWaveform, x::RealQuantity)\n\nConstructors:\n\nSignalEstimator(; fields...)\n\nFields:\n\ndni::DNIMethod: denoising and interpolation method\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.SimpleCSAFilter","page":"API","title":"RadiationDetectorDSP.SimpleCSAFilter","text":"struct SimpleCSAFilter <: AbstractRadIIRFilter\n\nSimulates the current-signal response of a charge-sensitive preamplifier with resistive reset, the output is a charge signal.\n\nIt is equivalent to the composition\n\nCRFilter(cr = tau_decay) ∘\nIntegrator(gain = gain) ∘\nRCFilter(rc = tau_rise)\n\nand maps to a single BiquadFilter.\n\nThis filter has an inverse, but the inverse is very unstable in the presence of additional noise if tau_rise is not zero (since the inverse of an RC-filter is unstable under noise). Even if tau_rise is zero the inverse will still amplify noise (since it differentiates), so it should be used very carefully when deconvolving signals in practical applications.\n\nConstructors:\n\nSimpleCSAFilter(fields...)\n\nFields:\n\ntau_rise::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: Rise time constant\ntau_decay::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: Decay time constant\ngain::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: Gain\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.TrapezoidalChargeFilter","page":"API","title":"RadiationDetectorDSP.TrapezoidalChargeFilter","text":"struct TrapezoidalChargeFilter <: AbstractRadNonlinearFilter\n\nFilter that responds to a step signal with a trapezoidal pulse.\n\nThe filter is equivalent to two moving averages separated by a gap.\n\nConstructors:\n\nTrapezoidalChargeFilter(; fields...)\n\nFields:\n\navgtime::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: pre-rise averaging time\ngaptime::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: gap time\navgtime2::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: post-rise averaging time\n\nA sharp step on the input will result in a trapezoid with rise time and fall time avgtime and a flat top of length gaptime.\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.TruncateFilter","page":"API","title":"RadiationDetectorDSP.TruncateFilter","text":"struct TruncateFilter <: AbstractRadSigFilter{LinearFiltering}\n\nFilter that truncates the input signal.\n\nConstructors:\n\nTruncateFilter(; fields...)\n\nFields:\n\ninterval::IntervalSets.ClosedInterval{<:Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real}: interval to keep\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.ZACChargeFilter","page":"API","title":"RadiationDetectorDSP.ZACChargeFilter","text":"struct ZACChargeFilter <: AbstractRadFIRFilter\n\nZero area cusp (ZAC) filter.\n\nFor the definition the filter and a discussion of the filter properties, see \"Improvement of the energy resolution via an optimized digital signal processing in GERDA Phase I\", Eur. Phys. J. C 75, 255 (2015).\n\nConstructors:\n\nZACChargeFilter(; fields...)\n\nFields:\n\ntau::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: equivalent of shaping time (τₛ) Default: (20, \"length of flat top (FT)\")\ntoplen::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: Default: 10\nlength::Union{Unitful.Quantity{<:var\"#s11\"}, var\"#s11\"} where var\"#s11\"<:Real: total length of the filter (L) Default: 100\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.MaybeWithUnits","page":"API","title":"RadiationDetectorDSP.MaybeWithUnits","text":"const MaybeWithUnits{T<:Number} = Union{T,Quantity{<:T}}\n\nA numerical value with or without units\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.RealQuantity","page":"API","title":"RadiationDetectorDSP.RealQuantity","text":"const RealQuantity = MaybeWithUnits{<:Real}\n\nA real value with or without units.\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.SamplesOrWaveform","page":"API","title":"RadiationDetectorDSP.SamplesOrWaveform","text":"const RadiationDetectorDSP.SamplesOrWaveform{T<:RealQuantity} = Union{AbstractSamples{T},RDWaveform{<:Any,T}}\n\nA vector of signal samples or a waveform.\n\n\n\n\n\n","category":"type"},{"location":"api/#RadiationDetectorDSP.adapt_memlayout","page":"API","title":"RadiationDetectorDSP.adapt_memlayout","text":"RadiationDetectorDSP::adapt_memlayout(\n    fi::AbstractRadSigFilterInstance,\n    backend::KernelAbstractions.Backend,\n    A::AbstractArray{<:Number}\n)\n\nAdapts the memory layout of A in a suitable fashion for fi on computing device backend.\n\nReturns a row-major version of A on all backends by default, filter instance types may specialize this behavior.\n\n\n\n\n\n","category":"function"},{"location":"api/#RadiationDetectorDSP.add_rect_pulse!","page":"API","title":"RadiationDetectorDSP.add_rect_pulse!","text":"add_rect_pulse!(samples::AbstractSamples, start::Integer, pulselen::Integer, amplitude::Real = 1.0)\n\nAdd a rectangular pulse to samples.\n\n\n\n\n\n","category":"function"},{"location":"api/#RadiationDetectorDSP.bc_rdfilt","page":"API","title":"RadiationDetectorDSP.bc_rdfilt","text":"bc_rdfilt(flt::AbstractRadSigFilter, input)\nbc_rdfilt(fi::AbstractRadSigFilterInstance, input)\n\nBroadcast filter instance fi over signals input, return the filtered signals.\n\n\n\n\n\n","category":"function"},{"location":"api/#RadiationDetectorDSP.bc_rdfilt!","page":"API","title":"RadiationDetectorDSP.bc_rdfilt!","text":"bc_rdfilt!(outputs, fi::AbstractRadSigFilterInstance, inputs)\n\nBroadcast filter flt or filter instance fi over signals inputs, storing the results in outputs.\n\ninputs and outputs must be of type AbstractVector{<:AbstractSamples}.\n\nReturns outputs.\n\n\n\n\n\n","category":"function"},{"location":"api/#RadiationDetectorDSP.charge_trapflt!-Tuple{AbstractVector{<:Union{var\"#s10\", var\"#s32\"} where {var\"#s10\"<:AbstractFloat, var\"#s32\"<:(SIMD.Vec{N, <:var\"#s10\"} where N)}}, Integer, Integer}","page":"API","title":"RadiationDetectorDSP.charge_trapflt!","text":"charge_trapflt!(samples::AbstractVector{<:RealOrSIMD{<:AbstractFloat}}, navg::Integer, ngap::Integer)\n\nApply a trapezoidal FIR filter to a charge signal in samples.\n\n\n\n\n\n","category":"method"},{"location":"api/#RadiationDetectorDSP.cr_filter-Tuple{Real}","page":"API","title":"RadiationDetectorDSP.cr_filter","text":"cr_filter(CR::Real)\n\nReturn a DSP.jl-compatible CR-filter.\n\n\n\n\n\n","category":"method"},{"location":"api/#RadiationDetectorDSP.create_zac_filter-NTuple{4, Any}","page":"API","title":"RadiationDetectorDSP.create_zac_filter","text":"create_zac_filter(Nₜ, FTₜ, τₜ, Δt)\n\ncreate the zac filter according to Eur. Phys. J. C (2015) 75:255], where  Nₜ is the length of the filter, FTₜ the length of the flat top, τₜ the the  exponential decay factor and Δt the sampling time\n\n\n\n\n\n","category":"method"},{"location":"api/#RadiationDetectorDSP.crmod_filter-Tuple{Real}","page":"API","title":"RadiationDetectorDSP.crmod_filter","text":"crmod_filter(CR::Real)\n\nReturn a DSP.jl-compatible modified CR-filter.\n\n\n\n\n\n","category":"method"},{"location":"api/#RadiationDetectorDSP.dfilt!","page":"API","title":"RadiationDetectorDSP.dfilt!","text":"rdfilt!(output, fi::AbstractRadSigFilterInstance, input)\n\nApply filter flt or filter instance fi to signal input and store the filtered signal in output. Return output.\n\n\n\n\n\n","category":"function"},{"location":"api/#RadiationDetectorDSP.differentiator_filter-Tuple{Real}","page":"API","title":"RadiationDetectorDSP.differentiator_filter","text":"differentiator_filter(gain::Real)\n\nReturn a DSP.jl-compatible differentiator filter.\n\n\n\n\n\n","category":"method"},{"location":"api/#RadiationDetectorDSP.elsmplinfo","page":"API","title":"RadiationDetectorDSP.elsmplinfo","text":"smplinfo(smpls::AbstractSamples)::SamplingInfo\nsmplinfo(wf::RDWaveform{T,U})::RDWaveform\n\nGet sampling information an array of vectors of samples, resp. an array of waveform. All elements must have equal samling information.\n\n\n\n\n\n","category":"function"},{"location":"api/#RadiationDetectorDSP.flt_input_length","page":"API","title":"RadiationDetectorDSP.flt_input_length","text":"RadiationDetectorDSP.flt_input_length(fi::AbstractRadSigFilterInstance)::Integer\n\nGet the output signal length of a filter instance fi.\n\nMust be implemented for all subtypes of AbstractRadSigFilterInstance.\n\n\n\n\n\n","category":"function"},{"location":"api/#RadiationDetectorDSP.flt_input_smpltype","page":"API","title":"RadiationDetectorDSP.flt_input_smpltype","text":"RadiationDetectorDSP.flt_input_smpltype(fi::AbstractRadSigFilterInstance)\n\nGet the input sample type of a filter instance fi.\n\nMust be implemented for all subtypes of AbstractRadSigFilterInstance.\n\n\n\n\n\n","category":"function"},{"location":"api/#RadiationDetectorDSP.flt_output_length","page":"API","title":"RadiationDetectorDSP.flt_output_length","text":"RadiationDetectorDSP.flt_output_length(fi::SomeFilterInstance)::Integer\n\nGet the output signal length of filter instance fi.\n\nMust be implemented for all subtypes of AbstractRadSigFilterInstance.\n\n\n\n\n\n","category":"function"},{"location":"api/#RadiationDetectorDSP.flt_output_smpltype","page":"API","title":"RadiationDetectorDSP.flt_output_smpltype","text":"RadiationDetectorDSP.flt_output_smpltype(fi::AbstractRadSigFilterInstance)\n\nGet the output sample type for\n\na filter flt given an input sample type input_smpltype\na filter instance fi\n\nMust be implemented for all subtypes of AbstractRadSigFilter.\n\n\n\n\n\n","category":"function"},{"location":"api/#RadiationDetectorDSP.flt_output_time_axis","page":"API","title":"RadiationDetectorDSP.flt_output_time_axis","text":"RadiationDetectorDSP.flt_output_time_axis(fi::SomeFilterInstance, time::AbstractVector{<:RealQuantity})::AbstractVector{<:RealQuantity}\n\nGet the output time axis of a filter instance fi, given an input time axis time.\n\nMust be implemented for subtypes of AbstractRadSigFilter and AbstractRadSigFilterInstance only if the filter's output time axis can be computed directly from the input time axis.\n\n\n\n\n\n","category":"function"},{"location":"api/#RadiationDetectorDSP.fltinstance","page":"API","title":"RadiationDetectorDSP.fltinstance","text":"fltinstance(flt::AbstractRadSigFilter, si::SamplingInfo)::AbstractRadSigFilterInstance\n\nCreate a filter instance of the filter flt, specialized for the given input, resp. input characteristics.\n\n\n\n\n\n","category":"function"},{"location":"api/#RadiationDetectorDSP.gen_rect_pulse","page":"API","title":"RadiationDetectorDSP.gen_rect_pulse","text":"gen_rect_pulse(tracelen::Integer, start::Integer, pulselen::Integer, amplitude::Real = 1.0)\n\nGenerate a rectangular pulse.\n\n\n\n\n\n","category":"function"},{"location":"api/#RadiationDetectorDSP.integrator_cr_filter-Tuple{Real, Real}","page":"API","title":"RadiationDetectorDSP.integrator_cr_filter","text":"integrator_cr_filter(gain::Real, CR::Real)\n\nReturn a DSP.jl-compatible integrator plus CR filter.\n\n\n\n\n\n","category":"method"},{"location":"api/#RadiationDetectorDSP.integrator_crmod_filter-Tuple{Real, Real}","page":"API","title":"RadiationDetectorDSP.integrator_crmod_filter","text":"integrator_crmod_filter(gain::Real, CR::Real)\n\nReturn a DSP.jl-compatible integrator plus modified CR filter.\n\n\n\n\n\n","category":"method"},{"location":"api/#RadiationDetectorDSP.integrator_filter-Tuple{Real}","page":"API","title":"RadiationDetectorDSP.integrator_filter","text":"integrator_filter(gain::Real)\n\nReturn a DSP.jl-compatible integrator filter.\n\n\n\n\n\n","category":"method"},{"location":"api/#RadiationDetectorDSP.inv_cr_filter-Tuple{Real}","page":"API","title":"RadiationDetectorDSP.inv_cr_filter","text":"inv_cr_filter(CR::Real)\n\nReturn a DSP.jl-compatible inverse CR-filter.\n\n\n\n\n\n","category":"method"},{"location":"api/#RadiationDetectorDSP.inv_crmod_filter-Tuple{Real}","page":"API","title":"RadiationDetectorDSP.inv_crmod_filter","text":"inv_crmod_filter(CR::Real)\n\nReturn a DSP.jl-compatible inverse modified CR-filter.\n\n\n\n\n\n","category":"method"},{"location":"api/#RadiationDetectorDSP.inv_rc_filter-Tuple{Real}","page":"API","title":"RadiationDetectorDSP.inv_rc_filter","text":"inv_rc_filter(RC::Real)\n\nReturn a DSP.jl-compatible RC-filter.\n\n\n\n\n\n","category":"method"},{"location":"api/#RadiationDetectorDSP.rc_filter-Tuple{Real}","page":"API","title":"RadiationDetectorDSP.rc_filter","text":"rc_filter(RC::Real)\n\nReturn a DSP.jl-compatible RC-filter.\n\n\n\n\n\n","category":"method"},{"location":"api/#RadiationDetectorDSP.rdfilt","page":"API","title":"RadiationDetectorDSP.rdfilt","text":"rdfilt(fi::AbstractRadSigFilterInstance, input)\n\nApply filter instance fi to signal input, return the filtered signal.\n\nReturns output.\n\n\n\n\n\n","category":"function"},{"location":"api/#RadiationDetectorDSP.shift_waveform","page":"API","title":"RadiationDetectorDSP.shift_waveform","text":"shift_waveform(signal::AbstractSamples, a::RealQuantity)\nshift_waveform(wf::RDWaveform, a::RealQuantity)\n\nShifts each sample of a waveform up by a.\n\n\n\n\n\n","category":"function"},{"location":"api/#RadiationDetectorDSP.simple_csa_response_filter","page":"API","title":"RadiationDetectorDSP.simple_csa_response_filter","text":"simplecsaresponsefilter(τrise::Real, τdecay::Real, gain::Real = one(τrise))\n\nReturn a DSP.jl-compatible filter that models the response of a typical charge-sensitive amplifier (CSA).\n\n\n\n\n\n","category":"function"},{"location":"api/#RadiationDetectorDSP.smplinfo","page":"API","title":"RadiationDetectorDSP.smplinfo","text":"smplinfo(smpls::AbstractSamples)::SamplingInfo\nsmplinfo(wf::RDWaveform{T,U})::RDWaveform\n\nGet sampling information from a vector of samples, resp. a waveform.\n\n\n\n\n\n","category":"function"},{"location":"api/#RadiationDetectorDSP.zac_charge_filter_coeffs-NTuple{4, Any}","page":"API","title":"RadiationDetectorDSP.zac_charge_filter_coeffs","text":"zac_charge_filter_coeffs(FT, τₛ, Δt, T)\n\nreturn a vector representing the zac filter applicaible on a charge  signal, where FT is the length of the flat top, τₛ the filter  shaping time, Δt the sampling time and T the total length of the  filter (see Eur. Phys. J. C (2015) 75:255).\n\n\n\n\n\n","category":"method"},{"location":"api/#RadiationDetectorDSP.zac_current_filter_integral-Tuple{Any, Any, Any}","page":"API","title":"RadiationDetectorDSP.zac_current_filter_integral","text":"zac_current_filter_integral(L, FT, τₛ)\n\nCompute the integral of the unnormalised zac filter using its analytical expression, where L is the length of the polynomial part, FT the length  of the flat top part and τₛ the filter shaping time\n\n\n\n\n\n","category":"method"},{"location":"api/#RadiationDetectorDSP.zac_current_filter_shape-Union{Tuple{T}, NTuple{5, T}} where T<:AbstractFloat","page":"API","title":"RadiationDetectorDSP.zac_current_filter_shape","text":"zac_current_filter_shape(t, τₛ, L, FT, A)\n\nvalue of zac filter evaluated at t for the given set of parameters  (see zac_charge_filter_coeffs)\n\n\n\n\n\n","category":"method"},{"location":"LICENSE/#LICENSE","page":"LICENSE","title":"LICENSE","text":"","category":"section"},{"location":"LICENSE/","page":"LICENSE","title":"LICENSE","text":"using Markdown\nMarkdown.parse_file(joinpath(@__DIR__, \"..\", \"..\", \"LICENSE.md\"))","category":"page"},{"location":"#RadiationDetectorDSP.jl","page":"Home","title":"RadiationDetectorDSP.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package provides DSP algorithms for the output signal waveforms of radiation detectors (e.g. semiconductor- or scintillator-based systems).","category":"page"}]
}
