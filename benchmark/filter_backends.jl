# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

# Benchmark of batched filtering on synthetic charge waveforms, on the CPU,
# with CUDA.jl and with Reactant, for the filter implementations as
# KernelAbstractions kernels and as whole-array operations.
#
# Usage: julia --threads=16 --project=<env> benchmark/filter_backends.jl [n_waveforms...]
# The environment needs RadiationDetectorDSP, ArraysOfArrays, CUDA, cuFFT and Reactant.

ENV["XLA_REACTANT_GPU_MEM_FRACTION"] = get(ENV, "XLA_REACTANT_GPU_MEM_FRACTION", "0.5")

using RadiationDetectorDSP, ArraysOfArrays
using CUDA, Reactant
import cuFFT
using LinearAlgebra, Random, Statistics, Printf

const RD = RadiationDetectorDSP

Reactant.set_default_backend("gpu")
CUDA.allowscalar(false)

const n_samples = 1024
const tau = 128f0

const filters = [
    "TrapezoidalChargeFilter (100, 60)" => TrapezoidalChargeFilter(100, 60),
    "InvCRFilter (first-order IIR)" => InvCRFilter(tau),
    "CRFilter (first-order IIR)" => CRFilter(tau),
    "BiquadFilter (second-order IIR)" => BiquadFilter((0.2f0, 0.15f0, 0.3f0), (-0.8f0, 0.4f0)),
    "SavitzkyGolayFilter (direct conv., 15 taps)" => SavitzkyGolayFilter(15, 3, 0),
    "moving average (direct conv., 64 taps)" => ConvolutionFilter(DirectConvolution(), fill(1f0 / 64, 64)),
    "CUSPChargeFilter (FFT conv., 400 taps)" => CUSPChargeFilter(sigma = 100f0, toplen = 20f0, tau = tau, length = 400f0, beta = 100f0),
    "ZACChargeFilter (FFT conv., 400 taps)" => ZACChargeFilter(sigma = 100f0, toplen = 20f0, tau = tau, length = 400f0, beta = 100f0),
]

kernels(flt, signals) = flatview(flt.(signals))
arrayops(flt, signals) = flatview(RD.bc_rdfilt_arrayops(flt, signals))

unwrapped(Y) = Y isa Transpose ? parent(Y) : Y
on_host(Y) = Y isa Transpose ? permutedims(Array(parent(Y))) : Array(Y)

# `Reactant.synchronize` silently does nothing for array wrappers:
function reactant_sync(Y)
    Y isa Reactant.ConcretePJRTArray || throw(ArgumentError("Can't synchronize a $(nameof(typeof(Y)))"))
    Reactant.synchronize(Y)
end

# name => (function, input, synchronize, compile options or nothing)
function approaches(X::Matrix)
    sig_cpu = ArrayOfSimilarVectors(X)
    sig_gpu = ArrayOfSimilarVectors(CuArray(X))
    sig_r = ArrayOfSimilarVectors(Reactant.to_rarray(X))
    cuda_sync = _ -> CUDA.synchronize()
    r_sync = Y -> reactant_sync(unwrapped(Y))
    return [
        "CPU, kernels" => (kernels, sig_cpu, identity, nothing),
        "CPU, array ops" => (arrayops, sig_cpu, identity, nothing),
        "CUDA.jl, kernels" => (kernels, sig_gpu, cuda_sync, nothing),
        "CUDA.jl, array ops" => (arrayops, sig_gpu, cuda_sync, nothing),
        "Reactant, CUDA kernels" => (kernels, sig_r, r_sync, (; raise = false)),
        "Reactant, raised kernels" => (kernels, sig_r, r_sync, (; raise = true)),
        "Reactant, array ops" => (arrayops, sig_r, r_sync, (; raise = false)),
    ]
end

function compiled(f, flt, signals, opts)
    opts === nothing && return f, 0.0
    t = @elapsed fc = opts.raise ? (Reactant.@compile raise = true f(flt, signals)) : (Reactant.@compile f(flt, signals))
    return fc, t
end

# CUDA.jl only reuses the memory of arrays that have been garbage-collected, and getting
# fresh device memory takes much longer than the filters themselves. So the garbage
# collector runs between the calls, outside of the timed region.
function timings(f, flt, signals, sync; n_warmup = 3, n = 20)
    for _ in 1:n_warmup
        sync(f(flt, signals))
        GC.gc(false)
    end
    map(1:n) do _
        t0 = time_ns()
        sync(f(flt, signals))
        t = (time_ns() - t0) / 1e6
        GC.gc(false)
        t
    end
end

# `skip` may contain approach names and `(filter name, approach name)` tuples.
function run_benchmarks(n_waveforms::Integer; skip = [], filters = filters)
    X = RD.rand_waveforms(Xoshiro(1), Float32, n_samples, n_waveforms; tau_fall = (tau, tau)).waveforms
    results = NamedTuple[]
    for (fname, flt) in filters
        Y_ref = nothing
        for (aname, (f, signals, sync, opts)) in approaches(X)
            (aname in skip || (fname, aname) in skip) && continue
            row = try
                fc, t_compile = compiled(f, flt, signals, opts)
                Y = on_host(fc(flt, signals))
                Y_ref === nothing && (Y_ref = Y)
                ts = timings(fc, flt, signals, sync)
                (filter = fname, approach = aname, n_waveforms = n_waveforms, t_min = minimum(ts), t_median = median(ts),
                 t_compile = t_compile, max_dev = Float64(maximum(abs, Y .- Y_ref)), scale = Float64(maximum(abs, Y_ref)), error = "")
            catch err
                (filter = fname, approach = aname, n_waveforms = n_waveforms, t_min = NaN, t_median = NaN,
                 t_compile = NaN, max_dev = NaN, scale = NaN, error = first(replace(sprint(showerror, err), '\n' => ' '), 200))
            end
            push!(results, row)
            @printf("%-44s %-26s %9.3f ms (median %9.3f ms)  dev %.1e  %s\n", fname, aname, row.t_min, row.t_median, row.max_dev / row.scale, row.error)
            flush(stdout)
            GC.gc()
        end
    end
    return results
end

function write_csv(path, results)
    open(path, "w") do io
        println(io, "filter,approach,n_samples,n_waveforms,t_min_ms,t_median_ms,t_compile_s,rel_max_dev,error")
        for r in results
            println(io, join(('"' * r.filter * '"', '"' * r.approach * '"', n_samples, r.n_waveforms, r.t_min, r.t_median,
                              r.t_compile, r.max_dev / r.scale, '"' * replace(r.error, '"' => '\'') * '"'), ","))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    sizes = isempty(ARGS) ? [2^15] : parse.(Int, ARGS)
    results = reduce(vcat, [run_benchmarks(n) for n in sizes])
    write_csv(joinpath(@__DIR__, "filter_backends.csv"), results)
end
