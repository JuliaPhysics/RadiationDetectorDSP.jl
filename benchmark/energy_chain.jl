# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

# Benchmark of an energy reconstruction chain (baseline subtraction, pole-zero
# correction, trapezoidal filter, maximum) on a batch of charge waveforms, on
# the CPU, with CUDA.jl and with Reactant.
#
# Usage: julia --threads=16 --project=<env> benchmark/energy_chain.jl [n_waveforms...]

ENV["XLA_REACTANT_GPU_MEM_FRACTION"] = get(ENV, "XLA_REACTANT_GPU_MEM_FRACTION", "0.5")

using RadiationDetectorDSP, ArraysOfArrays
using CUDA, Reactant
import cuFFT
using LinearAlgebra, Random, Statistics, Printf

const RD = RadiationDetectorDSP

Reactant.set_default_backend("gpu")
CUDA.allowscalar(false)

kernels(flt, signals) = flt.(signals)
arrayops(flt, signals) = RD.bc_rdfilt_arrayops(flt, signals)

# Filtered signals may come back in row-major memory layout:
# (as a matrix, since `vec` of a Reactant array is a lazy wrapper)
colmax(Y::AbstractMatrix) = maximum(Y, dims = 1)
colmax(Y::Transpose) = maximum(parent(Y), dims = 2)

struct EnergyChain{F1,F2}
    n_baseline::Int
    pz::F1
    trap::F2
end

# Batch formulation: `filt` is `kernels` or `arrayops`.
function energies(chain::EnergyChain, signals::ArrayOfSimilarVectors, filt)
    X = flatview(signals)
    baseline = sum(view(X, 1:chain.n_baseline, :), dims = 1) ./ chain.n_baseline
    shifted = ArrayOfSimilarVectors(X .- baseline)
    colmax(flatview(filt(chain.trap, filt(chain.pz, shifted))))
end

energies_kernels(chain, signals) = energies(chain, signals, kernels)
energies_arrayops(chain, signals) = energies(chain, signals, arrayops)

# Signal by signal with thread-local buffers, the cache-friendly way to run the chain on a CPU:
function energies_loop(chain::EnergyChain, signals::ArrayOfSimilarVectors)
    X = flatview(signals)
    T = float(eltype(X))
    si = SamplingInfo{T}(axes(X, 1))
    fi_pz = fltinstance(chain.pz, si)
    fi_trap = fltinstance(chain.trap, si)
    n_bl = chain.n_baseline
    E = Matrix{T}(undef, 1, size(X, 2))
    chunks = Iterators.partition(axes(X, 2), cld(size(X, 2), Threads.nthreads()))
    @sync for chunk in chunks
        Threads.@spawn begin
            x = Vector{T}(undef, size(X, 1))
            y = similar(x)
            z = Vector{T}(undef, RD.flt_output_length(fi_trap))
            for j in chunk
                x_j = view(X, :, j)
                x .= x_j .- sum(view(x_j, 1:n_bl)) / n_bl
                rdfilt!(y, fi_pz, x)
                rdfilt!(z, fi_trap, y)
                E[j] = maximum(z)
            end
        end
    end
    return E
end

# `Reactant.synchronize` silently does nothing for array wrappers:
function reactant_sync(Y)
    Y isa Reactant.ConcretePJRTArray || throw(ArgumentError("Can't synchronize a $(nameof(typeof(Y)))"))
    Reactant.synchronize(Y)
end

# See filter_backends.jl on `frameworks`.
function approaches(X::Matrix; frameworks = (:cpu, :cuda, :reactant))
    result = Pair{String,Any}[]
    if :cpu in frameworks
        sig_cpu = ArrayOfSimilarVectors(X)
        push!(result,
            "CPU, signal by signal" => (energies_loop, sig_cpu, identity, nothing),
            "CPU, kernels" => (energies_kernels, sig_cpu, identity, nothing),
        )
    end
    if :cuda in frameworks
        sig_gpu = ArrayOfSimilarVectors(CuArray(X))
        cuda_sync = _ -> CUDA.synchronize()
        push!(result,
            "CUDA.jl, kernels" => (energies_kernels, sig_gpu, cuda_sync, nothing),
            "CUDA.jl, array ops" => (energies_arrayops, sig_gpu, cuda_sync, nothing),
        )
    end
    if :reactant in frameworks
        sig_r = ArrayOfSimilarVectors(Reactant.to_rarray(X))
        push!(result,
            "Reactant, CUDA kernels" => (energies_kernels, sig_r, reactant_sync, (; raise = false)),
            "Reactant, raised kernels" => (energies_kernels, sig_r, reactant_sync, (; raise = true)),
            "Reactant, array ops" => (energies_arrayops, sig_r, reactant_sync, (; raise = false)),
        )
    end
    return result
end

function compiled(f, chain, signals, opts)
    opts === nothing && return f, 0.0
    t = @elapsed fc = opts.raise ? (Reactant.@compile raise = true f(chain, signals)) : (Reactant.@compile f(chain, signals))
    return fc, t
end

# See filter_backends.jl on the garbage collection between the calls.
function timings(f, chain, signals, sync; n_warmup = 3, n = 20)
    t_warmup = minimum(1:n_warmup) do _
        t = @elapsed sync(f(chain, signals))
        GC.gc(false)
        t
    end
    # Fewer repetitions for calls that take long:
    map(1:(t_warmup > 0.5 ? 5 : n)) do _
        t0 = time_ns()
        sync(f(chain, signals))
        t = (time_ns() - t0) / 1e6
        GC.gc(false)
        t
    end
end

function transfer_times(X::Matrix; n = 5, frameworks = (:cuda, :reactant))
    t_cuda = :cuda in frameworks ? minimum(_ -> (t = (@elapsed CUDA.@sync CuArray(X)); GC.gc(); t), 1:n) : NaN
    t_reactant = :reactant in frameworks ? minimum(_ -> (t = (@elapsed reactant_sync(Reactant.to_rarray(X))); GC.gc(); t), 1:n) : NaN
    return (cuda = 1e3 * t_cuda, reactant = 1e3 * t_reactant)
end

function run_benchmarks(X::Matrix, chain::EnergyChain; frameworks = (:cpu, :cuda, :reactant), kwargs...)
    run_benchmarks(approaches(X; frameworks = frameworks), X, chain; kwargs...)
end

# `backends` from `approaches`, `X` the host data for the reference results.
function run_benchmarks(backends::AbstractVector{<:Pair}, X::Matrix, chain::EnergyChain; skip = String[], label = "")
    results = NamedTuple[]
    n_samples, n_waveforms = size(X)
    E_ref = vec(energies_loop(chain, ArrayOfSimilarVectors(X)))
    # A failed kernel launch can leave Reactant's GPU context unusable:
    reactant_ok = true
    for (aname, (f, signals, sync, opts)) in backends
        aname in skip && continue
        row = try
            startswith(aname, "Reactant") && !reactant_ok && error("skipped after an earlier Reactant failure")
            fc, t_compile = compiled(f, chain, signals, opts)
            E = vec(Array(fc(chain, signals)))
            ts = timings(fc, chain, signals, sync)
            (label = label, approach = aname, n_samples = n_samples, n_waveforms = n_waveforms, t_min = minimum(ts), t_median = median(ts),
             t_compile = t_compile, rel_max_dev = Float64(maximum(abs, E .- E_ref) / maximum(abs, E_ref)), error = "")
        catch err
            (label = label, approach = aname, n_samples = n_samples, n_waveforms = n_waveforms, t_min = NaN, t_median = NaN,
             t_compile = NaN, rel_max_dev = NaN, error = first(replace(sprint(showerror, err), '\n' => ' '), 200))
        end
        push!(results, row)
        startswith(aname, "Reactant") && !(row.rel_max_dev < 1e-3) && (reactant_ok = false)
        @printf("%-26s %9.3f ms (median %9.3f ms)  %8.2f M waveforms/s  dev %.1e  compile %5.1f s  %s\n", aname, row.t_min, row.t_median,
                n_waveforms / row.t_min / 1e3, row.rel_max_dev, row.t_compile, row.error)
        flush(stdout)
        GC.gc()
    end
    return results, E_ref
end

function write_csv(path, results)
    open(path, "w") do io
        println(io, "data,approach,n_samples,n_waveforms,t_min_ms,t_median_ms,t_compile_s,rel_max_dev,error")
        for r in results
            println(io, join(('"' * r.label * '"', '"' * r.approach * '"', r.n_samples, r.n_waveforms, r.t_min, r.t_median,
                              r.t_compile, r.rel_max_dev, '"' * replace(r.error, '"' => '\'') * '"'), ","))
        end
    end
end

function synthetic(n_waveforms::Integer; n_samples = 1024, tau = 128f0)
    r = RD.rand_waveforms(Xoshiro(1), Float32, n_samples, n_waveforms; tau_fall = (tau, tau))
    chain = EnergyChain(r.i_trigger - 50, InvCRFilter(tau), TrapezoidalChargeFilter(100, 60))
    return r, chain
end

if abspath(PROGRAM_FILE) == @__FILE__
    sizes = isempty(ARGS) ? [2^15] : parse.(Int, ARGS)
    results = reduce(vcat, [first(run_benchmarks(synthetic(n)[1].waveforms, synthetic(n)[2]; label = "synthetic")) for n in sizes])
    write_csv(joinpath(@__DIR__, "energy_chain.csv"), results)
end
