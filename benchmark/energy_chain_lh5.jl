# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).

# Benchmark of the energy reconstruction chain of energy_chain.jl on the
# waveforms of one channel of a LEGEND raw data file.
#
# Usage: julia --threads=16 --project=<env> benchmark/energy_chain_lh5.jl <file.lh5> [channel]
# The environment additionally needs LegendHDF5IO. Without a channel name, the
# germanium detector channel with the most waveforms is used.

include(joinpath(@__DIR__, "energy_chain.jl"))

using LegendHDF5IO

function read_waveforms(filename::AbstractString, channel = nothing; column::Symbol = :waveform_presummed)
    lh5open(filename) do f
        ch = if channel === nothing
            detectors = filter(k -> occursin(r"^[BCPV]\d", string(k)), collect(keys(f)))
            detectors[argmax([length(f[k].raw) for k in detectors])]
        else
            channel
        end
        wfs = getproperty(f[ch].raw, column)[:]
        (waveforms = Float32.(flatview(wfs.signal)), dt = step(wfs.time[1]), channel = string(ch))
    end
end

# Waveforms of the germanium detector channels of several files, `n_waveforms` in total:
function read_waveforms(filenames::AbstractVector{<:AbstractString}, n_waveforms::Integer; column::Symbol = :waveform_presummed)
    X, dt, n = nothing, nothing, 0
    for filename in filenames
        n >= n_waveforms && break
        lh5open(filename) do f
            for ch in filter(k -> occursin(r"^[BCPV]\d", string(k)), collect(keys(f)))
                col = getproperty(f[ch].raw, column)
                m = min(length(col), n_waveforms - n)
                m > 0 || continue
                wfs = col[1:m]
                A = flatview(wfs.signal)
                if X === nothing
                    X, dt = Matrix{Float32}(undef, size(A, 1), n_waveforms), step(wfs.time[1])
                end
                size(A, 1) == size(X, 1) || continue
                X[:, (n + 1):(n + m)] .= A
                n += m
            end
        end
    end
    n == n_waveforms || throw(ArgumentError("Files contain only $n waveforms"))
    return (waveforms = X, dt = dt)
end

# Mean of the large pulses, with the baseline subtracted, and the position of its half maximum:
function mean_pulse(X::AbstractMatrix)
    baseline = mean(view(X, 1:100, :), dims = 1)
    heights = vec(maximum(X, dims = 1) .- baseline)
    lo, hi = quantile(heights, (0.90, 0.99))
    large = findall(h -> lo < h < hi, heights)
    pulse = vec(mean(X[:, large] .- baseline[:, large], dims = 2))
    return pulse, findfirst(>(maximum(pulse) / 2), pulse)
end

# Baseline window and decay time (in samples) from the mean of the large pulses:
function chain_parameters(X::AbstractMatrix; rise = 78, flat = 23)
    pulse, i_rise = mean_pulse(X)
    i1, i2 = i_rise + 60, length(pulse) - 5
    tau = (i2 - i1) / log(pulse[i1] / pulse[i2])
    EnergyChain(i_rise - 60, InvCRFilter(Float32(tau)), TrapezoidalChargeFilter(rise, flat))
end

if abspath(PROGRAM_FILE) == @__FILE__
    data = read_waveforms(ARGS[1], length(ARGS) >= 2 ? ARGS[2] : nothing)
    chain = chain_parameters(data.waveforms)
    println("channel ", data.channel, ": ", size(data.waveforms, 2), " waveforms of ", size(data.waveforms, 1), " samples, sampling step ", data.dt)
    println("baseline window 1:", chain.n_baseline, ", decay time ", chain.pz.cr, " samples")
    results, _ = run_benchmarks(data.waveforms, chain; label = data.channel)
    write_csv(joinpath(@__DIR__, "energy_chain_lh5.csv"), results)
end
