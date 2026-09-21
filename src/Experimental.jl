"""
Experimental algorithms whose numerical behavior is useful but whose public
interface is not yet stable.

Access these functions through `PIDicke.Experimental`; they are deliberately
not exported by the top-level `PIDicke` module.
"""
module Experimental

using LinearAlgebra
using SparseArrays

using ..PIDicke: DickeModel,
                 ErgodicDickeSamples,
                 PopulationEventBuffer,
                 PopulationState,
                 _build_population_events!,
                 _require_even_N,
                 PopulationRepresentation

export padded_sample_basis,
       active_population_generator,
       population_from_samples,
       steady_by_implicit_euler,
       stationary_diagnostics,
       active_region_steady_state

"""
    padded_sample_basis(model, samples; padding)
    padded_sample_basis(model, S_samples, M_samples; padding)

Construct the active population basis used in the working notebook
implementation. The sample extrema define a rectangle in `(S,M)`, the same
integer `padding` is added on every side, and the result is intersected with
the physical Dicke triangle `0 <= S <= N/2`, `|M| <= S`.

This is intentionally an extrema-based rectangle, not an outgoing-transition
shell construction.
"""
function padded_sample_basis(
    model::DickeModel,
    S_samples,
    M_samples;
    padding::Integer,
)
    _require_even_N(model, PopulationRepresentation())
    length(S_samples) == length(M_samples) ||
        throw(DimensionMismatch("S_samples and M_samples must have equal length"))
    isempty(S_samples) && throw(ArgumentError("at least one sample is required"))
    padding >= 0 || throw(ArgumentError("padding must be nonnegative"))

    sampled_S = round.(Int, S_samples)
    sampled_M = round.(Int, M_samples)
    Jmax = model.N ÷ 2
    all(eachindex(sampled_S)) do index
        S, M = sampled_S[index], sampled_M[index]
        0 <= S <= Jmax && abs(M) <= S
    end || throw(ArgumentError("samples must lie inside the Dicke triangle"))

    Smin = clamp(minimum(sampled_S) - padding, 0, Jmax)
    Smax = clamp(maximum(sampled_S) + padding, 0, Jmax)
    Mmin = minimum(sampled_M) - padding
    Mmax = maximum(sampled_M) + padding

    return [
        (S, M)
        for S in Smin:Smax
        for M in max(-S, Mmin):min(S, Mmax)
    ]
end

padded_sample_basis(model::DickeModel, samples::ErgodicDickeSamples; padding::Integer) =
    padded_sample_basis(model, samples.S, samples.M; padding)

"""
    active_population_generator(model, basis; close_boundary=true)

Build `dp/dt = G*p` on a supplied subset of the Dicke triangle. Columns are
source states and rows are destination states. `leakage[j]` is the total rate
from source `j` to states omitted from `basis`.

With `close_boundary=true`, omitted transitions are removed from both the
off-diagonal and diagonal entries, so the restricted generator conserves
probability. With `false`, the diagonal retains the omitted rates and the
restricted dynamics loses probability at the measured leakage rate.
"""
function active_population_generator(
    model::DickeModel,
    active_basis;
    close_boundary::Bool = true,
)
    _require_even_N(model, PopulationRepresentation())
    basis = Tuple{Int,Int}[(Int(S), Int(M)) for (S, M) in active_basis]
    isempty(basis) && throw(ArgumentError("active basis must not be empty"))
    length(unique(basis)) == length(basis) ||
        throw(ArgumentError("active basis contains duplicate states"))

    Jmax = model.N ÷ 2
    all(state -> 0 <= state[1] <= Jmax && abs(state[2]) <= state[1], basis) ||
        throw(ArgumentError("active basis contains a state outside the Dicke triangle"))

    indices = Dict(state => index for (index, state) in enumerate(basis))
    rows, columns, values = Int[], Int[], Float64[]
    leakage = zeros(Float64, length(basis))
    buffer = PopulationEventBuffer()

    for (source_index, (S, M)) in enumerate(basis)
        state = PopulationState(S, M, 0.0)
        _build_population_events!(buffer, model, state)
        internal_rate = 0.0
        total_rate = 0.0

        for event_index in 1:buffer.n
            rate = buffer.rates[event_index]
            target = (S + buffer.dS[event_index], M + buffer.dM[event_index])
            total_rate += rate
            if haskey(indices, target)
                push!(rows, indices[target])
                push!(columns, source_index)
                push!(values, rate)
                internal_rate += rate
            else
                leakage[source_index] += rate
            end
        end

        diagonal_rate = close_boundary ? internal_rate : total_rate
        if diagonal_rate > 0
            push!(rows, source_index)
            push!(columns, source_index)
            push!(values, -diagonal_rate)
        end
    end

    generator = sparse(rows, columns, values, length(basis), length(basis))
    dropzeros!(generator)
    return (; basis, indices, generator, leakage, close_boundary)
end

"""Create a normalized population histogram on an active basis from samples."""
function population_from_samples(basis, indices::AbstractDict, S_samples, M_samples)
    length(S_samples) == length(M_samples) ||
        throw(DimensionMismatch("S_samples and M_samples must have equal length"))
    population = zeros(Float64, length(basis))

    for (S, M) in zip(S_samples, M_samples)
        state = (round(Int, S), round(Int, M))
        haskey(indices, state) && (population[indices[state]] += 1)
    end

    normalization = sum(population)
    normalization > 0 || throw(ArgumentError("no samples lie inside the active basis"))
    population ./= normalization
    return population
end

population_from_samples(basis, indices::AbstractDict, samples::ErgodicDickeSamples) =
    population_from_samples(basis, indices, samples.S, samples.M)

"""
    steady_by_implicit_euler(generator, initial; dt=1e4, nsteps=5, ...)

Approach the stationary population by repeated pseudo-time steps
`p <- (I - dt*G) \\ p`. The factorization is reused for every step. This is the
steady-state iteration used by the active-space notebooks.
"""
function steady_by_implicit_euler(
    generator::AbstractMatrix,
    initial::AbstractVector;
    dt::Real = 1e4,
    nsteps::Integer = 5,
    renormalize::Bool = true,
    clip_negative::Bool = true,
)
    size(generator, 1) == size(generator, 2) == length(initial) ||
        throw(DimensionMismatch("generator and initial population sizes do not match"))
    dt > 0 || throw(ArgumentError("dt must be positive"))
    nsteps > 0 || throw(ArgumentError("nsteps must be positive"))

    n = length(initial)
    identity_matrix = spdiagm(0 => ones(Float64, n))
    factorization = factorize(identity_matrix - Float64(dt) * generator)
    population = collect(Float64, initial)

    for _ in 1:nsteps
        population = factorization \ population
        clip_negative && (population[population .< 0] .= 0)
        if renormalize
            normalization = sum(population)
            normalization > 0 || error("implicit Euler iteration lost all population")
            population ./= normalization
        end
    end
    return population
end

"""Return conservation, stationarity, and active-boundary leakage diagnostics."""
function stationary_diagnostics(generator, population; leakage = nothing)
    diagnostics = (
        normalization = sum(population),
        minimum = minimum(population),
        maximum = maximum(population),
        residual_l1 = norm(generator * population, 1),
        residual_linf = norm(generator * population, Inf),
        column_sum_linf = norm(vec(sum(generator; dims = 1)), Inf),
    )
    leakage === nothing && return diagnostics
    return merge(diagnostics, (;
        mean_leakage = dot(leakage, population),
        maximum_leakage = maximum(leakage),
    ))
end

"""
    active_region_steady_state(model, samples; padding, dt=1e3, nsteps=5)

Run the notebook active-region pipeline: construct the padded sample rectangle,
close its population boundary, initialize it with the sample histogram, and
apply the implicit-Euler stationary iteration.
"""
function active_region_steady_state(
    model::DickeModel,
    samples::ErgodicDickeSamples;
    padding::Integer,
    dt::Real = 1e3,
    nsteps::Integer = 5,
    renormalize::Bool = true,
    clip_negative::Bool = true,
)
    basis = padded_sample_basis(model, samples; padding)
    restricted = active_population_generator(model, basis; close_boundary = true)
    initial = population_from_samples(restricted.basis, restricted.indices, samples)
    steady = steady_by_implicit_euler(
        restricted.generator,
        initial;
        dt,
        nsteps,
        renormalize,
        clip_negative,
    )
    diagnostics = stationary_diagnostics(
        restricted.generator,
        steady;
        leakage = restricted.leakage,
    )
    return merge(restricted, (; initial, steady, diagnostics, samples))
end

end
