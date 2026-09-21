# Convention: `(S,M)` represents the operator |S,M-1><S,M|.

@inline _valid_coherence_state(S::Int, M::Int, Jmax::Int) =
    0 <= S <= Jmax && abs(M) <= S && abs(M - 1) <= S

function first_order_coherence_basis(model::DickeModel)
    _require_even_N(model, FirstOrderCoherenceRepresentation())
    Jmax = model.N ÷ 2
    return [(S, M) for S in 0:Jmax for M in -S:S if _valid_coherence_state(S, M, Jmax)]
end

@inline _rate_collective_decay(S, M, model) =
    model.collective_decay * A_JM_minus2(S, M)
@inline _rate_collective_pump(S, M, model) =
    model.collective_pump * A_JM_plus2(S, M)
@inline _rate_local_decay_plus(S, M, model) =
    model.local_decay * P_JM_minus_plus2(S, M, model.N)
@inline _rate_local_decay_zero(S, M, model) =
    model.local_decay * P_JM_minus_02(S, M, model.N)
@inline _rate_local_decay_minus(S, M, model) =
    model.local_decay * P_JM_minus_minus2(S, M, model.N)
@inline _rate_local_pump_plus(S, M, model) =
    model.local_pump * P_JM_plus_plus2(S, M, model.N)
@inline _rate_local_pump_zero(S, M, model) =
    model.local_pump * P_JM_plus_02(S, M, model.N)
@inline _rate_local_pump_minus(S, M, model) =
    model.local_pump * P_JM_plus_minus2(S, M, model.N)
@inline _rate_dephasing_plus(S, M, model) =
    model.local_dephasing * P_JM_z_plus2(S, M, model.N)
@inline _rate_dephasing_zero(S, M, model) =
    model.local_dephasing * P_JM_z_02(S, M, model.N)
@inline _rate_dephasing_minus(S, M, model) =
    model.local_dephasing * P_JM_z_minus2(S, M, model.N)

@inline function _coherence_internal_and_loss(rate_M::Real, rate_Mminus1::Real)
    internal = sqrt(rate_M * rate_Mminus1)
    loss = (rate_M + rate_Mminus1) / 2
    return internal, loss
end

function _coherence_branches(model::DickeModel, S::Int, M::Int)
    branches = Tuple{Int,Int,Float64,Float64}[]
    Jmax = model.N ÷ 2

    function add!(dS, dM, rate_function)
        internal, loss = _coherence_internal_and_loss(
            rate_function(S, M, model),
            rate_function(S, M - 1, model),
        )
        target_valid = _valid_coherence_state(S + dS, M + dM, Jmax)
        push!(branches, (dS, dM, target_valid ? internal : 0.0, loss))
        return nothing
    end

    model.collective_decay > 0 && add!(0, -1, _rate_collective_decay)
    model.collective_pump > 0 && add!(0, 1, _rate_collective_pump)
    if model.local_decay > 0
        add!(1, -1, _rate_local_decay_plus)
        add!(0, -1, _rate_local_decay_zero)
        add!(-1, -1, _rate_local_decay_minus)
    end
    if model.local_pump > 0
        add!(1, 1, _rate_local_pump_plus)
        add!(0, 1, _rate_local_pump_zero)
        add!(-1, 1, _rate_local_pump_minus)
    end
    if model.local_dephasing > 0
        add!(1, 0, _rate_dephasing_plus)
        add!(0, 0, _rate_dephasing_zero)
        add!(-1, 0, _rate_dephasing_minus)
    end
    return branches
end

"""Construct the sparse generator for `|S,M-1><S,M|` coherences."""
function first_order_coherence_generator(model::DickeModel)
    states = first_order_coherence_basis(model)
    indices = Dict(state => i for (i, state) in enumerate(states))
    entries = DefaultDict{Tuple{Int,Int},ComplexF64}(0.0)

    for (source_index, (S, M)) in enumerate(states)
        total_loss = 0.0
        for (dS, dM, internal, loss) in _coherence_branches(model, S, M)
            total_loss += loss
            if internal > 0
                target = (S + dS, M + dM)
                entries[(indices[target], source_index)] += internal
            end
        end
        entries[(source_index, source_index)] -= total_loss
    end

    rows, columns, values = Int[], Int[], ComplexF64[]
    for ((row, column), value) in entries
        value == 0 && continue
        push!(rows, row)
        push!(columns, column)
        push!(values, value)
    end
    generator = sparse(rows, columns, values, length(states), length(states))
    return states, indices, generator
end

function coherence_initial_state(
    population_states,
    populations::AbstractVector,
    coherence_indices::AbstractDict,
)
    length(population_states) == length(populations) || throw(DimensionMismatch("population vector has the wrong length"))
    initial = zeros(ComplexF64, length(coherence_indices))
    for (index, (S, M)) in enumerate(population_states)
        if haskey(coherence_indices, (S, M))
            initial[coherence_indices[(S, M)]] = A_JM_minus(S, M) * populations[index]
        end
    end
    return initial
end

coherence_readout_vector(states) = ComplexF64[A_JM_minus(S, M) for (S, M) in states]

function solve_first_order_coherences(
    model::DickeModel,
    initial::AbstractVector,
    times;
    solver = Tsit5(),
    kwargs...,
)
    states, indices, generator = first_order_coherence_generator(model)
    length(initial) == length(states) || throw(DimensionMismatch("initial coherence vector has the wrong length"))
    rhs! = (du, u, _, _) -> mul!(du, generator, u)
    problem = ODEProblem(rhs!, ComplexF64.(initial), (first(times), last(times)))
    solution = solve(problem, solver; saveat = times, kwargs...)
    return (; states, indices, generator, solution)
end

"""Compute the first-order correlation from a Dicke population distribution."""
function first_order_correlation(
    model::DickeModel,
    populations::AbstractVector,
    times;
    solver = Tsit5(),
    kwargs...,
)
    population_states, _, _ = population_generator(model)
    coherence_states, coherence_indices, _ = first_order_coherence_generator(model)
    initial = coherence_initial_state(population_states, populations, coherence_indices)
    result = solve_first_order_coherences(model, initial, times; solver, kwargs...)
    readout = coherence_readout_vector(coherence_states)
    correlation = [dot(readout, state) for state in result.solution.u]
    return merge(result, (; correlation, initial))
end

"""
    simulate_first_order_coherence_trajectory(model, S, M, times; ...)

Stochastic unraveling of the reduced first-order coherence generator. The
initial label `(S,M)` represents `|S,M-1><S,M|`. This is the useful core of the
original `CoherenceMC.jl`; stationary burn-in sampling is intentionally left
out.
"""
function simulate_first_order_coherence_trajectory(
    model::DickeModel,
    S0::Int,
    M0::Int,
    times;
    initial_weight::Real = A_JM_minus(S0, M0),
    rng::AbstractRNG = Random.default_rng(),
)
    _require_even_N(model, FirstOrderCoherenceRepresentation())
    grid = collect(times)
    contribution = zeros(Float64, length(grid))
    _valid_coherence_state(S0, M0, model.N ÷ 2) || return contribution

    S, M, time = S0, M0, 0.0
    absorbed = false
    index = 1
    while index <= length(grid)
        value = absorbed ? 0.0 : initial_weight * A_JM_minus(S, M)
        branches = _coherence_branches(model, S, M)
        total_loss = sum(branch[4] for branch in branches)
        if total_loss <= 0
            contribution[index:end] .= value
            break
        end

        time += -log(rand(rng)) / total_loss
        while index <= length(grid) && grid[index] <= time
            contribution[index] = value
            index += 1
        end
        index > length(grid) && break

        threshold = rand(rng) * total_loss
        accumulated = 0.0
        jumped = false
        for (dS, dM, internal, _) in branches
            accumulated += internal
            if threshold <= accumulated
                S += dS
                M += dM
                jumped = true
                break
            end
        end
        if !jumped
            absorbed = true
            contribution[index:end] .= 0.0
            break
        end
    end
    return contribution
end

"""Average reduced-coherence trajectories initialized from supplied `(S,M)` samples."""
function simulate_first_order_correlation_from_samples(
    model::DickeModel,
    S_samples::AbstractVector{Int},
    M_samples::AbstractVector{Int},
    times;
    rng::AbstractRNG = Random.default_rng(),
)
    length(S_samples) == length(M_samples) || throw(DimensionMismatch("sample vectors must have equal length"))
    isempty(S_samples) && throw(ArgumentError("at least one sample is required"))
    correlation = zeros(Float64, length(times))
    for (S, M) in zip(S_samples, M_samples)
        correlation .+= simulate_first_order_coherence_trajectory(model, S, M, times; rng)
    end
    correlation ./= length(S_samples)
    return collect(times), correlation
end
