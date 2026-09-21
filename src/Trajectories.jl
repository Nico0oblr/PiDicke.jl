# This is the tested Gillespie implementation from PopulationMC.jl, with the
# model type unified and the burn-in/stationary-sampling helpers intentionally
# omitted.

mutable struct PopulationState
    S::Int
    M::Int
    t::Float64
end

mutable struct PopulationEventBuffer
    rates::Vector{Float64}
    dS::Vector{Int}
    dM::Vector{Int}
    n::Int
end

PopulationEventBuffer(max_events::Int = 11) =
    PopulationEventBuffer(zeros(max_events), zeros(Int, max_events), zeros(Int, max_events), 0)

@inline _valid_population_state(S::Int, M::Int, Jmax::Int) =
    0 <= S <= Jmax && abs(M) <= S

@inline function _push_population_event!(buf::PopulationEventBuffer, rate::Real, dS::Int, dM::Int)
    rate <= 0 && return nothing
    buf.n += 1
    buf.rates[buf.n] = rate
    buf.dS[buf.n] = dS
    buf.dM[buf.n] = dM
    return nothing
end

function _build_population_events!(buf::PopulationEventBuffer, model::DickeModel, state::PopulationState)
    _require_even_N(model, TrajectoryRepresentation())
    buf.n = 0
    S, M = state.S, state.M
    N, Jmax = model.N, model.N ÷ 2

    if model.collective_decay > 0 && _valid_population_state(S, M - 1, Jmax)
        _push_population_event!(buf, model.collective_decay * A_JM_minus2(S, M), 0, -1)
    end
    if model.collective_pump > 0 && _valid_population_state(S, M + 1, Jmax)
        _push_population_event!(buf, model.collective_pump * A_JM_plus2(S, M), 0, 1)
    end

    if model.local_decay > 0
        _valid_population_state(S + 1, M - 1, Jmax) &&
            _push_population_event!(buf, model.local_decay * P_JM_minus_plus2(S, M, N), 1, -1)
        _valid_population_state(S, M - 1, Jmax) &&
            _push_population_event!(buf, model.local_decay * P_JM_minus_02(S, M, N), 0, -1)
        _valid_population_state(S - 1, M - 1, Jmax) &&
            _push_population_event!(buf, model.local_decay * P_JM_minus_minus2(S, M, N), -1, -1)
    end

    if model.local_pump > 0
        _valid_population_state(S + 1, M + 1, Jmax) &&
            _push_population_event!(buf, model.local_pump * P_JM_plus_plus2(S, M, N), 1, 1)
        _valid_population_state(S, M + 1, Jmax) &&
            _push_population_event!(buf, model.local_pump * P_JM_plus_02(S, M, N), 0, 1)
        _valid_population_state(S - 1, M + 1, Jmax) &&
            _push_population_event!(buf, model.local_pump * P_JM_plus_minus2(S, M, N), -1, 1)
    end

    if model.local_dephasing > 0
        _valid_population_state(S + 1, M, Jmax) &&
            _push_population_event!(buf, model.local_dephasing * P_JM_z_plus2(S, M, N), 1, 0)
        _valid_population_state(S - 1, M, Jmax) &&
            _push_population_event!(buf, model.local_dephasing * P_JM_z_minus2(S, M, N), -1, 0)
    end
    return nothing
end

function _gillespie_step!(
    rng::AbstractRNG,
    state::PopulationState,
    buf::PopulationEventBuffer,
    model::DickeModel,
)
    _build_population_events!(buf, model, state)
    buf.n == 0 && return false
    total_rate = sum(@view buf.rates[1:buf.n])
    total_rate <= 0 && return false

    state.t += -log(rand(rng)) / total_rate
    threshold = rand(rng) * total_rate
    accumulated = 0.0
    selected = buf.n
    @inbounds for k in 1:buf.n
        accumulated += buf.rates[k]
        if threshold <= accumulated
            selected = k
            break
        end
    end
    state.S += buf.dS[selected]
    state.M += buf.dM[selected]
    return true
end

"""Simulate one Gillespie trajectory on the Dicke triangle."""
function run_population_trajectory(
    model::DickeModel,
    observable = (S, M, _) -> (S, M);
    S0::Int = model.N ÷ 2,
    M0::Int = model.N ÷ 2,
    tmax::Real = Inf,
    max_steps::Int = typemax(Int),
    stop_condition = nothing,
    rng::AbstractRNG = Random.default_rng(),
    sizehint_steps::Int = 1024,
)
    _require_even_N(model, TrajectoryRepresentation())
    _valid_population_state(S0, M0, model.N ÷ 2) || throw(ArgumentError("invalid initial state"))
    state = PopulationState(S0, M0, 0.0)
    buf = PopulationEventBuffer()
    times = Float64[]
    values = Vector{typeof(observable(S0, M0, model))}()
    sizehint!(times, sizehint_steps)
    sizehint!(values, sizehint_steps)
    push!(times, 0.0)
    push!(values, observable(S0, M0, model))

    steps = 0
    while state.t < tmax && steps < max_steps
        if stop_condition !== nothing && stop_condition(state.S, state.M, state.t, model)
            break
        end
        _gillespie_step!(rng, state, buf, model) || break
        push!(times, state.t)
        push!(values, observable(state.S, state.M, model))
        steps += 1
    end
    return times, values
end

"""Average an observable over independent population trajectories on a fixed grid."""
function simulate_population_ensemble(
    model::DickeModel,
    observable,
    ntrajectories::Int;
    S0::Int = model.N ÷ 2,
    M0::Int = model.N ÷ 2,
    tmax::Real,
    ngrid::Int = 1000,
    rng::AbstractRNG = Random.default_rng(),
)
    grid = collect(range(0.0, Float64(tmax), length = ngrid))
    average = zeros(Float64, ngrid)
    for _ in 1:ntrajectories
        times, values = run_population_trajectory(
            model,
            observable;
            S0,
            M0,
            tmax,
            rng,
        )
        index = 1
        @inbounds for k in eachindex(grid)
            while index < length(times) && times[index + 1] <= grid[k]
                index += 1
            end
            average[k] += values[index]
        end
    end
    average ./= ntrajectories
    return grid, average
end
