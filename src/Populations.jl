"""
    population_generator(model)

Construct the sparse Markov generator for diagonal Dicke populations. It uses
the same tested event builder as the Gillespie trajectories.
"""
function population_generator(model::DickeModel)
    _require_even_N(model, PopulationRepresentation())
    states = dicke_triangle(model)
    indices = Dict(state => i for (i, state) in enumerate(states))
    rows, columns, values = Int[], Int[], Float64[]
    buffer = PopulationEventBuffer()

    for (source_index, (S, M)) in enumerate(states)
        state = PopulationState(S, M, 0.0)
        _build_population_events!(buffer, model, state)
        total_rate = 0.0
        for k in 1:buffer.n
            target = (S + buffer.dS[k], M + buffer.dM[k])
            rate = buffer.rates[k]
            push!(rows, indices[target])
            push!(columns, source_index)
            push!(values, rate)
            total_rate += rate
        end
        if total_rate > 0
            push!(rows, source_index)
            push!(columns, source_index)
            push!(values, -total_rate)
        end
    end

    generator = sparse(rows, columns, values, length(states), length(states))
    dropzeros!(generator)
    return states, indices, generator
end

function population_initial_state(indices::AbstractDict, state::Tuple)
    haskey(indices, state) || throw(ArgumentError("state $state is not in the basis"))
    initial = zeros(Float64, length(indices))
    initial[indices[state]] = 1.0
    return initial
end

function population_initial_state(model::DickeModel, state::Tuple = (model.N ÷ 2, model.N ÷ 2))
    _, indices, _ = population_generator(model)
    return population_initial_state(indices, state)
end

"""Solve `G*p=0` with unit normalization."""
function population_steady_state(generator::AbstractMatrix)
    matrix = Matrix(generator)
    matrix[end, :] .= 1
    rhs = zeros(eltype(matrix), size(matrix, 1))
    rhs[end] = 1
    return matrix \ rhs
end

population_steady_state(model::DickeModel) = population_steady_state(last(population_generator(model)))

function solve_populations(
    model::DickeModel,
    initial::AbstractVector,
    times;
    solver = Tsit5(),
    kwargs...,
)
    states, indices, generator = population_generator(model)
    length(initial) == length(states) || throw(DimensionMismatch("initial population has the wrong length"))
    rhs! = (du, u, _, _) -> mul!(du, generator, u)
    problem = ODEProblem(rhs!, collect(initial), (first(times), last(times)))
    solution = solve(problem, solver; saveat = times, kwargs...)
    return (; states, indices, generator, solution)
end
