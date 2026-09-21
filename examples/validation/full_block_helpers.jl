using LinearAlgebra
using PIDicke

"""Map every integer `(S,M)` label to its index in a Dicke `SumBasis`."""
function full_dicke_indices(space)
    indices = Dict{Tuple{Int,Int},Int}()
    for (sector_index, basis) in enumerate(space.bases)
        S = Int(basis.spinnumber)
        block = block_ranges(space)[sector_index]
        for (local_index, full_index) in enumerate(block)
            M = S - (local_index - 1)
            indices[(S, M)] = full_index
        end
    end
    return indices
end

"""Extract the diagonal-population block from a full PI Liouvillian."""
function population_block_from_full(liouvillian, space, states)
    full_indices = full_dicke_indices(space)
    dimension = length(space)
    block = zeros(ComplexF64, length(states), length(states))

    for (source_column, source) in enumerate(states)
        operator = zeros(ComplexF64, dimension, dimension)
        source_index = full_indices[source]
        operator[source_index, source_index] = 1
        derivative = reshape(liouvillian.data * vec(operator), dimension, dimension)

        for (target_row, target) in enumerate(states)
            target_index = full_indices[target]
            block[target_row, source_column] = derivative[target_index, target_index]
        end
    end
    return block
end

"""Extract the `|S,M-1><S,M|` block from a full PI Liouvillian."""
function coherence_block_from_full(liouvillian, space, states)
    full_indices = full_dicke_indices(space)
    dimension = length(space)
    block = zeros(ComplexF64, length(states), length(states))

    for (source_column, (S, M)) in enumerate(states)
        operator = zeros(ComplexF64, dimension, dimension)
        operator[full_indices[(S, M - 1)], full_indices[(S, M)]] = 1
        derivative = reshape(liouvillian.data * vec(operator), dimension, dimension)

        for (target_row, (target_S, target_M)) in enumerate(states)
            block[target_row, source_column] = derivative[
                full_indices[(target_S, target_M - 1)],
                full_indices[(target_S, target_M)],
            ]
        end
    end
    return block
end

function population_operator(population, states, space)
    full_indices = full_dicke_indices(space)
    operator = zeros(ComplexF64, length(space), length(space))
    for (weight, state) in zip(population, states)
        index = full_indices[state]
        operator[index, index] = weight
    end
    return operator
end

function population_from_operator(operator, states, space)
    full_indices = full_dicke_indices(space)
    return ComplexF64[operator[full_indices[state], full_indices[state]] for state in states]
end

function coherence_operator(coherences, states, space)
    full_indices = full_dicke_indices(space)
    operator = zeros(ComplexF64, length(space), length(space))
    for (weight, (S, M)) in zip(coherences, states)
        operator[full_indices[(S, M - 1)], full_indices[(S, M)]] = weight
    end
    return operator
end

function coherences_from_operator(operator, states, space)
    full_indices = full_dicke_indices(space)
    return ComplexF64[
        operator[full_indices[(S, M - 1)], full_indices[(S, M)]]
        for (S, M) in states
    ]
end
