using LinearAlgebra
using PIDicke

include("full_block_helpers.jl")

# Rates are normalized by Gamma0; tau = Gamma0*t.
Gamma0 = 1.0
model = DickeModel(
    4;
    collective_decay = 0.7 * Gamma0,
    collective_pump = 0.11 * Gamma0,
    local_decay = 0.13 * Gamma0,
    local_pump = 0.19 * Gamma0,
    local_dephasing = 0.17 * Gamma0,
)

coherence_states, coherence_indices, coherence_generator_matrix =
    first_order_coherence_generator(model)
full = full_liouvillian(model)
full_coherence_block = coherence_block_from_full(
    full.liouvillian,
    full.space,
    coherence_states,
)

generator_error = maximum(abs, full_coherence_block - coherence_generator_matrix)

population_states, _, _ = population_generator(model)
steady_population = population_steady_state(model)
initial = coherence_initial_state(
    population_states,
    steady_population,
    coherence_indices,
)
readout = coherence_readout_vector(coherence_states)

tau = 0.8
coherence_at_tau = exp(Matrix(coherence_generator_matrix) * tau / Gamma0) * initial
initial_operator = coherence_operator(initial, coherence_states, full.space)
full_operator_at_tau = reshape(
    exp(Matrix(full.liouvillian.data) * tau / Gamma0) * vec(initial_operator),
    length(full.space),
    length(full.space),
)
full_coherence_at_tau = coherences_from_operator(
    full_operator_at_tau,
    coherence_states,
    full.space,
)

coherence_error = maximum(abs, full_coherence_at_tau - coherence_at_tau)
correlation_error = abs(
    dot(readout, full_coherence_at_tau) - dot(readout, coherence_at_tau),
)

println("maximum generator-entry error: ", generator_error)
println("maximum coherence error at tau = ", tau, ": ", coherence_error)
println("correlation error at tau = ", tau, ": ", correlation_error)
