using LinearAlgebra
using PIDicke

include("full_block_helpers.jl")

# All rates are measured relative to an arbitrary reference rate Gamma0.
# Consequently tau = Gamma0*t is dimensionless time.
Gamma0 = 1.0
model = DickeModel(
    4;
    collective_decay = 0.7 * Gamma0,
    collective_pump = 0.11 * Gamma0,
    local_decay = 0.13 * Gamma0,
    local_pump = 0.19 * Gamma0,
    local_dephasing = 0.17 * Gamma0,
)

states, indices, population_generator_matrix = population_generator(model)
full = full_liouvillian(model)
full_population_block = population_block_from_full(
    full.liouvillian,
    full.space,
    states,
)

generator_error = maximum(abs, full_population_block - population_generator_matrix)

initial = population_initial_state(indices, (2, 0))
tau = 0.8
population_at_tau = exp(Matrix(population_generator_matrix) * tau / Gamma0) * initial

initial_operator = population_operator(initial, states, full.space)
full_operator_at_tau = reshape(
    exp(Matrix(full.liouvillian.data) * tau / Gamma0) * vec(initial_operator),
    length(full.space),
    length(full.space),
)
full_population_at_tau = population_from_operator(full_operator_at_tau, states, full.space)
evolution_error = maximum(abs, full_population_at_tau - population_at_tau)

println("maximum generator-entry error: ", generator_error)
println("maximum population error at tau = ", tau, ": ", evolution_error)
