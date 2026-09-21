using OrdinaryDiffEq
using PIDicke

model = DickeModel(
    20;
    collective_decay = 1.0,
    local_pump = 0.2,
    local_dephasing = 0.05,
)

states, indices, generator = population_generator(model)
initial = population_initial_state(indices, (10, 10))
times = range(0.0, 5.0, length = 201)
result = solve_populations(model, initial, times)

excitation = [sum((M + model.N / 2) * p[i] for (i, (_, M)) in enumerate(states)) for p in result.solution.u]
println("final excitation = ", last(excitation))
