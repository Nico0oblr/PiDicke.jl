using PIDicke
using Random

# Both calculations use the same event/rate implementation. This example is a
# lightweight version of the original MC-versus-ODE comparison scripts.
model = DickeModel(
    12;
    collective_decay = 1.0,
    local_pump = 0.15,
    local_dephasing = 0.03,
)

times = range(0.0, 2.0, length = 101)
states, indices, _ = population_generator(model)
initial = population_initial_state(indices, (6, 6))
deterministic = solve_populations(model, initial, times)

excitation(S, M, model) = M + model.N / 2
ode_excitation = [
    sum(excitation(S, M, model) * probabilities[i] for (i, (S, M)) in enumerate(states))
    for probabilities in deterministic.solution.u
]

_, mc_excitation = simulate_population_ensemble(
    model,
    excitation,
    500;
    tmax = last(times),
    ngrid = length(times),
    rng = Xoshiro(1234),
)

println("final ODE excitation: ", last(ode_excitation))
println("final Monte Carlo excitation: ", last(mc_excitation))
println("grid RMS difference: ", sqrt(sum(abs2, ode_excitation .- mc_excitation) / length(times)))
