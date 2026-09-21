using PIDicke
using Random

# Ensemble version of the trajectory workflow used throughout the original
# MonteCarloMF and PopulationMC calculations.
model = DickeModel(
    20;
    collective_decay = 1.0,
    local_pump = 0.25,
    local_dephasing = 0.05,
)

excitation(S, M, model) = M + model.N / 2

times, average_excitation = simulate_population_ensemble(
    model,
    excitation,
    250;
    tmax = 3.0,
    ngrid = 151,
    rng = Xoshiro(1234),
)

println("initial excitation: ", first(average_excitation))
println("final ensemble-averaged excitation: ", last(average_excitation))
println("number of output times: ", length(times))
