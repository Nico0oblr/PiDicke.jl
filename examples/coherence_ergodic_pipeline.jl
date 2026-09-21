using PIDicke
using Random

# Full stochastic workflow from the original simulate_correlator_ergodic:
# one long population trajectory -> stationary (S,M) seeds -> coherence paths.
model = DickeModel(
    20;
    collective_decay = 1.0,
    local_pump = 0.4,
    local_dephasing = 0.02,
)

coherence_times = range(0.0, 10.0, length = 501)
result = ergodic_first_order_correlation(
    model,
    coherence_times;
    nsamples = 500,
    equilibration_time = 10.0,
    sampling_window = 30.0,
    rng = Xoshiro(1234),
)

println("number of stationary seeds: ", length(result.samples))
println("sampling interval: ", extrema(result.samples.sample_times))
println("C(0): ", first(result.correlation))
println("C(tmax): ", last(result.correlation))

# The two stages can also be called separately:
# samples = ergodic_population_samples(model, 500; ...)
# times, correlation = simulate_first_order_correlation_from_samples(
#     model, samples, coherence_times,
# )
