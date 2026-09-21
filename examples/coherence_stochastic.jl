using PIDicke
using Random

model = DickeModel(
    20;
    collective_decay = 1.0,
    local_pump = 0.4,
    local_dephasing = 0.02,
)

population_states, _, _ = population_generator(model)
steady = population_steady_state(model)

# Draw starting Dicke states directly from the exact stationary distribution;
# no burn-in trajectory is needed.
rng = Xoshiro(1234)
cumulative = cumsum(steady)
sample_indices = [searchsortedfirst(cumulative, rand(rng)) for _ in 1:500]
S_samples = [population_states[index][1] for index in sample_indices]
M_samples = [population_states[index][2] for index in sample_indices]

times = range(0.0, 10.0, length = 501)
_, stochastic_correlation = simulate_first_order_correlation_from_samples(
    model,
    S_samples,
    M_samples,
    times;
    rng,
)
deterministic = first_order_correlation(model, steady, times)

println("deterministic C(0): ", real(first(deterministic.correlation)))
println("stochastic C(0): ", first(stochastic_correlation))
println("deterministic C(tmax): ", real(last(deterministic.correlation)))
println("stochastic C(tmax): ", last(stochastic_correlation))
