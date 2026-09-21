using LinearAlgebra
using PIDicke
using Random

# Rates are normalized by Gamma0; tau = Gamma0*t.
Gamma0 = 1.0
model = DickeModel(
    20;
    collective_decay = 0.3 * Gamma0,
    local_decay = 0.05 * Gamma0,
    local_pump = 0.4 * Gamma0,
    local_dephasing = 0.02 * Gamma0,
)

states, indices, _ = population_generator(model)
exact = population_steady_state(model)
rng = Xoshiro(1234)

println("nsamples    L1 error    largest-bin error")
for nsamples in (1_000, 10_000, 100_000)
    # Increase the observed trajectory duration together with the number of
    # sample times. Densifying samples in one fixed window would not add new
    # ergodic information once that path is well resolved.
    sampling_window = nsamples / (10 * Gamma0)
    samples = ergodic_population_samples(
        model,
        nsamples;
        equilibration_time = 100 / Gamma0,
        sampling_window,
        rng,
    )
    histogram = zeros(length(states))
    for state in zip(samples.S, samples.M)
        histogram[indices[state]] += 1
    end
    histogram ./= sum(histogram)

    l1_error = norm(histogram - exact, 1)
    maximum_error = maximum(abs, histogram - exact)
    println(
        rpad(string(nsamples), 12),
        rpad(string(round(l1_error; sigdigits = 6)), 12),
        round(maximum_error; sigdigits = 6),
    )
end

println("Samples come from one trajectory and are correlated; errors are convergence diagnostics, not independent-sample confidence intervals.")
