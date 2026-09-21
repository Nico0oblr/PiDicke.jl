using PIDicke
using Random

const Active = PIDicke.Experimental

model = DickeModel(
    20;
    collective_decay = 1.0,
    local_decay = 0.1,
    local_pump = 0.4,
    local_dephasing = 0.05,
)

samples = ergodic_population_samples(
    model,
    2_000;
    equilibration_time = 10.0,
    sampling_window = 30.0,
    rng = Xoshiro(1234),
)

result = Active.active_region_steady_state(
    model,
    samples;
    padding = 2,
    dt = 1e3,
    nsteps = 5,
)

println("active basis size: ", length(result.basis))
println("normalization: ", result.diagnostics.normalization)
println("stationary residual (L1): ", result.diagnostics.residual_l1)
println("mean omitted outgoing rate: ", result.diagnostics.mean_leakage)
