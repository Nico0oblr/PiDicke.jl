using PIDicke
using Random

# Port of the basic population trajectory in PopulationMC.jl. With collective
# decay only, the state remains on the S=N/2 edge of the Dicke triangle.
model = DickeModel(40; collective_decay = 1.0)

observable(S, M, _) = (
    S = S,
    M = M,
    excitation = M + model.N / 2,
    intensity = A_JM_minus(S, M)^2,
)

times, values = run_population_trajectory(
    model,
    observable;
    tmax = 1.0,
    stop_condition = (S, M, _, _) -> M == -S,
    rng = Xoshiro(1234),
)

intensities = getproperty.(values, :intensity)
peak_intensity, peak_index = findmax(intensities)

println("number of jumps: ", length(times) - 1)
println("peak intensity: ", peak_intensity, " at t = ", times[peak_index])
println("final Dicke state: ", (last(values).S, last(values).M))
