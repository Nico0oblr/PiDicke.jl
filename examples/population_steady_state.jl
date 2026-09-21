using LinearAlgebra
using PIDicke

model = DickeModel(
    30;
    collective_decay = 1.0,
    local_pump = 0.4,
    local_decay = 0.02,
    local_dephasing = 0.05,
)

states, _, generator = population_generator(model)
steady = population_steady_state(generator)

mean_S = sum(S * steady[i] for (i, (S, _)) in enumerate(states))
mean_M = sum(M * steady[i] for (i, (_, M)) in enumerate(states))
mean_intensity = sum(A_JM_minus(S, M)^2 * steady[i] for (i, (S, M)) in enumerate(states))

println("normalization: ", sum(steady))
println("stationary residual: ", norm(generator * steady))
println("<S>: ", mean_S)
println("<M>: ", mean_M)
println("<S+S->: ", mean_intensity)
