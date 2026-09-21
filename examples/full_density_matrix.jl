using LinearAlgebra
using PIDicke

model = DickeModel(
    6;
    collective_decay = 1.0,
    local_decay = 0.1,
    local_dephasing = 0.05,
)

full = full_liouvillian(model)
println("Dicke Hilbert-space dimension: ", length(full.space))
println("Liouvillian size: ", size(full.liouvillian.data))

times = range(0.0, 2.0, length = 101)
result = master_evolution(times, fully_excited_state(model), model)
println("final trace: ", tr(last(result.states)))
