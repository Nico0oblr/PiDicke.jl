using LinearAlgebra
using PIDicke

# The efficient integer Dicke triangle is currently even-N only, but the full
# QuantumOptics SumBasis naturally handles half-integer spin sectors.
model = DickeModel(5; collective_decay = 1.0, local_decay = 0.05)
space = dicke_space(model)
full = full_liouvillian(model)

times = range(0.0, 1.0, length = 51)
result = master_evolution(times, fully_excited_state(model), model)

println("spin sectors: ", getproperty.(space.bases, :spinnumber))
println("Hilbert-space dimension: ", length(space))
println("Liouvillian size: ", size(full.liouvillian.data))
println("final trace: ", tr(last(result.states)))
