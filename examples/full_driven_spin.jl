using PIDicke
using QuantumOptics: expect

model = DickeModel(
    6;
    collective_decay = 1.0,
    local_decay = 0.05,
    local_dephasing = 0.02,
)

space = dicke_space(model)
Sx = collective_operator(space, :x)
Sz = collective_operator(space, :z)
hamiltonian = 0.2 * Sx

times = range(0.0, 3.0, length = 151)
result = master_evolution(
    times,
    fully_ground_state(model),
    model;
    hamiltonian,
)

inversion = real.(expect(Sz, result.states))
println("initial <Sz>: ", first(inversion))
println("final <Sz>: ", last(inversion))
