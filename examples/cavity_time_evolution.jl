using PIDicke
using QuantumOptics: expect, timeevolution

# Small executable version of the cavity construction in
# LasingQuantumOpticsImplementation.ipynb.
model = DickeModel(
    4;
    local_decay = 0.02,
    local_pump = 0.3,
    local_dephasing = 0.01,
)

system = cavity_system(
    model;
    cutoff = 3,
    g = 0.15,
    kappa = 1.0,
    cavity_drive = 0.02,
)

times = range(0.0, 5.0, length = 101)
tout, states = timeevolution.master(
    times,
    system.initial,
    system.hamiltonian,
    system.jumps;
    rates = system.rates,
)

photon_number = real.(expect(system.adag * system.a, states))
inversion = real.(expect(system.Sz, states))

println("final photon number: ", last(photon_number))
println("final inversion: ", last(inversion))
println("number of output times: ", length(tout))
