using LinearAlgebra
using OrdinaryDiffEq: Rodas5P
using PIDicke
using QuantumOptics: expect, fockstate, timeevolution, ⊗

# PIDicke uses kappa*D[a], so photon number decays at kappa, the field
# amplitude decays at kappa/2, and bad-cavity elimination gives
# Gamma_c = 4g^2/kappa.
N = 2
g = 1.0
tau = collect(range(0.0, 4.0, length = 161)) # tau = Gamma_c*t

println("kappa/g    max inversion error after cavity transient")
for kappa_over_g in (10.0, 20.0, 40.0)
    kappa = kappa_over_g * g
    Gamma_c = 4g^2 / kappa

    cavity_spin = DickeModel(N)
    cavity = cavity_system(cavity_spin; cutoff = N, g, kappa)
    cavity_initial = fully_excited_state(cavity_spin) ⊗ fockstate(cavity.cavity, 0)

    times = tau ./ Gamma_c
    _, cavity_states = timeevolution.master(
        times,
        cavity_initial,
        cavity.hamiltonian,
        cavity.jumps;
        rates = cavity.rates,
        alg = Rodas5P(autodiff = false),
        reltol = 1e-9,
        abstol = 1e-11,
    )
    cavity_inversion = real.(expect(cavity.Sz, cavity_states))

    eliminated_model = DickeModel(N; collective_decay = Gamma_c)
    eliminated = master_evolution(
        times,
        fully_excited_state(eliminated_model),
        eliminated_model;
        alg = Rodas5P(autodiff = false),
        reltol = 1e-9,
        abstol = 1e-11,
    )
    eliminated_Sz = collective_operator(eliminated_model, :z)
    eliminated_inversion = real.(expect(eliminated_Sz, eliminated.states))

    # Do not assess the Markov approximation during the initial cavity memory
    # time. In scaled time, five cavity lifetimes are 5*Gamma_c/kappa.
    first_comparison = searchsortedfirst(tau, 5Gamma_c / kappa)
    error = maximum(abs,
        cavity_inversion[first_comparison:end] -
        eliminated_inversion[first_comparison:end],
    ) / N
    println(rpad(kappa_over_g, 11), error)
end
