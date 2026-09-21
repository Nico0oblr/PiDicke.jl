"""
    cavity_system(model; cutoff, g, kappa, cavity_drive=0, spin_drive=0)

Compose the full PI spin `SumBasis` with a truncated cavity. This is the
refactored two-level construction from `LasingQuantumOpticsImplementation.ipynb`.
The returned named tuple contains the bases, Hamiltonian, jumps, rates, initial
state, and common observables.
"""
function cavity_system(
    model::DickeModel;
    cutoff::Int,
    g::Real,
    kappa::Real,
    cavity_drive::Real = 0.0,
    spin_drive::Real = 0.0,
)
    cutoff >= 1 || throw(ArgumentError("cutoff must be at least one"))
    kappa >= 0 || throw(ArgumentError("kappa must be nonnegative"))

    spin = dicke_space(model)
    cavity = FockBasis(cutoff)
    id_spin = collective_operator(spin, :identity)
    id_cavity = identityoperator(cavity)
    lowering = collective_operator(spin, :minus)
    raising = collective_operator(spin, :plus)

    a = id_spin ⊗ destroy(cavity)
    adag = id_spin ⊗ create(cavity)
    Sm = lowering ⊗ id_cavity
    Sp = raising ⊗ id_cavity

    hamiltonian = g * (Sp * a + Sm * adag)
    if !iszero(cavity_drive)
        hamiltonian += cavity_drive * (a + adag)
    end
    if !iszero(spin_drive)
        hamiltonian += spin_drive * (Sp + Sm)
    end

    _, spin_jumps, spin_rates, spin_labels = collapse_operators(model)
    jumps = Any[jump ⊗ id_cavity for jump in spin_jumps]
    rates = copy(spin_rates)
    labels = copy(spin_labels)
    if kappa > 0
        push!(jumps, a)
        push!(rates, Float64(kappa))
        push!(labels, :cavity_decay)
    end

    initial = fully_ground_state(model) ⊗ coherentstate(cavity, 0.0)
    Sz = collective_operator(spin, :z) ⊗ id_cavity
    return (;
        spin,
        cavity,
        basis = spin ⊗ cavity,
        hamiltonian,
        jumps,
        rates,
        labels,
        initial,
        a,
        adag,
        Sm,
        Sp,
        Sz,
    )
end

function cavity_liouvillian(model::DickeModel; kwargs...)
    system = cavity_system(model; kwargs...)
    L = QuantumOptics.liouvillian(system.hamiltonian, system.jumps; rates = system.rates)
    return merge(system, (; liouvillian = L))
end
