# Full multiplicity-free density-matrix representation ported from FullPerm.jl
# and LasingQuantumOpticsImplementation.ipynb.

"""Operational multiplicity convention used by the original implementation."""
d_N_J(::Int, ::Real) = 1.0

"""
    build_dicke_operator(N, matrix_element; space=dicke_space(N))

Build a sparse operator on the Dicke `SumBasis`. `matrix_element` has signature
`(N, S, M, Sp, Mp)` and returns `<Sp,Mp|O|S,M>`.
"""
function build_dicke_operator(
    N::Int,
    matrix_element::Function;
    space::SumBasis = dicke_space(N),
)
    subspaces = space.bases
    ranges = block_ranges(space)
    matrix = spzeros(ComplexF64, length(space), length(space))
    for (source_sector, source_basis) in enumerate(subspaces)
        S = source_basis.spinnumber
        for (target_sector, target_basis) in enumerate(subspaces)
            Sp = target_basis.spinnumber
            for source_index in 1:length(source_basis)
                M = S - (source_index - 1)
                for target_index in 1:length(target_basis)
                    Mp = Sp - (target_index - 1)
                    value = matrix_element(N, S, M, Sp, Mp)
                    value == 0 && continue
                    matrix[ranges[target_sector][target_index], ranges[source_sector][source_index]] = value
                end
            end
        end
    end
    return SparseOperator(space, matrix)
end

_me_collective_decay(N, S, M, Sp, Mp) =
    S == Sp && Mp == M - 1 ? sqrt(d_N_J(N, S)) * A_JM_minus(S, M) : 0.0
_me_collective_pump(N, S, M, Sp, Mp) =
    S == Sp && Mp == M + 1 ? sqrt(d_N_J(N, S)) * A_JM_plus(S, M) : 0.0

_me_local_decay_zero(N, S, M, Sp, Mp) =
    S == Sp && Mp == M - 1 ? sqrt(d_N_J(N, S)) * P_JM_minus_0(S, M, N) : 0.0
_me_local_decay_minus(N, S, M, Sp, Mp) =
    Sp == S - 1 && Mp == M - 1 ? sqrt(d_N_J(N, S)) * P_JM_minus_minus(S, M, N) : 0.0
_me_local_decay_plus(N, S, M, Sp, Mp) =
    Sp == S + 1 && Mp == M - 1 ? sqrt(d_N_J(N, S)) * P_JM_minus_plus(S, M, N) : 0.0

_me_local_pump_zero(N, S, M, Sp, Mp) =
    S == Sp && Mp == M + 1 ? sqrt(d_N_J(N, S)) * P_JM_plus_0(S, M, N) : 0.0
_me_local_pump_minus(N, S, M, Sp, Mp) =
    Sp == S - 1 && Mp == M + 1 ? sqrt(d_N_J(N, S)) * P_JM_plus_minus(S, M, N) : 0.0
_me_local_pump_plus(N, S, M, Sp, Mp) =
    Sp == S + 1 && Mp == M + 1 ? sqrt(d_N_J(N, S)) * P_JM_plus_plus(S, M, N) : 0.0

_me_dephasing_zero(N, S, M, Sp, Mp) =
    S == Sp && Mp == M ? sqrt(d_N_J(N, S)) * P_JM_z_0(S, M, N) : 0.0
_me_dephasing_minus(N, S, M, Sp, Mp) =
    Sp == S - 1 && Mp == M ? sqrt(d_N_J(N, S)) * P_JM_z_minus(S, M, N) : 0.0
_me_dephasing_plus(N, S, M, Sp, Mp) =
    Sp == S + 1 && Mp == M ? sqrt(d_N_J(N, S)) * P_JM_z_plus(S, M, N) : 0.0

function collective_operator(space::SumBasis, component::Symbol)
    operators = if component === :minus
        sigmam.(space.bases)
    elseif component === :plus
        sigmap.(space.bases)
    elseif component === :x
        sigmax.(space.bases)
    elseif component === :y
        sigmay.(space.bases)
    elseif component === :z
        sigmaz.(space.bases)
    elseif component === :identity
        identityoperator.(space.bases)
    else
        throw(ArgumentError("unknown collective component $component"))
    end
    if length(operators) == 1
        return SparseOperator(space, sparse(first(operators).data))
    end
    return directsum(operators...)
end

collective_operator(model::DickeModel, component::Symbol) = collective_operator(dicke_space(model), component)

"""Return `(space, jumps, rates, labels)` for the full density representation."""
function collapse_operators(model::DickeModel; include_zero_rates::Bool = false)
    N = model.N
    space = dicke_space(model)
    specifications = (
        (:collective_decay, _me_collective_decay, model.collective_decay),
        (:collective_pump, _me_collective_pump, model.collective_pump),
        (:local_decay_zero, _me_local_decay_zero, model.local_decay),
        (:local_decay_minus, _me_local_decay_minus, model.local_decay),
        (:local_decay_plus, _me_local_decay_plus, model.local_decay),
        (:local_pump_zero, _me_local_pump_zero, model.local_pump),
        (:local_pump_minus, _me_local_pump_minus, model.local_pump),
        (:local_pump_plus, _me_local_pump_plus, model.local_pump),
        (:local_dephasing_zero, _me_dephasing_zero, model.local_dephasing),
        (:local_dephasing_minus, _me_dephasing_minus, model.local_dephasing),
        (:local_dephasing_plus, _me_dephasing_plus, model.local_dephasing),
    )

    jumps, rates, labels = Any[], Float64[], Symbol[]
    for (label, matrix_element, rate) in specifications
        !include_zero_rates && iszero(rate) && continue
        push!(jumps, build_dicke_operator(N, matrix_element; space))
        push!(rates, Float64(rate))
        push!(labels, label)
    end
    return space, jumps, rates, labels
end

function _zero_hamiltonian(space)
    return SparseOperator(space, spzeros(ComplexF64, length(space), length(space)))
end

function full_liouvillian(model::DickeModel; hamiltonian = nothing)
    space, jumps, rates, labels = collapse_operators(model)
    H = isnothing(hamiltonian) ? _zero_hamiltonian(space) : hamiltonian
    L = QuantumOptics.liouvillian(H, jumps; rates)
    return (; space, hamiltonian = H, jumps, rates, labels, liouvillian = L)
end

function master_evolution(
    times,
    initial,
    model::DickeModel;
    hamiltonian = nothing,
    kwargs...,
)
    space, jumps, rates, labels = collapse_operators(model)
    H = isnothing(hamiltonian) ? _zero_hamiltonian(space) : hamiltonian
    tout, states = timeevolution.master(times, initial, H, jumps; rates, kwargs...)
    return (; space, hamiltonian = H, jumps, rates, labels, times = tout, states)
end

function _polarized_state(N::Int, direction::Symbol)
    space = dicke_space(N)
    state = direction === :up ? spinup(first(space.bases)) : spindown(first(space.bases))
    data = zeros(ComplexF64, length(space))
    data[first(block_ranges(space))] .= state.data
    return Ket(space, data)
end

fully_excited_state(N::Int) = _polarized_state(N, :up)
fully_excited_state(model::DickeModel) = fully_excited_state(model.N)
fully_ground_state(N::Int) = _polarized_state(N, :down)
fully_ground_state(model::DickeModel) = fully_ground_state(model.N)
