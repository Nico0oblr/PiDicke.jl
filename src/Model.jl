"""
    DickeModel(N; collective_decay=0, collective_pump=0,
               local_decay=0, local_pump=0, local_dephasing=0)

Rates for a permutation-invariant ensemble of `N` two-level systems. The same
model is accepted by the trajectory, population, coherence-band, and full
density-matrix representations.
"""
struct DickeModel{T<:Real}
    N::Int
    collective_decay::T
    collective_pump::T
    local_decay::T
    local_pump::T
    local_dephasing::T
end

function DickeModel(
    N::Integer;
    collective_decay::Real = 0.0,
    collective_pump::Real = 0.0,
    local_decay::Real = 0.0,
    local_pump::Real = 0.0,
    local_dephasing::Real = 0.0,
)
    N > 0 || throw(ArgumentError("N must be positive"))
    rates = promote(
        collective_decay,
        collective_pump,
        local_decay,
        local_pump,
        local_dephasing,
    )
    all(>=(zero(first(rates))), rates) || throw(ArgumentError("rates must be nonnegative"))
    return DickeModel(Int(N), rates...)
end

jmax(model::DickeModel) = model.N // 2

function model_parameters(model::DickeModel)
    return (
        collective_decay = model.collective_decay,
        collective_pump = model.collective_pump,
        local_decay = model.local_decay,
        local_pump = model.local_pump,
        local_dephasing = model.local_dephasing,
    )
end

function _require_even_N(model::DickeModel, representation::AbstractDickeRepresentation)
    iseven(model.N) || throw(ArgumentError(
        "$(typeof(representation)) currently uses integer (S,M) labels and requires even N; " *
        "FullDensityRepresentation supports both even and odd N",
    ))
    return nothing
end
