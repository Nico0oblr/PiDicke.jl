"""Return the integer `(S,M)` Dicke triangle used by the efficient backends."""
function dicke_triangle(N::Integer)
    iseven(N) || throw(ArgumentError("the integer Dicke triangle currently requires even N"))
    Jmax = Int(N) ÷ 2
    return [(S, M) for S in 0:Jmax for M in -S:S]
end

dicke_triangle(model::DickeModel) = dicke_triangle(model.N)

"""One representative `SpinBasis` for every allowed total-spin sector."""
function dicke_subspaces(N::Integer)
    js = iseven(N) ? collect(0:(Int(N) ÷ 2)) : collect(1//2:1:(Int(N)//2))
    return reverse(SpinBasis.(js))
end

"""Multiplicity-free direct sum of all two-level Dicke sectors."""
dicke_space(N::Integer) = directsum(dicke_subspaces(N)...)
dicke_space(model::DickeModel) = dicke_space(model.N)

function block_ranges(space::SumBasis)
    ranges = Vector{UnitRange{Int}}(undef, length(space.shape))
    start = 1
    for (i, len) in enumerate(space.shape)
        ranges[i] = start:(start + len - 1)
        start += len
    end
    return ranges
end

representation_basis(model::DickeModel, ::TrajectoryRepresentation) =
    (_require_even_N(model, TrajectoryRepresentation()); dicke_triangle(model))
representation_basis(model::DickeModel, ::PopulationRepresentation) =
    (_require_even_N(model, PopulationRepresentation()); dicke_triangle(model))
representation_basis(model::DickeModel, ::FirstOrderCoherenceRepresentation) =
    first_order_coherence_basis(model)
representation_basis(model::DickeModel, ::FullDensityRepresentation) = dicke_space(model)
representation_basis(model::DickeModel, rep::CavityDensityRepresentation) =
    dicke_space(model) ⊗ FockBasis(rep.cutoff)
