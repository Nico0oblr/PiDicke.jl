"""Marker type for a numerical representation of the two-level PI dynamics."""
abstract type AbstractDickeRepresentation end

"""A single stochastic trajectory on the Dicke triangle `(S,M)`."""
struct TrajectoryRepresentation <: AbstractDickeRepresentation end

"""Diagonal Dicke populations `p[S,M]`."""
struct PopulationRepresentation <: AbstractDickeRepresentation end

"""The band `|S,M-1><S,M|` used for first-order correlations."""
struct FirstOrderCoherenceRepresentation <: AbstractDickeRepresentation end

"""The full density operator on the multiplicity-free Dicke `SumBasis`."""
struct FullDensityRepresentation <: AbstractDickeRepresentation end

"""The full Dicke `SumBasis` tensored with a truncated cavity basis."""
struct CavityDensityRepresentation <: AbstractDickeRepresentation
    cutoff::Int
end
