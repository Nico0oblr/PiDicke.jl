# PIDicke.jl

Two-level permutation-invariant open-system calculations in the Dicke basis.
The package collects the tested implementations that previously lived across
`SuperradiantLasing`, `SpinMax/FullPerm.jl`, and the cavity notebooks.

The same physical model is available at several computational levels:

| Representation | Stored object | Intended use |
|---|---|---|
| `TrajectoryRepresentation` | one `(S,M)` state | Gillespie trajectories |
| `PopulationRepresentation` | `p[S,M]` | sparse diagonal evolution and steady states |
| `FirstOrderCoherenceRepresentation` | `|S,M-1><S,M|` | `g¹`, spectra, and linewidths |
| `FullDensityRepresentation` | operator on all `SpinBasis(S)` sectors | full PI master equations |
| `CavityDensityRepresentation` | Dicke `SumBasis ⊗ FockBasis` | spin-cavity master equations |

The efficient population and coherence representations currently use integer
`(S,M)` labels and therefore require even `N`. The full `SumBasis`
representation supports even and odd `N`.

## Installation

From the Julia package manager:

```julia
pkg> activate /path/to/PIDicke.jl
pkg> instantiate
```

## Model

```julia
using PIDicke

model = DickeModel(
    100;
    collective_decay = 1.0,
    collective_pump = 0.0,
    local_decay = 0.1,
    local_pump = 0.2,
    local_dephasing = 0.05,
)
```

The rate convention is shared by every representation.

## Dicke populations

```julia
states, indices, G = population_generator(model)
p0 = population_initial_state(indices, (50, 50))
times = range(0.0, 5.0, length = 501)
result = solve_populations(model, p0, times)

pss = population_steady_state(G)
```

For stochastic dynamics:

```julia
times, states = run_population_trajectory(model; tmax = 5.0)
```

Burn-in and ergodic stationary samplers from the research scripts are
deliberately not part of this package. A steady distribution is obtained from
the sparse population generator.

## First-order coherence sector

```julia
pss = population_steady_state(model)
times = range(0.0, 20.0, length = 2001)
result = first_order_correlation(model, pss, times)

ω, spectrum = spectrum_even_fft(times, result.correlation)
linewidth = spectrum_fwhm(ω, spectrum)
```

The coherence label `(S,M)` means `|S,M-1><S,M|`.

The original stochastic coherence unraveling is also available without the
sample-preparation layer:

```julia
sample_times, stochastic_g1 = simulate_first_order_correlation_from_samples(
    model,
    S_samples,
    M_samples,
    times,
)
```

For the original one-long-trajectory ergodic preparation followed by the
coherence ensemble:

```julia
result = ergodic_first_order_correlation(
    model,
    times;
    nsamples = 10_000,
    equilibration_time = 20.0,
    sampling_window = 60.0,
    rng,
)

result.samples.S
result.samples.M
result.samples.sample_times
result.correlation
```

The stages remain independently accessible through
`ergodic_population_samples` and
`simulate_first_order_correlation_from_samples`.

## Full density matrix

The full representation is the multiplicity-free direct sum of all allowed
two-level spin sectors:

```julia
space = dicke_space(model)
full = full_liouvillian(model)

ψ0 = fully_excited_state(model)
times = range(0.0, 2.0, length = 201)
result = master_evolution(times, ψ0, model)
```

`d_N_J = 1` is intentional: one representative of every
multiplicity-equivalent spin sector is retained. The transition amplitudes are
the ones used by the original `FullPerm.jl` implementation.

Low-level operator construction is public:

```julia
Sm = collective_operator(space, :minus)
Sp = collective_operator(space, :plus)
Sz = collective_operator(space, :z)

space, jumps, rates, labels = collapse_operators(model)
```

## Cavity composition

```julia
system = cavity_system(
    model;
    cutoff = 6,
    g = 0.1,
    kappa = 1.0,
    cavity_drive = 0.0,
    spin_drive = 0.0,
)

L = cavity_liouvillian(
    model;
    cutoff = 6,
    g = 0.1,
    kappa = 1.0,
).liouvillian
```

The returned system exposes `hamiltonian`, `jumps`, `rates`, `initial`, `a`,
`Sm`, `Sp`, and `Sz`, so it can be passed directly to QuantumOptics.jl or
modified by a downstream project.

## Source provenance

- Dicke amplitudes: `SuperradiantLasing/PermRates.jl` and `SpinMax/FullPerm.jl`
- Population trajectories: `SuperradiantLasing/PopulationMC.jl`
- Population and coherence generators: `SuperradiantLasing/ODE.jl`
- Full operators: `SpinMax/FullPerm.jl`
- Cavity construction: `LasingQuantumOpticsImplementation.ipynb`
- Spectral utilities: `SuperradiantLasing/Spectrum.jl`

Plotting and mean-field comparisons remain outside the core package.

## Examples

The [examples directory](examples/README.md) contains executable programs for
all five representation levels: individual trajectories, sparse populations,
the first-order coherence band, full PI density matrices, and PI spin-cavity
systems.
