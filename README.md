# PIDicke.jl

Two-level permutation-invariant open-system calculations in the Dicke basis.
This is the group reference implementation: it collects the tested code that
previously lived across `SuperradiantLasing`, `SpinMax/FullPerm.jl`, and the
cavity notebooks behind one set of conventions.

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

## Conventions

- `N` is the number of identical two-level systems and `S` is total spin. For
  the efficient even-`N` representations, `S = 0,1,...,N/2` and
  `M = -S,-S+1,...,S`.
- Population vectors obey `dp/dt = G*p`. A column of `G` is a source state and
  a row is a destination state; consequently a closed generator has zero
  column sums.
- The five fields of `DickeModel` are nonnegative Lindblad rates. The full
  representation uses the QuantumOptics.jl convention
  `rate * (J*rho*J' - (J'*J*rho + rho*J'*J)/2)`. The population and trajectory
  event rates are the same rate multiplied by the squared Dicke amplitude.
- `collective_decay` and `collective_pump` multiply `S-` and `S+`. The three
  local rates are represented by the PI transitions between neighboring total
  spin sectors used in the original code.
- The first-order coherence label `(S,M)` denotes
  `|S,M-1><S,M|`. Correlations are read out with the matrix element of `S-`.
- The full density representation keeps one representative of each total-spin
  sector. Thus `d_N_J = 1` is an operational multiplicity convention, not the
  physical combinatorial multiplicity of that sector.
- Units are arbitrary but must be consistent. With rates in inverse time, all
  evolution times use the corresponding time unit. QuantumOptics.jl uses
  `hbar = 1`.
- `FockBasis(cutoff)` contains occupations `0:cutoff`, so the cavity dimension
  is `cutoff + 1`.
- `kappa` is the coefficient of the standard dissipator `D[a]`. Photon number
  therefore decays at `kappa`, while the empty-cavity field amplitude decays at
  `kappa/2`. On resonance, eliminating a bad cavity with Hamiltonian coupling
  `g * (S+*a + S-*a')` gives the collective spin decay rate
  `Gamma_c = 4g^2/kappa`.

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

The exact sparse steady distribution is available through
`population_steady_state`. For large systems, the original one-trajectory
ergodic sampler is available through `ergodic_population_samples`.

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

## Experimental active population region

The large-`N` active-region steady-state calculation is available under the
explicitly experimental namespace `PIDicke.Experimental`. Its region shape is
the implementation that worked in `ActiveSpacePlots.ipynb`: the extrema of the
ergodic samples define a rectangle in `(S,M)`, the same padding is added on all
sides, and that rectangle is intersected with the Dicke triangle. It does not
use outgoing-transition shells.

```julia
const Active = PIDicke.Experimental

samples = ergodic_population_samples(
    model,
    10_000;
    equilibration_time = 20.0,
    sampling_window = 60.0,
    rng,
)

result = Active.active_region_steady_state(
    model,
    samples;
    padding = 50,
    dt = 1e3,
    nsteps = 5,
)
```

The restricted boundary is closed, the initial population is the normalized
sample histogram, and the steady state is approached through repeated implicit
Euler steps. Inspect `result.diagnostics.mean_leakage` when deciding whether
the padding is large enough. `dt` uses the model's time units: the notebooks
used `1e3` when the reference decay rate was one and `1e3 / Γ` when rates
were expressed using an explicit `Γ`. This API is intentionally not exported
from the top-level module while the region geometry is still being developed.

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
- Experimental active region and implicit-Euler iteration:
  `ActiveSpacePlots.ipynb` and `SteadyState.jl`

Plotting and mean-field comparisons remain outside the core package.

## Numerical validation

Small deterministic regression tests extract the population and first-order
coherence blocks from the full PI Liouvillian and compare them element by
element with the efficient generators. Both comparisons currently agree at
machine precision. A separate `N=1` regression checks that the explicit cavity
approaches collective spin decay with `Gamma_c = 4g^2/kappa` in the bad-cavity
limit.

The executable programs in `examples/validation` expose the same comparisons,
including propagated states, a sweep of `kappa/g`, and convergence of an
ergodic histogram toward the exact population steady state.

## Examples

The [examples directory](examples/README.md) contains executable programs for
all five representation levels: individual trajectories, sparse populations,
the first-order coherence band, full PI density matrices, and PI spin-cavity
systems.
