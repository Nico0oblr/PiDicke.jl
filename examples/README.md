# Examples by representation

Every example is a small executable Julia program. Run one from the package
directory with, for example:

```sh
julia --project=. examples/trajectory_single_decay.jl
```

Unless a file states otherwise, rates are normalized to an arbitrary reference
rate `Gamma0 = 1` and times are expressed as the dimensionless variable
`tau = Gamma0*t`. These are algorithm demonstrations, not experimental
parameter sets. The validation examples spell this normalization out in code.

## 1. Dicke-state trajectories

- `trajectory_single_decay.jl` — one superradiant decay trajectory and its
  emitted intensity.
- `trajectory_ensemble.jl` — ensemble-averaged excitation with collective
  decay, local pump, and local dephasing.

These use only a single current `(S,M)` label per trajectory.

## 2. Dicke populations

- `population_dynamics.jl` — deterministic sparse-generator evolution.
- `population_steady_state.jl` — direct stationary distribution and moments.
- `population_vs_trajectories.jl` — deterministic evolution compared with a
  Gillespie ensemble using the same model.

## 3. First-order coherence band

- `spectrum.jl` — deterministic `|S,M-1><S,M|` evolution, spectrum, and
  linewidth.
- `coherence_stochastic.jl` — stochastic coherence trajectories initialized
  from samples of the exact stationary population distribution.
- `coherence_ergodic_pipeline.jl` — one long population trajectory, ergodic
  stationary sampling, and the resulting stochastic coherence ensemble.

## 4. Full PI density matrix

- `full_density_matrix.jl` — construct and evolve the complete reduced master
  equation.
- `full_driven_spin.jl` — add a collective coherent drive to the PI master
  equation.
- `full_odd_N.jl` — demonstrate that the `SumBasis` backend also supports odd
  atom numbers.

## 5. PI spin plus cavity

- `cavity.jl` — construct the combined basis and explicit Liouvillian.
- `cavity_time_evolution.jl` — evolve a pumped spin ensemble coupled to a
  lossy cavity and evaluate photon/spin observables.

## Experimental: active population region

- `experimental/active_region_steady_state.jl` — prepare ergodic population
  samples, construct their padded rectangular region inside the Dicke triangle,
  and find its closed-boundary steady state with implicit-Euler iteration.

The experimental API is accessed as `PIDicke.Experimental` and is not part of
the top-level exported interface.

## Validation examples

- `validation/population_vs_full_liouvillian.jl` — extract the diagonal block
  of the full Liouvillian and compare its generator and propagated population
  with the efficient population representation.
- `validation/coherence_vs_full_liouvillian.jl` — perform the corresponding
  comparison for the `|S,M-1><S,M|` first-order coherence block and its
  correlation readout.
- `validation/bad_cavity_elimination.jl` — sweep `kappa/g` and demonstrate
  convergence of the explicit cavity model to collective decay with
  `Gamma_c = 4g^2/kappa`.
- `validation/ergodic_vs_exact_steady_state.jl` — compare one-trajectory
  ergodic histograms with the exact sparse population steady state as both the
  number of samples and observed trajectory duration increase.

`validation/full_block_helpers.jl` contains basis-extraction helpers shared by
the two deterministic full-Liouvillian comparisons.

The examples deliberately avoid plotting dependencies. Their returned arrays
can be passed directly to a plotting package or used in a notebook.
