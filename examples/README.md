# Examples by representation

Every example is a small executable Julia program. Run one from the package
directory with, for example:

```sh
julia --project=. examples/trajectory_single_decay.jl
```

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

The examples deliberately avoid plotting dependencies. Their returned arrays
can be passed directly to a plotting package or used in a notebook.
