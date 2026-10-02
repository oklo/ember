# Measuring computational cost

Report stellar evolution and physics-table generation separately. Prepared
EOS, opacity and atmosphere tables are reused across time steps and runs.

- **CPU time:** summed user and system time across the process and its workers.
  One CPU-hour is one core working for one hour.
- **Wall time:** elapsed time to completion. State whether machine suspension
  is included and record competing jobs and thread counts.
- **Worker-hours:** use CPU-hours instead unless measuring allocated capacity;
  an allocated worker may be idle.
- **FLOPs:** label estimates, state precision and the measured throughput used.
  Hardware peak throughput is not the achieved rate of the stellar solver.

Compare the same starting state, physical inputs, target age and accuracy.
Include rejected trials and startup cost. Report changes in global quantities,
fuel inventories, profiles and convective boundaries alongside the speed ratio.
A short gradual-burning benchmark does not measure a full lifetime or flash cost.

A moving convective boundary can make the estimated composition error jump
as it crosses a mesh cell. A [matched 1 Gyr comparison](results/warm_species_time_tolerance_oct2.json)
uses a stated larger time-error allowance and reduces 16 accepted intervals
to three, with a 4.816-fold wall-time saving. Conservation tolerances remain
unchanged; local profile differences and the limited evolutionary range of
the test are recorded alongside the timing.

The solver reduces repeated work with structure prediction, composition warm
starts, local burning feedback and bounded reuse of EOS and collision responses.
Independent zone work uses CPU threads. Source atmosphere calculations can run
in parallel; the stellar solver does not currently use the GPU.

When an absent species needs a positive Newton guess, the solver mixes complete
compositions. It limits increases of existing species to 1%, avoiding a large
temporary fuel supply in a depleted core. CN mass accounting allows only
relative floating-point roundoff in a zero inert-metal remainder. These fixes
recover a failed 40 Myr cooling interval with unchanged time, structure and
conservation tolerances.
[Regression and stellar checks](results/cold_composition_start_oct1_v1.json).

Degenerate-electron collision integrals concentrate quadrature points near the
Fermi surface. A saved-state batch uses about half the CPU time, with transport
changes below 1.447e-9. A 200 Myr cooling comparison confirms negligible changes
in structure and composition. This measures the collision calculation, not a
factor-of-two speedup of the whole track.
[Numerical checks](results/fermi_surface_quadrature_oct1_v1.json).

Configuration and validation records:
[driver options](LIFETIME_DRIVER.md),
[starting guesses](results/solver_starting_guesses_sept27_v1.json),
[collision reuse](results/collision_reuse_sept27_v1.json),
[time integration](results/richardson_cooling_sept27_v1.json).
