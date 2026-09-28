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

The solver reduces repeated work with structure prediction, composition warm
starts, local burning feedback and bounded reuse of EOS and collision responses.
Independent zone work uses CPU threads. Source atmosphere calculations can run
in parallel; the stellar solver does not currently use the GPU.

Configuration and validation records:
[driver options](LIFETIME_DRIVER.md),
[starting guesses](results/solver_starting_guesses_sept27_v1.json),
[collision reuse](results/collision_reuse_sept27_v1.json),
[time integration](results/richardson_cooling_sept27_v1.json).
