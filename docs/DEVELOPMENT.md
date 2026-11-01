# Working on Ember

Start with the equations and the scientific question. Keep physical
approximations visible in the configuration and in the resulting record.
Numerical convergence, agreement with an independent calculation and physical
validity are separate claims.

| Responsibility | Entry point |
| --- | --- |
| Structural block elimination | `src/henyey.cpp`, [method](HENYEY.md) |
| Evolution and time-step selection | `src/evolution.cpp`, `src/controller.cpp` |
| Lifetime configuration and restart | `apps/lifetime_driver.cpp`, [configuration](LIFETIME_DRIVER.md) |
| Material table readers and joins | `src/eos_*.cpp`, `src/opacity_table.cpp`, `include/ember/opacity_radiative.hpp` |
| Retained study verification and figures | `research/record/reproduce.py` |
| Offline source generation | `scripts/`, indexed by [physics inputs](DATA_REPRODUCTION.md) |

Use the smallest function that expresses a repeated operation. Name quantities
for their physical or numerical role, with units where the convention is not
obvious. Comments should explain assumptions, equations, coordinate changes
and reasons for numerical choices. Avoid restating assignments or promising
future features. Preserve arithmetic order when a change is intended only to
improve readability.

Run the existing tests that exercise the changed behavior. For material joins,
check values, derivatives, abundance conventions and domain boundaries. For the
Henyey solver, retain the independent dense solve and original-equation residual
checks. Do not replace scientific checks with tests that merely repeat the
implementation.

The dated scripts under `out/` are experiment records. They sometimes depend on
specific source copies or continuation chains. Keep those identities intact;
promote a workflow into maintained code only when its inputs and behavior are
understood. The compact research-record renderer deliberately reads saved data
without importing campaign supervisors.
