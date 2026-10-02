# Reusing envelope material evaluations

The implicit envelope integrator holds pressure, temperature and composition
fixed while solving for radius and enclosed mass. It now evaluates the EOS and
opacity once per temperature trial. Gravity and the radiative/convective
gradient are still recomputed at every geometry iteration. Equations, source
tables, tolerances and error checks are unchanged; nothing is reused between
different temperature trials.

Twelve envelope controls spanning Hayashi contraction, main-sequence evolution
and cooling retain identical boundary values. Integration-only CPU time falls
from 1.044 to 0.3028 seconds for five warmer cases, and from 2.722 to 0.7745
seconds for seven cooler cases. Table-loading time is excluded from those sums.

A matched, fully checked stellar interval from 1.892 to 1.900 million years
uses 30.80 versus 10.61 CPU seconds, and 24.95 versus 9.069 wall seconds,
at two threads on an M4 Max. Both accept one interval without rejection;
the final checkpoint differs only in its executable fingerprint. This short
comparison does not establish a speed ratio for the complete stellar lifetime.

The envelope Gibbs, transport, initialization and runtime-packaging tests pass.
Local comparison records: `out/envelope-material-reuse-oct2-v1/`. Restarting
with the optimized executable requires an explicit, recorded executable
identity update; the saved physical state and all other input identities remain
unchanged. No checkpoint-compatibility check was removed.
