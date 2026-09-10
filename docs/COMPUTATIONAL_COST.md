# What is consuming the computation time?

## Further reuse of nuclear calculations

The evolution driver now keeps up to **8192 complete nuclear results per
worker**, compared with 16 previously. It reuses a result only when temperature,
density, every abundance and both composition labels agree exactly. This avoids
repeating calculations across the stellar mesh. Reaction formulas and precision
are unchanged.

In a 512-point calculation through 100 billion years, an alternating comparison
measured **42.21% less CPU time**: 87.57 seconds before versus 50.61 seconds after.
All four full output files are byte-identical. Elapsed means were 121.9 and 65.50
seconds, but concurrent source jobs changed during the comparison; these wall
times are not an isolated-machine benchmark. Both CPU-time pairs were consistent.
[Measurements](results/nuclear_response_cache_benchmark_v4.json).

All 32 test suites pass for the candidate. The production build reproduces the
short output exactly and passes its native-restart test. Separate longer runs
with the old and new implementation are checking the complete covered trajectory.
The gain during shell burning and cooling remains unmeasured.
[Validation](results/nuclear_response_cache_validation_v4.json).

## Fixed-input production measurements

Fresh calculations with unchanged physical inputs completed the first 3.600
trillion years of the 0.1 Msun star:

| Calculation | Awake elapsed | Total CPU | Accepted steps |
|---|---:|---:|---:|
| 512 points, fourfold tighter time control | 39.78 min | 64.06 min | 5165 |
| 1024 points, original time control | 59.26 min | 96.15 min | 2782 |

Both use two CPU workers. Their UTC elapsed times include an additional 40.88
minutes of computer suspension. Source calculations and other jobs are excluded
from their CPU times. The [receipts and accuracy comparison](results/evolution_accuracy_3600gyr_v1.json)
record input and executable identities. Different mesh and time controls make
these useful production measurements, not an isolated mesh-scaling benchmark.
The full track through helium-remnant cooling has not yet been timed.

The subsequent 512-point calculation reached 3.685 trillion years before the
atmosphere's 3400 K boundary stopped it. That is an input-coverage limit.
Generating and validating hotter atmosphere models extends reusable inputs;
ordinary stellar steps interpolate those tables. The new trajectory uses fixed
132-model atmosphere and hot-core opacity tables, with the same timing receipt.

A further attempt to share the electron calculation between pp reactions gave
identical outputs but **did not establish a speedup**. An alternating comparison
measured 73.31 seconds for the candidate versus 54.31 seconds for the unchanged
code; CPU times were 94.82 versus 68.27 seconds. Other jobs and changing machine
conditions contributed timing variation. The candidate remains isolated and is
not installed. [Measurements](results/screening_reuse_benchmark_v2.json).

## Changes measured on September 10

Two changes now reduce repeated work and use an additional CPU core:

- A conduction mixture reuses each source-ion interpolation at its fixed
  temperature and electron density. The mixture sums and their derivatives
  retain their original order. A replay of the actual 3.600-trillion-year
  profile runs **2.773 times faster in the conduction calculation alone**,
  with identical opacity, derivative and accumulated-sum output.
- `ember-evolve --step-workers 2` computes the independent full-step estimate
  concurrently with the two sequential half steps. Each worker has its own
  nuclear and atmosphere caches; immutable physical tables are shared. The
  default is one worker. Accuracy tests and accepted-state ordering are unchanged.

In an alternating before/after/after/before comparison, a 512-point evolution
through 100 billion years took **40.40 seconds before and 27.16 seconds after**
on average: **32.76% less elapsed time**, or **1.487 times faster**. Total CPU
use fell by **6.371%**. Every complete output file matched byte for byte,
including all 48 accepted steps, final structure and abundances. Startup and
initial-model calculation are included; this is not a measured full-lifetime
speedup. The short earlier one-billion-year benchmark spends much of its time
starting the calculation and found only a small improvement.

All 32 test suites pass. Added restart checks compare one- and two-worker
histories and checkpoint bytes exactly, including a change of worker count
when resuming the same saved model. The physical precision and source-domain
limits are unchanged. Measurements:
[whole stellar calculation](results/duplication_parallel_benchmark_v1.json),
[conduction replay](results/conduction_profile_reuse_v1.json).

This machine is an Apple M4 Max with 10 performance and four efficiency CPU
cores, a 32-core GPU and 36 GB of memory. The two-worker mode uses additional
CPU capacity within one star. Independent atmosphere models and convergence
controls can use other cores. Ages within one trajectory remain sequential;
using all cores requires additional work on independent spatial calculations
and their caches. No all-core speedup has yet been measured.

Ember does not currently use Metal. Apple's current
[Metal language specification](https://developer.apple.com/metal/Metal-Shading-Language-Specification.pdf)
excludes the `double` type used by the stellar solver. A GPU implementation
would require a validated precision strategy or software arithmetic. Apple
also explains why small, frequently synchronized workloads can scale poorly
in [its compute guidance](https://developer.apple.com/videos/play/wwdc2022/10159/).
The existing C++ build targets the M4 CPU, permits fused arithmetic, and links
Accelerate; its small structure-solver blocks use the portable kernel.

The electron-screening experiments above have not justified an additional
change. Reuse between neighboring zones and parallel spatial calculations
remain candidates, each requiring an actual trajectory comparison.

## Earlier stellar measurement — September 10

A five-second sample of the running 512-point star after it developed a stable
central region attributes about 48.81% of sampled time to nuclear calculations,
including 40.46% in the electron integrals used for reaction screening.
Opacity and conduction account for about 36.75%, and the equation of state
for 9.650%. These are inclusive sampled call counts; electron work is already
included in nuclear work. This short interval does not establish a lifetime
speedup. The [measurement](results/stellar_cpu_profile_20260910_v1.json)
records the source sample, executable identity and aggregation definitions.

At the progress check the run had advanced from 3.560 to 3.590 trillion years
in 4.896 minutes, accepting 216 steps and rejecting eight attempts. It averaged
1.360 seconds per accepted step including retries; recent steps advanced about
156.2 million years. The complete 3.560–3.600-trillion-year segment subsequently finished in
6.250 CPU minutes, with 277 accepted steps and ten rejected attempts.
The process used essentially one full CPU core. It was
neither waiting for atmosphere tables nor stuck in repeated failures.

The earlier lack of progress beyond 3.560 trillion years included an avoidable
scheduling gap: both requested segments had finished and no continuation was
launched during the PDF update. That elapsed delay must not be attributed to
the stellar solver. The new continuation was launched after the user asked.

One spatial dimension keeps the structure equations small. The present cost
comes from repeatedly evaluating material properties across the mesh while
burning, mixing and structure are iterated together. Every step also computes
one full step and two half steps to estimate its time error. The recent core
transition needs much smaller time steps than the long fully convective phase.
Useful optimization targets are repeated electron screening and conductive
interpolation. A previously tried electron cache changed one output bit and
was rejected; no new speedup is claimed or installed from this CPU sample.

The following notes preserve older timings and their different scopes.

Validated atmosphere tables are reusable inputs. Stellar timesteps interpolate
their matching states; they do not rerun wavelength-dependent radiative transfer.
The source calculations recur when coverage is extended, the physical prescription
changes, or independent resolution and interpolation controls are needed. A new
stellar track inside the same supported domain reuses the existing tables.

## Observed jobs on 2026-09-10

At 11:29 UTC, eight source/stellar processes each used approximately one CPU
core: one stellar evolution, two progressing XH=.1 atmosphere calculations,
four older atmosphere initializers repeating two-cycles, and one slow direct
canonical atmosphere control. The four repeated cycles persisted through
113--135 completed iterations. They were consuming resources without improving
convergence. The separate direct control's corrections were still changing;
it was not established to be in the same cycle. Parent-owned controllers were
left untouched. Cancellation receipts on the old plans prevent future source
launches, not completion of an already active child.

All twelve XH=.1 warm nodes subsequently completed using the finer opacity
derivative initializer and the unchanged canonical final solver. Two additional
high-gravity depth controls were launched separately because the 3200 K source
models ended at actual Rosseland optical depths165/170, while matching at100.
These controls measure a specific unresolved numerical effect; they are not
duplicate production runs.

The fresh512-zone stellar run reached3.30 trillion years in2351 child CPU
seconds (39.18 minutes) and4088 UTC elapsed seconds (68.14 minutes), with1822
accepted steps and3 rejected attempts. Concurrent jobs and possible suspension
make elapsed time unsuitable as an isolated-machine speed benchmark. Its output
is `out/evolution-cold-remnant-forward-512-3300gyr-v1.json`, with an adjacent
receipt; the copied executable is retained in
`/tmp/ember-cold-remnant-forward-3300gyr-v1`. Physical input hashes remained
unchanged. The receipt records source-tree edits from a side conversation, but
the running copied executable remained independent of those edits.

## Comparison with the fast F77 calculation

The local `low_mass_stars/f77/PLAN.md` records an18-second0.1-Msun calculation
from Hayashi contraction through exhaustion and cooling to about1251 K. That
is a real recorded baseline, although the note explicitly calls the cold
cooling clock provisional. It is not a100 K calculation. The local source
and the retained `/tmp/ember-f77-baseline-henyey77.f` show these differences:

- F77 obtains its surface condition through a grey atmosphere/envelope
  integration and uses existing mean-opacity tables. Ember's offline gas
  sources solve wavelength-dependent transfer at300 depths and20000
  frequencies, with composition-specific material opacities. Line-list
  synthesis is a separate preparation stage; it is not repeated at every
  stellar timestep or every atmosphere iteration.
- The F77 arrays allow201 interior points and include adaptive insertion and
  removal. The present Ember track uses512 fixed mass points.
- F77 burns and mixes composition once before each structure solve. Ember
  iterates burning/mixing and structure to consistency, then compares one
  full timestep with two half timesteps before acceptance. See
  `apps/evolve.cpp` and `src/evolution.cpp`.
- The selected nuclear rates, EOS, opacity and conductive inputs differ.
  Ember includes SFII rate quadratures and finite-degeneracy screening;
  the F77 code uses fitted rates. F77 already has explicit He3, CNO, conduction
  and adaptive mesh capabilities, so Ember is not a strict physics superset.

These differences establish additional work, but do not quantitatively explain
the entire speed gap. A matched-input F77/Ember benchmark remains missing. The newer sampled
stellar profile above locates costs in Ember without isolating the full
cross-code speed difference. Higher cost is not itself evidence of greater
physical accuracy.

## Measured optimization evidence

The archived [atmosphere replay profile](results/nongrey_cpu_profile.json)
attributes54% of its measured phase CPU time to molecular chemistry. Its dense
300-by-300 solve used only0.0033 seconds. This is one older source replay,
not a profile of every atmosphere or of the stellar evolution process.

Two independent ABBA stellar benchmarks measured12.52% and10.84% CPU reductions
from exact nuclear-response and bare-rate caches, respectively. Each comparison
reproduced the complete output bytes. The product corresponds to about22.00%
less CPU for those benchmark conditions; it is not an isolated512-zone
full-lifetime measurement. Reports:
[response cache](results/nuclear_cache_benchmark_v1.json) and
[bare-rate cache](results/nuclear_rate_cache_benchmark_v1.json).

The next performance diagnosis should separate the stellar integrator from
source-table generation, use immutable binaries and identical inputs, and
measure EOS, nuclear, structure and coupling costs on an owned short run.
Source convergence failures also need bounded termination and retained receipts;
continuing an established two-cycle does not improve the table.

At11:47 UTC the old cycling processes and the slow direct control were no
longer listed. This turn did not signal them, and their disappearance alone
does not establish successful completion. The new interpolation controller
then accounted for only two active source workers. The3.30->3.40T exact
stellar restart had completed in97.00 CPU seconds with unchanged inputs.

The separate F77 reproduction used the retained source,151 initial points,
alpha1.36 and the AJR opacity option. Its unmodified control used27.37 CPU
seconds through the cooling stop, with original output bytes exactly matching
the print-only diagnostic. This is a new recorded timing, distinct from the
historical18-second note. The central transition and same-hydrogen comparisons
are in [the onset audit](results/f77_core_onset_comparison_v1.json) and
[the matched-state report](results/f77_ember_matched_hydrogen_v1.json).

At XH~.203 the AJR F77 and Ember radii and central temperatures are similar,
but F77 is 76.4% more luminous and 15.5% hotter at the surface after converting
luminosity to common solar units. The original comparison of each code's native
solar units gave 74.9%; the source records retain those original numbers and
now state both conventions. Within F77,
changing only the opacity selector to Ferguson raises luminosity54.5% at
this composition. The two pp implementations differ by8.8--9.9% in heating
at identical central state coordinates, while the evolved He3 inventories
also differ. These comparisons support atmosphere/opacity and the resulting
composition evolution as major candidates; they do not uniquely attribute
the F77/Ember discrepancy or establish either model's accuracy.

The direction of the surface difference is consistent with the grey-boundary
overestimates of Teff and luminosity discussed in section2 of
[Baraffe et al. (1998)](https://arxiv.org/pdf/astro-ph/9805009). That literature
comparison is context, not validation of Ember's current atmosphere family.
