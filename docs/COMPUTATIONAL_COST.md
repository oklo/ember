# What is consuming the computation time?

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

The fresh512-zone stellar run reached3.30 trillion years in2350.57 child CPU
seconds (39.18 minutes) and4088.19 UTC elapsed seconds (68.14 minutes), with1822
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
the entire speed gap. No matched-input F77/Ember benchmark or current whole-star
CPU profile has been obtained. Higher cost is not itself evidence of greater
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
but F77 is74.9% more luminous and15.5% hotter at the surface. Within F77,
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
