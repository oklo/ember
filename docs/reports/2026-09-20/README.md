# Ember working paper — September 20, 2026

[Read the PDF](ember_status_and_future.pdf).

The calculation with finite convective transport from the helium-3 flash onset
reaches **16.85 thousand years**, **4794 K** and nuclear power **4.061e34 erg/s**
on **1899** mass points. Its final **535.7 years** release **5.878e44 erg**;
all **91** accepted intervals pass independent timestep and conservation checks.
One mass cell supplies **20.03%** of nuclear power at the endpoint, requiring
further local refinement before continuation. The calculation has passed a
local nuclear-power maximum, a decline and renewed burning; its largest peak
and complete pulse remain unresolved.

An alternative flash history reaches **3.054 million years**, **6205 K** and
nuclear power **2.231e31 erg/s** on **2115** mass points. Nuclear burning supplies
**33.92%** of photon luminosity. Radius and surface luminosity have declined
**65.58%** and **87.71%** from their maxima. Effective temperature passes a
maximum of **6399 K** near **1.088 million years**, then falls **193.7 K** over
**1.966 million years**. The final **400,000 years** use **17.69 CPU minutes**
and **6.257 minutes of active elapsed time**, including timestep comparisons.

**These are alternative histories from a common model, not consecutive pieces
of one verified track.** Their mixing treatments during the flash differ.
Checks of later evolution do not establish the accuracy of the earlier flash.
Both histories inherit the supplied hydrogen atmosphere reference; whether a
flash occurs under a different atmosphere prescription remains untested.
Its appearance about **0.7546 Gyr** after the switch distinguishes it from the
immediate surface readjustment, but does not establish independence from the
subsequent envelope and fuel history. The short atmosphere comparison is now
explicitly described as exploratory because one calculation exceeds its
timestep accuracy limit.

The late atmosphere grid contains **54** mixture columns over **6000–6500 K**
and **log g = 4.20–5.40**. Source, depth, native and matched stellar checks pass.
A **273,200-year** alternative-reference comparison changes temperature by
**0.09741 K**, surface luminosity by **0.004235%** and nuclear heat by **0.01126%**,
with identical final convection. Its shared **126,800-year** prefix is reused.
This measures local sensitivity after ignition; the absolute atmosphere
prescription remains an uncertainty.

The retained pre-flash trajectory accounts for **15.13 elapsed hours** and
**36.31 CPU-hours**, excluding source construction and inter-run pauses.
A complete full-physics cooling runtime has not yet been measured.

This directory contains the PDF, LaTeX inputs and included figure PDFs.
The [late cooling figure](late_cooling.pdf), its [data](late_cooling.csv) and
plotting script are included. The LBA97 opacity comparison is retained.
The [standalone HR diagram of both histories](../../figures/2026-09-20/ember_two_sequences_hr.png)
shows the alternatives through 16.31 thousand years and 2.654 million years;
standalone flash diagrams remain separate from this draft.

Rebuild here with `tectonic --keep-logs ember_status_and_future.tex`.
