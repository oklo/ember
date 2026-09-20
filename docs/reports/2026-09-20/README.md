# Ember working paper — September 20, 2026

[Read the PDF](ember_status_and_future.pdf).

The calculation with finite convective transport from the helium-3 flash onset
reaches **19.36 thousand years**, **4818 K** and nuclear power **1.069e34 erg/s**
on **1983** mass points. Nuclear power declines through the final **2000 years**,
which release **9.335e44 erg** through **173** checked intervals. A matched
**100-year** comparison on **1927/1983** points differs by **8.261%** in nuclear
heat and **15.24%** in final nuclear power. The first year after refinement has
an integrated energy residual of **1.844e43 erg**, associated with a mechanical
readjustment. Correcting only the initial radii leaves the one-year nuclear
heat and final power unchanged within **1e-12** relative, but retains the
explicit initial energy change. The complete flash and its spatial accuracy
remain unresolved.

An alternative flash history reaches **3.454 million years**, **6157 K** and
nuclear power **1.948e31 erg/s** on **2115** mass points. Nuclear burning supplies
**33.40%** of photon luminosity. Radius and surface luminosity have declined
**67.08%** and **89.10%** from their maxima. Effective temperature passes a
maximum of **6399 K** near **1.088 million years**, then falls **242.0 K** over
**2.366 million years**. The final **400,000 years** use **15.30 CPU minutes**
and **5.511 minutes of active elapsed time**, including timestep comparisons.

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
shows the alternatives through 19.36 thousand years and 3.454 million years;
standalone flash diagrams remain separate from this draft.

Rebuild here with `tectonic --keep-logs ember_status_and_future.tex`.
