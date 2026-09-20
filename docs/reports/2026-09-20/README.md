# Ember working paper — September 20, 2026

[Read the PDF](ember_status_and_future.pdf).

The calculation with finite convective transport from the helium-3 flash onset
reaches **69.76 thousand years**, **5024 K** and nuclear power **5.733e32 erg/s**
on **1983** mass points. Nuclear power continues to decline through independently checked intervals. A matched
**100-year** comparison on **1927/1983** points differs by **8.261%** in nuclear
heat and **15.24%** in final nuclear power. The first year after refinement has
an integrated energy residual of **1.844e43 erg**, associated with a mechanical
readjustment. Correcting only the initial radii leaves the one-year nuclear
heat and final power unchanged within **1e-12** relative, but retains the
explicit initial energy change. The complete flash and its spatial accuracy
remain unresolved.

An alternative flash history reaches **8.086 million years**, **5699 K** and
nuclear power **8.327e30 erg/s** on **2115** mass points. Nuclear burning supplies
**33.37%** of photon luminosity. Radius and surface luminosity have declined
**74.86%** and **95.34%** from their maxima. Effective temperature passes a
maximum of **6399 K** near **1.088 million years**, then falls **699.8 K** over
**6.999 million years**. Completed native calls over the final **2.1 million years** account for **89.89 CPU minutes**;
active elapsed time is **33.70 minutes**, including timestep comparisons.

A matched **1000-year** test near **log g = 5.42** changes surface temperature by **5.918 K** and deposited nuclear energy by **0.00005299%** between the two inferred hydrogen gravity responses, with identical final convective regions. This local comparison does not validate the absolute atmosphere or establish that ignition is independent of it.

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

The cooler atmosphere table contains **18** mixture columns over **5500–6500 K**
and **log g = 5.40–5.55**. A matched **800,000-year** comparison of inferred
hydrogen temperature responses changes the final temperature by **2.084 K**
and nuclear heat by **0.002024%**, with identical final convection. Independent
5750 K columns check the composition response between temperature and gravity
nodes. The absolute atmosphere prescription remains an uncertainty.

The gravity extension uses **18** mixture columns over **5500–6000 K** and
**log g = 5.40–5.80**. A matched **2.1-million-year** comparison of two hydrogen
references changes endpoint temperature by **2.292 K**, surface luminosity by
**0.02832%** and nuclear heat by **0.002680%**, with identical final convection.
Independent **5750 K**, **log g = 5.675** columns check the interpolated
composition response to **0.001026%** in temperature and **0.03026%** in pressure.
These local checks do not establish the absolute atmosphere prescription.


The retained pre-flash trajectory accounts for **15.13 elapsed hours** and
**36.31 CPU-hours**, excluding source construction and inter-run pauses.
A complete full-physics cooling runtime has not yet been measured.

This directory contains the PDF, LaTeX inputs and included figure PDFs.
The [late cooling figure](late_cooling.pdf), its [data](late_cooling.csv) and
plotting script are included. The LBA97 opacity comparison is retained.
The [standalone HR diagram of both histories](../../figures/2026-09-20/ember_two_sequences_hr.png)
shows the alternatives through 69.76 thousand years and 8.086 million years;
standalone flash diagrams remain separate from this draft.

Rebuild here with `tectonic --keep-logs ember_status_and_future.tex`.
