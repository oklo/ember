# Ember working paper — September 20, 2026

[Read the PDF](ember_status_and_future.pdf).

The calculation with finite convective transport from the helium-3 flash onset
reaches **128.3 thousand years**, **5594 K** and nuclear power **2.393e32 erg/s**
on **1983** mass points. Nuclear power continues to decline through independently checked intervals. A matched
**100-year** comparison on **1927/1983** points differs by **8.261%** in nuclear
heat and **15.24%** in final nuclear power. The first year after refinement has
an integrated energy residual of **1.844e43 erg**, associated with a mechanical
readjustment. Correcting only the initial radii leaves the one-year nuclear
heat and final power unchanged within **1e-12** relative, but retains the
explicit initial energy change. The complete flash and its spatial accuracy
remain unresolved.

An alternative flash history reaches **17.91 million years**, **5183 K** and
nuclear power **4.608e30 erg/s** on **2115** mass points. Nuclear burning supplies
**40.25%** of photon luminosity. Radius and surface luminosity have declined
**79.42%** and **97.86%** from their maxima. Effective temperature passes a
maximum of **6399 K** near **1.088 million years**, then falls **1216 K** over
**16.82 million years**. Completed native calls over the final **3.740 million years** account for **159.3 CPU minutes**;
active elapsed time is **57.01 minutes**, including timestep comparisons.

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

The cooler atmosphere table contains **27** source columns over **5000–6000 K**
and **log g = 5.40–5.80**. Independent **5250 K** and **5750 K** columns test
interpolation between temperature and gravity nodes. At **5250 K**, the
composition response agrees to **0.000007587%** in temperature and **0.1582%**
in pressure.

A matched **2.1-million-year** cooling comparison over **5500–6000 K** changes
endpoint temperature by **2.292 K** and nuclear heat by **0.002680%** between
two hydrogen reference prescriptions. The detailed-onset comparison over
**10,000 years** changes the endpoint by **4.187 K** and nuclear heat by
**0.009269%**. Both comparisons retain identical final convective regions.
These local tests do not establish the absolute hydrogen boundary condition.

A further **10,000-year** comparison entering the **5500–6000 K** atmosphere
interval changes the detailed history's endpoint by **9.515 K**, surface
luminosity by **0.6680%**, and nuclear heat by **0.001667%**, with identical
final convective regions. This measures local sensitivity to the inferred
hydrogen reference; it does not validate its absolute normalization.


A further matched **2.000-million-year** comparison in the **5000–6000 K**
grid changes the endpoint by **6.020 K**, surface luminosity by **0.1435%**
and integrated nuclear heat by **0.009124%**. The surface convective boundary
differs by one mass interface. This tests local atmosphere sensitivity.


The retained pre-flash trajectory accounts for **15.13 elapsed hours** and
**36.31 CPU-hours**, excluding source construction and inter-run pauses.
A complete full-physics cooling runtime has not yet been measured.

This directory contains the PDF, LaTeX inputs and included figure PDFs.
The [late cooling figure](late_cooling.pdf), its [data](late_cooling.csv) and
plotting script are included. The LBA97 opacity comparison is retained.
The [standalone HR diagram of both histories](../../figures/2026-09-20/ember_two_sequences_hr.png)
shows the alternatives through 128.3 thousand years and 17.91 million years;
standalone flash diagrams remain separate from this draft.

Rebuild here with `tectonic --keep-logs ember_status_and_future.tex`.
