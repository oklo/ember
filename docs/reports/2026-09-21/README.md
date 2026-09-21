# Ember working paper — September 21, 2026

[Read the PDF](ember_status_and_future.pdf).

The calculation with finite convective transport from the helium-3 flash onset
reaches **201.1 thousand years**, **6108 K** and nuclear power **1.436e32 erg/s**
on **1983** mass points. Nuclear power continues to decline through independently checked intervals. A matched
**100-year** comparison on **1927/1983** points differs by **8.261%** in nuclear
heat and **15.24%** in final nuclear power. The first year after refinement has
an integrated energy residual of **1.844e43 erg**, associated with a mechanical
readjustment. Correcting only the initial radii leaves the one-year nuclear
heat and final power unchanged within **1e-12** relative, but retains the
explicit initial energy change. The complete flash and its spatial accuracy
remain unresolved.

An alternative flash history reaches **23.48 million years**, **5023 K** and
nuclear power **3.981e30 erg/s** on **2115** mass points. Nuclear burning supplies
**44.44%** of photon luminosity. Radius and surface luminosity have declined
**80.61%** and **98.33%** from their maxima. Effective temperature passes a
maximum of **6399 K** near **1.088 million years**, then falls **1376 K** over
**22.39 million years**. The final **2.171 million years** use **31.46 CPU minutes** in 63 completed native calls, including timestep comparisons and rejected trials. A preceding **32,400-year** check of stable microscopic arithmetic and solved species-flux conservation differs by **0.02234%** in nuclear heat, with unchanged physical tolerances and convection. These are local timings; the full cooling-track runtime remains unmeasured.

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

The selected cooler atmosphere grid has **24** source columns over **5000–5500 K** and **log g = 5.40–6.00**, with six new depth-checked columns and 18 preserved rows. A matched **million-year** reference comparison gives **0.2632 K** and **0.0002222%** in nuclear heat, with the same final convection. Atmosphere coverage below **5000 K** is the next cooling requirement. Independent **5250 K** and **5750 K** columns test interpolation within the broader source collection. At **5250 K**, the composition response agrees to **0.000007587%** in temperature and **0.1582%** in pressure.

A matched **2.1-million-year** cooling comparison over **5500–6000 K** changes
endpoint temperature by **2.292 K** and nuclear heat by **0.002680%** between
two hydrogen reference prescriptions. The detailed-onset comparison over
**10,000 years** changes the endpoint by **4.187 K** and nuclear heat by
**0.009269%**. Both comparisons retain identical final convective regions.
These local tests do not establish the absolute hydrogen boundary condition.

A matched **10,000-year** numerical comparison reduces native CPU time from
**46.28 to 27.65 minutes**, a **40.25%** saving, with unchanged accuracy
criteria. Temperature differs by **0.07287 K**, nuclear heat by **0.03629%**,
and the convective regions agree. A separate **640-year** timestep comparison
passes. These measurements support the selected detailed continuation;
they do not establish the runtime or accuracy of a complete cooling track.

A matched **7500-year** comparison entering the **6000–6500 K** atmosphere interval gives differences of **1.008 K**, **0.06316%** in surface luminosity and **0.00002583%** in nuclear heat, with identical final convection. Each accepted interval passes the original full-step/two-half-step checks.

A matched **25,000-year** comparison of the corrected hot reference tables gives **1.366 K**, **0.06345%** in surface luminosity and **0.01279%** in nuclear heat, with identical final convection. Adaptive full intervals may reach **2560 years**; every accepted interval retains the original checks.

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
The [standalone HR diagram of both histories](../../figures/2026-09-21/ember_two_sequences_hr.png)
shows the alternatives through 201.1 thousand years and 23.48 million years;
standalone flash diagrams remain separate from this draft.

Rebuild here with `tectonic --keep-logs ember_status_and_future.tex`.
