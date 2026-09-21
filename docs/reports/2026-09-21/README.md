# Ember working paper — September 21, 2026

[Read the PDF](ember_status_and_future.pdf).

The calculation with finite convective transport from the helium-3 flash onset
reaches **4.619 million years**, **6025 K** and nuclear power **1.371e31 erg/s**
on **1983** mass points. Nuclear power continues to decline through independently checked intervals. A matched
**100-year** comparison on **1927/1983** points differs by **8.261%** in nuclear
heat and **15.24%** in final nuclear power. The first year after refinement has
an integrated energy residual of **1.844e43 erg**, associated with a mechanical
readjustment. Correcting only the initial radii leaves the one-year nuclear
heat and final power unchanged within **1e-12** relative, but retains the
explicit initial energy change. The complete flash and its spatial accuracy
remain unresolved.

An alternative flash history reaches **42.48 million years**, **4771 K** and
nuclear power **3.166e30 erg/s** on **2115** mass points. Nuclear burning supplies
**55.85%** of photon luminosity. Radius and surface luminosity have declined
**82.91%** and **98.94%** from their maxima. Effective temperature passes a
maximum of **6399 K** near **1.088 million years**, then falls **1628 K** over
**41.39 million years**. The **24.48–29.48 million year** interval uses **74.91 CPU minutes** and **40.23 wall minutes**, including timestep comparisons and rejected trials. Nuclear power falls **7.869%**, while the whole-star helium-3 inventory grows **0.2052%**. At **29.48 million years**, helium-3 production is **1.263e12 g/s** and destruction is **4.106e11 g/s**. Helium-3 fusion supplies **23.12%** of nuclear power, down from **26.97%** five million years before. The discrete source sum reproduces the inventory change to better than **1e-10** relative. This establishes the balance of production and destruction; it does not establish stability against a later flash. The full cooling-track runtime remains unmeasured. Both histories now turn toward lower effective temperature. The finite-mixing calculation reaches **6422 K** near **1.139 million years** before declining.

A provisional pre-main-sequence calculation now follows **100,000 years** of contraction with its initial deuterium and helium-3 inventories intact. It uses **74.79 CPU seconds**, excluding atmosphere preparation, and ends at **2876 K** with radius **1.389 solar radii**. An independent atmosphere inside the starting cell differs from interpolation by **−1.304%** in matching temperature and **+2.156%** in pressure. A local stellar sensitivity test changes the initial effective temperature by **10.89 K** and radius by **0.5213%**. The tested trace-isotope effect is much smaller. Local grid refinement remains in progress; this is not yet the accepted physical reference track.


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

The selected cooling atmosphere window has **18** columns over **4750–5500 K** and **log g = 5.80–6.00** at three hydrogen fractions. Twelve occupied rows are preserved exactly; six new columns pass independent source and depth checks. The combined table passes interpolation, EOS and full stellar-response comparisons, followed by **6 million years** of reviewed cooling. Higher-gravity columns through **log g = 6.20** pass source and response checks; the missing **4500 K** coverage currently limits further cooling. A matched **million-year** reference comparison gives **0.2632 K** and **0.0002222%** in nuclear heat, with the same final convection. Independent **5250 K** and **5750 K** columns test interpolation within the broader source collection; at **5250 K**, the composition response agrees to **0.000007587%** in temperature and **0.1582%** in pressure.

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

A further **25,000-year** comparison below **log g = 4.50** gives **0.4703 K** difference in endpoint temperature and **0.006257%** in nuclear heat. The final convective regions agree. The selected extension uses the measured gravity response of solved near-hydrogen atmospheres, preserving the occupied reference values. This checks local sensitivity, not the absolute atmosphere.

Microscopic species transport extends to **1.000 MK** under a locally checked ionization approximation. A **1600/800/400-year** comparison gives **0.2686%** difference in nuclear heat with unchanged accuracy and conservation criteria. The heat-transport blend remains at **2–3 MK**.

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
shows the alternatives through 4.619 million years and 42.48 million years;
standalone flash diagrams remain separate from this draft.

Rebuild here with `tectonic --keep-logs ember_status_and_future.tex`.
