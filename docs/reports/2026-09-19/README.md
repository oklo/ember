# Ember working paper — September 19, 2026

[Read the PDF](ember_status_and_future.pdf).

The diagnostic late helium-3 calculation reaches **201.2 thousand years**,
**6115 K** and nuclear power **1.811e32 erg/s** on **2115** mass points.
Nuclear burning supplies **36.84%** of photon luminosity. The star is still
expanding. The final **20,000 yr** use **15.40 CPU minutes**, including
timestep comparisons; all original accuracy and conservation checks pass.
A matched **2000 yr** mass-grid comparison changes deposited energy by
**0.5220%**. The complete pulse remains uncertain.

The calculation resolving the onset reaches **9736 yr**, **4583 K** and
nuclear power **1.149e34 erg/s**. It shows a rise, decline and renewed burning.
The later sequence and this onset calculation have not yet established one
continuous, spatially converged flash followed by cooling.

The draft describes microscopic species transport down to **1 MK** in a
newly radiative envelope layer, supported by source ionization, actual flux
and **400/200/100 yr** comparisons. The local mixed-atmosphere grid has
**18 columns** over **6000–6500 K** and **log g = 4.35–4.65**. A matched
**9600 yr** comparison of two inferred hydrogen references changes endpoint
temperature by **0.04684 K**, with the same convection. This does not establish
absolute atmosphere accuracy.

The pre-flash main trajectory accounts for **15.13 elapsed hours** and
**36.31 CPU-hours**, excluding source construction and inter-run pauses.
The runtime of a complete full-physics cooling track is not yet measured.

This directory contains the PDF, LaTeX inputs and included figure PDFs.
Existing figures, including the LBA97 opacity comparison, are retained.
Standalone pulse diagrams are separate from this draft.

Rebuild from this directory with `tectonic --keep-logs ember_status_and_future.tex`.
