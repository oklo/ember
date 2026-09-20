# Ember working paper — September 19, 2026

[Read the PDF](ember_status_and_future.pdf).

The diagnostic late helium-3 calculation reaches **171.2 thousand years**,
**6079 K** and nuclear power **2.085e32 erg/s** on **2115** mass points.
Nuclear burning supplies **51.93%** of photon luminosity. A matched
**2000 yr** comparison of two mass grids changes deposited nuclear energy by
**0.5220%**. A checked **800 yr** trial-step limit speeds the smooth later
sequence; the **17,200 yr** segment uses **31.21 CPU minutes**, including
its timestep comparisons and rejected trials. The final **16,800 yr**
use **22.69 CPU minutes**, with unchanged accuracy checks.

The calculation resolving the onset reaches **6736 yr** and **4594 K**.
It shows a rise, decline and renewed burning. Its latest **1.147 yr**
convective adjustment passes a **2/4-step** comparison, with deposited
energy differing by **0.1756%**. The full pulse's peak and spatial
convergence remain uncertain.

The draft includes the checked diffusion extension at the newly radiative
envelope layer and a **27-column** local mixed-atmosphere grid covering
**5500–6500 K** and **log g = 4.50–4.80**. A complete white-dwarf cooling
sequence is not yet established.

A matched **3000 yr** atmosphere comparison changes the endpoint temperature by
**0.1460 K**, with the same convection. The nominal segment takes
**8.605 CPU minutes**, including time-step checks.

This directory contains the PDF, LaTeX inputs and included figure PDFs.
The existing figures, including the LBA97 opacity comparison, are retained.
Standalone pulse diagrams are separate from this draft.

Rebuild from this directory with `tectonic --keep-logs ember_status_and_future.tex`.
