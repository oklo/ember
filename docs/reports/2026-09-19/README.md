# Ember working paper — September 19, 2026

[Read the PDF](ember_status_and_future.pdf).

The diagnostic helium-3 pulse calculation with **2010** mass points reaches
**116.7 thousand years**, **5753 K** and surface **X = 0.8158**. Nuclear
power falls **4.729%** during the latest **4400 yr**, to **3.238e32 erg/s**;
shell burning continues. All **11** intervals in this segment pass independent checks.
A matched local refinement changes deposited nuclear energy by **1.468%**
and endpoint nuclear power by **1.980%** over **978.3 yr**. This supports
the late shell burning; onset and peak remain uncertain.
A tested larger timestep ceiling uses **81** versus **120** structure solves
over **2000 yr**, with released nuclear energy differing by **0.1353%**.
A matched **12,000 yr** atmosphere comparison reaches **5701 K** and changes
surface temperature by **0.2296 K** and luminosity by **0.01569%**, with the same
convective regions. The final **4400 yr** segment uses **11.27 CPU minutes**.
The draft also describes a checked **1619 yr** onset calculation. A **3.573 yr**
convective adjustment passes an eight-step/sixteen-step comparison, with nuclear
energy differing by **0.1222%** and the same final convection.
Atmosphere source coverage reaches **6000 K** and **log g = 4.80**, with independent
depth comparisons. A continuous, converged passage through the whole pulse remains unresolved.

This directory contains the PDF, LaTeX inputs and included figure PDFs.
The existing figures, including the LBA97 opacity comparison, are retained.
Standalone pulse diagrams are separate from this draft.

Rebuild from this directory with `tectonic --keep-logs ember_status_and_future.tex`.
