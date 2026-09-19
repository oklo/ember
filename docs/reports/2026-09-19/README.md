# Ember working paper — September 19, 2026

[Read the PDF](ember_status_and_future.pdf).

The diagnostic helium-3 pulse calculation with **2010** mass points reaches
**100.3 thousand years**, **5553 K** and surface **X = 0.8158**. Nuclear
power falls **10.49%** during the latest **8000 yr**, to **3.935e32 erg/s**;
shell burning continues. All **20** intervals in this segment pass independent checks.
A matched local refinement changes deposited nuclear energy by **1.468%**
and endpoint nuclear power by **1.980%** over **978.3 yr**. This supports
the late shell burning; onset and peak remain uncertain.
A tested larger timestep ceiling uses **81** versus **120** structure solves
over **2000 yr**, with released nuclear energy differing by **0.1353%**.
A matched **8000 yr** atmosphere comparison entering the hotter interval
changes temperature by **4.301 K** and luminosity by **0.3070%**, with the same
convective regions. The nominal segment uses **10.89 CPU minutes**.
The draft also describes a checked **1569 yr** onset calculation: a **2.548 yr**
comparison with **16** and **32** steps agrees in nuclear energy to **0.01692%**.
Atmosphere source coverage reaches **6000 K** and **log g = 4.80**, with independent
depth comparisons. A continuous, converged passage through the whole pulse remains unresolved.

This directory contains the PDF, LaTeX inputs and included figure PDFs.
The existing figures, including the LBA97 opacity comparison, are retained.
Standalone pulse diagrams are separate from this draft.

Rebuild from this directory with `tectonic --keep-logs ember_status_and_future.tex`.
