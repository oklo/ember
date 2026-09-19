# Ember working paper — September 19, 2026

[Read the PDF](ember_status_and_future.pdf).

The diagnostic helium-3 pulse calculation with **2010** mass points reaches
**92.27 thousand years**, **5453 K** and surface **X = 0.8158**. Nuclear
power falls **23.4%** during the latest **1.640e4 yr**, to **4.396e32 erg/s**;
shell burning continues. All **41** intervals in this segment pass independent checks.
A matched local refinement changes deposited nuclear energy by **1.468%**
and endpoint nuclear power by **1.980%** over **978.3 yr**. This supports
the late shell burning; onset and peak remain uncertain.
A tested larger timestep ceiling uses **81** versus **120** structure solves
over **2000 yr**, with released nuclear energy differing by **0.1353%**.
A matched **10,000 yr** atmosphere comparison inside the lower-gravity interval
changes temperature by **0.2991 K** and luminosity by **0.02206%**, with the same
convective regions. The nominal segment uses **14.93 CPU minutes**.
The draft also describes a checked **1431 yr** onset calculation and the
**6000 K** atmosphere sources and independent depth comparisons.

This directory contains the PDF, LaTeX inputs and included figure PDFs.
The existing figures, including the LBA97 opacity comparison, are retained.
Standalone pulse diagrams are separate from this draft.

Rebuild from this directory with `tectonic --keep-logs ember_status_and_future.tex`.
