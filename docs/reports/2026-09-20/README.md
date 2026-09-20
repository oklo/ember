# Ember working paper — September 20, 2026

[Read the PDF](ember_status_and_future.pdf).

The diagnostic late helium-3 calculation reaches **749.6 thousand years**,
**6372 K** and nuclear power **9.184e31 erg/s** on **2115** mass points.
Nuclear burning supplies **35.16%** of photon luminosity. Radius and
surface luminosity have passed maxima and declined by **34.95%** and
**51.19%**; effective temperature still rises. The star is contracting,
without an established white-dwarf cooling sequence.

The calculation resolving the onset reaches **12.65 thousand years**,
**4721 K** and nuclear power **1.571e35 erg/s**. It has passed a local maximum
of **3.719e35 erg/s**, an **81.64%** decline and renewed burning. The overall
maximum, complete pulse and connection to the later calculation remain unresolved.
A checked **0.6656-year** mixing adjustment changes surface hydrogen from
**X = 0.9984** to **X = 0.8085** and helium-3 from **X3 = 0.001368** to
**X3 = 0.007264**. Two-step/four-step paths pass the original physical checks.

A **3200/1600/800 yr** comparison supports the late timestep ceiling with all
original physical and conservation criteria. A **170,400 yr** segment uses
**60.49 CPU minutes**, including timestep comparisons. The higher-gravity
atmosphere interval passes a separate **3200/1600 yr** comparison.

The mixed-atmosphere grid has **30 columns** over **6000–6500 K** and
**log g = 4.20–4.80**, with corresponding lower-density EOS coverage.
The absolute hydrogen atmosphere below the published gravity range remains
inferred. Microscopic species transport extends locally to **1 MK**, supported
by ionization and finite-step checks; this is not neutral-fluid transport.

The retained pre-flash trajectory accounts for **15.13 elapsed hours** and
**36.31 CPU-hours**, excluding source construction and inter-run pauses.
A complete full-physics cooling runtime has not yet been measured.

This directory contains the PDF, LaTeX inputs and included figure PDFs.
Existing figures, including the LBA97 opacity comparison, are retained.
Standalone pulse diagrams remain separate from this draft.

Rebuild here with `tectonic --keep-logs ember_status_and_future.tex`.
