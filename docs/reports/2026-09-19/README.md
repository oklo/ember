# Ember working paper — September 19, 2026

[Read the PDF](ember_status_and_future.pdf).

The diagnostic late helium-3 calculation reaches **251.2 thousand years**,
**6150 K** and nuclear power **1.552e32 erg/s** on **2115** mass points.
Nuclear burning supplies **29.01%** of photon luminosity. The star is still
expanding. The final **50,000 yr** use **33.63 CPU minutes**, including
timestep comparisons; all original accuracy and conservation checks pass.

The calculation resolving the onset reaches **12.50 thousand years**,
**4557 K** and nuclear power **1.213e35 erg/s**. It has passed a local maximum
of **3.719e35 erg/s**, an **81.64%** decline and renewed burning. The current
rise is not its first. The overall maximum, complete pulse and connection
to the later calculation remain unresolved.

A **12.26 yr** convection merger passes a **16/32-step** comparison;
integrated nuclear energies differ by **0.02667%**. A subsequent **2.587 yr**
adjustment passes a **4/8-step** comparison, differing by **0.02665%**.
These support the endpoint structures, without resolving every transient.

The mixed-atmosphere grid has **18 columns** over **6000–6500 K** and
**log g = 4.20–4.50**, with the corresponding lower-density EOS coverage.
The absolute hydrogen atmosphere below the published gravity range remains
inferred. Microscopic species transport extends locally to **1 MK**, supported
by ionization and finite-step checks; this is not neutral-fluid transport.

The retained pre-flash trajectory accounts for **15.13 elapsed hours** and
**36.31 CPU-hours**, excluding source construction and inter-run pauses.
A complete full-physics cooling runtime has not yet been measured.

This directory contains the PDF, LaTeX inputs and included figure PDFs.
Existing figures, including the LBA97 opacity comparison, are retained.
Standalone pulse diagrams are separate from this draft.

Rebuild here with `tectonic --keep-logs ember_status_and_future.tex`.
