# Ember working paper — September 20, 2026

[Read the PDF](ember_status_and_future.pdf).

The diagnostic late helium-3 calculation reaches **1.171 million years**,
**6398 K** and nuclear power **6.418e31 erg/s** on **2115** mass points.
Nuclear burning supplies **38.39%** of photon luminosity. Radius and
surface luminosity have passed maxima and declined by **48.37%** and
**68.75%**; effective temperature passes a shallow maximum of **6399 K** near
**1.088 million years**, then falls **1.151 K** over **83,200 years**.
Both atmosphere references retain this local turn toward cooling. The full
connection through the flash and a long white-dwarf cooling sequence remain
unverified.

The calculation resolving the onset reaches **13.13 thousand years**,
**4741 K** and nuclear power **5.285e34 erg/s** on **1899** mass points. The 1695-point calculation shows a local maximum
of **3.719e35 erg/s**, an **81.64%** decline and renewed burning. The overall
maximum, complete pulse and connection to the later calculation remain unresolved.
A checked **0.6656-year** mixing adjustment changes surface hydrogen from
**X = 0.9984** to **X = 0.8085** and helium-3 from **X3 = 0.001368** to
**X3 = 0.007264**. Two-step/four-step paths pass the original physical checks.

Over the final **296.7 years** of the onset calculation, nuclear power falls
**42.09%** and deposited nuclear energy is **6.334e44 erg**. All **80** retained
intervals pass the unchanged timestep and conservation checks.

A matched **100-year** comparison on **1763/1899 points** changes deposited
nuclear energy by **4.026%** and endpoint nuclear power by **5.589%**.
Surface temperatures differ by **0.03544 K**. The finer mesh spreads half
the burning over **31 cells** and is used for continued calculation with
explicit spatial uncertainty; this does not establish convergence of the full pulse.

A **3200/1600/800 yr** comparison tests time resolution. The final
**128,000 years** use trial steps up to **12,800 years**, with every interval
passing the same accuracy and conservation criteria. This segment uses
**21.56 CPU minutes**, including timestep comparisons. The higher-gravity
atmosphere interval passes a separate **3200/1600 yr** comparison.

The mixed-atmosphere grid has **42 columns** over **6000–6500 K** and
**log g = 4.20–5.10**, with corresponding lower-density EOS coverage.
The absolute hydrogen atmosphere below the published gravity range remains
inferred. Microscopic species transport extends locally to **1 MK**, supported
by ionization and finite-step checks; this is not neutral-fluid transport.

A matched **25,600-year** comparison in the **log g = 4.80–4.95** interval
changes effective temperature by **0.01540 K**, surface luminosity by
**0.0004895%**, and integrated nuclear heat by **2.359e-6%**, with identical
final convection. The preceding **38,400 years** are reused because both
references have identical table values there. This tests local sensitivity,
while the absolute hydrogen reference remains inferred.
A **76,800-year** comparison in the next gravity interval changes temperature
by **0.1187 K**, surface luminosity by **0.002568%**, and nuclear heat by
**3.438e-5%**, with identical convection. Its **51,200-year** common prefix
is reused; both prescriptions retain the temperature decline.

The retained pre-flash trajectory accounts for **15.13 elapsed hours** and
**36.31 CPU-hours**, excluding source construction and inter-run pauses.
A complete full-physics cooling runtime has not yet been measured.

This directory contains the PDF, LaTeX inputs and included figure PDFs.
Existing figures, including the LBA97 opacity comparison, are retained.
Standalone pulse diagrams remain separate from this draft.

Rebuild here with `tectonic --keep-logs ember_status_and_future.tex`.
