# Ember working paper — September 20, 2026

[Read the PDF](ember_status_and_future.pdf).

The diagnostic late helium-3 calculation reaches **2.254 million years**,
**6301 K** and nuclear power **3.143e31 erg/s** on **2115** mass points.
Nuclear burning supplies **35.66%** of photon luminosity. Radius and
surface luminosity have passed maxima and declined by **61.36%** and
**83.53%**; effective temperature passes a shallow maximum of **6399 K** near
**1.088 million years**, then falls **97.80 K** over **1.166 million years**.
Both atmosphere references retain the initial turn toward cooling; the nominal calculation continues through the full decline quoted here. The full
connection through the flash and a long white-dwarf cooling sequence remain
unverified.

The calculation resolving the onset reaches **14.77 thousand years**,
**4765 K** and nuclear power **4.177e34 erg/s** on **1899** mass points. The 1695-point calculation shows a local maximum
of **3.719e35 erg/s**, an **81.64%** decline and renewed burning. The overall
maximum, complete pulse and connection to the later calculation remain unresolved.
A checked **0.6656-year** mixing adjustment changes surface hydrogen from
**X = 0.9984** to **X = 0.8085** and helium-3 from **X3 = 0.001368** to
**X3 = 0.007264**. Two-step/four-step paths pass the original physical checks.

Over the final **707.6 years** of the onset calculation, nuclear power rises
**16.54%** and deposited nuclear energy is **8.515e44 erg**. All **99** retained
intervals pass the unchanged timestep and conservation checks.

A matched **100-year** comparison on **1763/1899 points** changes deposited
nuclear energy by **4.026%** and endpoint nuclear power by **5.589%**.
Surface temperatures differ by **0.03544 K**. The finer mesh spreads half
the burning over **31 cells** and is used for continued calculation with
explicit spatial uncertainty; this does not establish convergence of the full pulse.

A **3200/1600/800 yr** comparison tests time resolution. The final
**178,100 years** use trial steps up to **25,600 years**, with every interval
passing the same accuracy and conservation criteria. This segment uses
**18.40 CPU minutes**, including timestep comparisons. The higher-gravity
atmosphere interval passes a separate **3200/1600 yr** comparison.

The mixed-atmosphere grid has **48 columns** over **6000–6500 K** and
**log g = 4.20–5.25**, with corresponding lower-density EOS coverage.
The extension through **log g = 5.25** passes source, depth, native and matched stellar checks. A **102,400-year** alternative-reference comparison changes temperature by **0.07474 K**, surface luminosity by **0.001740%**, and nuclear heat by **2.291e-5%**, with identical final convection.
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
The added [late cooling figure](late_cooling.pdf) shows the checked temperature turn and luminosity decline. Its [data](late_cooling.csv) and plotting script are included. The LBA97 opacity comparison is retained.
Standalone pulse diagrams remain separate from this draft.

Rebuild here with `tectonic --keep-logs ember_status_and_future.tex`.
