# Ember working paper — September 26, 2026

[Read the PDF](ember_status_and_future.pdf).

The common evolution program carries the star from a Hayashi start to
**1 Gyr**, **2765 K**, **0.1264 solar radii**, and **99.17%** nuclear support.
The [PMS figure](pms_main_sequence.pdf) contains **1020** saved rows, including
the starting model; its [data](pms_main_sequence.csv) and
[plotting script](plot_pms_main_sequence.py) are included. Initial D/pp/CN,
plasma losses and checked whole-star convection share one history.
Composition heat is included after initial D exhaustion. This sequence
remains separate from the plotted late-evolution histories.

The same star has continued to **4 Gyr**, with sustained hydrogen burning
and no rejected steps. The figure emphasizes contraction through 1 Gyr.

The paper describes the common burning/mixing/energy solver and the smooth composition-dependent EOS extension. All **512** main-sequence reference zones support composition derivatives. A bounded **1-million-year** finite-mixing and material-heat control converges on the actual 0.1-solar-mass structure, with a **3.856e-8** full-step/two-half-step difference in logarithmic structure variables. That control omits microscopic drift and does not extend the production track.

All **1019** accepted intervals pass isotope and energy checks, with no
rejected intervals, using **8.166 CPU minutes** including material setup.
Exact restart replay passes. The paper explains the convective approximation
and stellar sensitivity to uncertain kinetic heat and conduction. Radiative
microscopic transport and continuous atmosphere coverage remain unfinished.

The long-term endpoints are unchanged: **4.002 Tyr** before the atmosphere adjustment; detailed flash/cooling history **5.744 Myr / 5903 K**, alternative history **42.48 Myr / 4771 K**, measured from their common pre-flash reference. These are separate histories. The late helium-3 pulse is not yet established in a continuous calculation from contraction.

The LBA97 opacity-map comparison and lifetime timeline are retained. Some cooler enriched EOS states and hotter states remain unsupported; cool microscopic transport and consistent atmosphere coverage remain unfinished. The complete PMS-to-white-dwarf calculation is not yet available.

Rebuild with `tectonic --keep-logs ember_status_and_future.tex`.
