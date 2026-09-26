# Ember working paper — September 26, 2026

[Read the PDF](ember_status_and_future.pdf).

The contraction calculation reaches **1 Gyr**, **2765 K**, **0.1264 solar radii**, and **99.23%** nuclear support. The new [PMS figure](pms_main_sequence.pdf) contains all **3882** saved model rows, including the starting model; its [data](pms_main_sequence.csv) and [plotting script](plot_pms_main_sequence.py) are included. This is a homogeneous initial-D/pp calculation, separate from the plotted long-term trajectories.

The paper describes the common burning/mixing/energy solver and the smooth composition-dependent EOS extension. All **512** main-sequence reference zones support composition derivatives. A bounded **1-million-year** finite-mixing and material-heat control converges on the actual 0.1-solar-mass structure, with a **3.856e-8** full-step/two-half-step difference in logarithmic structure variables. That control omits microscopic drift and does not extend the production track.

The long-term endpoints are unchanged: **4.002 Tyr** before the atmosphere adjustment; detailed flash/cooling history **5.744 Myr / 5903 K**, alternative history **42.48 Myr / 4771 K**, measured from their common pre-flash reference. These are separate histories. The late helium-3 pulse is not yet established in a continuous calculation from contraction.

The LBA97 opacity-map comparison and lifetime timeline are retained. Some cooler enriched EOS states and hotter states remain unsupported; initial-D transport and consistent atmosphere coverage remain unfinished. The complete PMS-to-white-dwarf calculation is not yet available.

Rebuild with `tectonic --keep-logs ember_status_and_future.tex`.
