# Does the current F77 calculation reproduce the 1997 paper?

The broad sequence agrees, but the quantitative agreement is poor for the
0.1-solar-mass star. The local F77 program is a recent reconstruction, with
documented substitutions. It is not the archived executable used by Laughlin,
Bodenheimer and Adams in [The End of the Main Sequence](https://doi.org/10.1086/304125).

| Quantity | LBA97, published | Current F77, AJR opacity | Current F77, Ferguson opacity |
| --- | ---: | ---: | ---: |
| Age of maximum helium-3 abundance, trillion years | 1.38 | 0.4153 | 0.3799 |
| Maximum helium-3 mass fraction | 0.0995 | 0.07964 | 0.07743 |
| Age when the center becomes radiative, trillion years | 5.742 | 3.224 | 2.946 |
| Central hydrogen mass fraction then | about 0.16 | 0.01756 | 1.064e-5 |
| Highest surface temperature, K | 5807 | 5936 | 6282 |

Published values come from section 3.2, pages 422–424, and Figure 1. The
helium-3 maximum in the text is 9.95%; the figure labels it 9.96%. Several prose
ages on page 423 are printed as Gyr, but the figure and surrounding discussion
establish the trillion-year scale. The radiative-core label reads 5742 Gyr.

The current AJR star forms its central radiative region about 43.9% earlier,
with much less hydrogen remaining. The peak temperature is within about 2.2%,
but that single agreement does not validate the complete track. The change in
composition at core formation means a uniform rescaling of time cannot bring
the calculations into agreement.

The reconstruction's high-temperature opacity uses the Cox–Tabor King IVa
table with composition dependence supplied by the Bodenheimer et al. formula.
The published model instead interpolates between the King IVa and Ross–Aller 2
mixtures using the Weiss, Keady and Magee formulations. The reconstruction
also uses Caughlan–Fowler nuclear rates in place of the cited Bahcall rates.
Selecting Ferguson introduces another change in low-temperature opacity.
These are differences identified in the source, not an established explanation
of the discrepancy. No controlled experiment has yet assigned their separate
effects on the published milestones.

F77's central convection flag changes only after a finite margin beyond
marginal stability. The [recorded transition brackets](results/f77_core_onset_comparison_v1.json)
retain that distinction. The Ferguson calculation later fails during cooling;
its last accepted state does not establish a completed cooling calculation.
No new F77 evolution or numerical convergence study was needed for this
comparison. The [comparison data](results/f77_lba97_comparison_v1.json) were
extracted from the same saved calculations used for the working-paper graphs.

`scripts/summarize_f77_lba97.py` checks the published figure CSV against its
extraction record and checks the retained F77 source hash before selecting
milestones. It also records the original paper's PDF hash. The early luminosity
minimum included in that data file is explicitly distinguished from the paper's
definition of the start of the main sequence, which requires nuclear heating
to supply all of the surface luminosity.
