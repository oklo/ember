# Ember HR overview — October 5

`ember_hr.png` and `ember_hr.pdf` show 23 sequences across 0.005–0.235 solar masses. The comparison sequences reach 10 Gyr; the 0.1-solar-mass sequence continues into white-dwarf cooling. There are 4530 saved states in those sequences, including 850 on the retained 0.1-solar-mass track. An additional 11 states show a separate revised-envelope comparison.

Six complete cloud-free tracks have been added at 0.005, 0.010, 0.020, 0.030, 0.040 and 0.050 solar masses. Their cool-envelope assumptions and checks are described in `docs/research/fable/results/COOL_MASS_TRACKS_20261005.md`. Numerical completion is not observational validation. The provisional main-sequence boundary brackets are unchanged.

The retained cooling sequence reaches 373.3 K at 4.190 Tyr. The revised CMS hydrogen-envelope calculation reaches 353.0 K at 4.205 Tyr and is shown with separate outlined squares. Its initial EOS readjustment is omitted rather than presented as stellar evolution. Cold physics was introduced through continuations; this is not a single fixed-physics calculation from the Hayashi track.

Core crystallization marks equality of the central liquid and solid free energies, interpolated between saved models. It does not mean half the star is solid. Nuclear peak symbols require at least 1% depletion of the corresponding fuel; negligible mathematical rate maxima are omitted. Lithium is a passive diagnostic restricted to valid fully convective histories.

The hexagons now include 59 additional completed atmosphere-job receipts, adding 2.283 CPU-h. The displayed table-work estimate is 310.4 CPU-h. It excludes work without recoverable timing and is not the project's total compute. Hexagon placement uses the existing 0.1-solar-mass reference; it is an allocation of shared table work, not a map of interpolation coverage.

Run `out/science-env-sept18-v1/bin/python out/headline-hr-oct5-v1/make_figure.py` to render the saved inputs. `gather_update.py` reconstructs the selected histories and checks restart joins, complete age coverage and accepted-step audits. `update_compute.py` updates the receipt accounting. Inputs, source paths and checksums are recorded in the JSON companions. No private Jupiter trial or new no-fingering-flux cloudy comparison is included pending review.
