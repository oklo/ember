# Interior evolution at 0.05835 solar masses

The retained cloudy calculation has initial Z = 0.04945, Y = 0.27 and
f_sed = 2. The diagram follows it to 897 Gyr, when nuclear support falls
below 0.1%. It does not extrapolate to 1 Tyr.

- Top: nuclear luminosity divided by surface luminosity, with central
  temperature on the right axis. Nuclear support peaks at 99.72%.
- Middle: convection and the dominant diffusive heat carrier outside
  convective regions. Dashed contours give temperature. Electron conduction
  exceeding photon diffusion does not imply that it exceeds convective heat
  transport in a convective region.
- Bottom: hydrogen mass fraction, dotted metal-abundance contours and white
  convective boundaries. The two convective regions do not mix with one another.

The right column enlarges the final decline. Vertical dotted lines mark
the prescribed cloud clearing between 1200 and 900 K. They do not establish
cloud clearing as the cause: central cooling starts much earlier.

At the endpoint, the center is at 0.6103 MK, the surface at 825.4 K, and
the outer 12.20% of the mass remains convective. About 98.24% of the initial
hydrogen remains. Sampled nuclear-rate decompositions show that almost all
heating comes from pp reactions producing helium-3, rather than its destruction.

The calculation omits semiconvection and uses additive neutral metals in
the core EOS. Pointwise comparisons identify missing metal-electron pressure
and a late heat-capacity interpolation error. The pressure correction grows
as the core accumulates metals. Correcting these issues can change the
trajectory; the precise mass and failure age are not validated predictions.

Run `python plot_structure.py` with NumPy and Matplotlib to recreate the PNG
and PDF. The NPZ contains selected calculated profiles, including both sides
of convective-boundary changes; the JSONL contains the retained global track.
`provenance.json` records source fingerprints and comparison checks for the
profile reconstruction. The atmosphere is not spatially resolved here.
