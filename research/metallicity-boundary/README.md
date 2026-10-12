# Hydrogen burning near the stellar boundary

Ember compares clear and cloudy atmospheres from one third to three times
Z_ref = 0.01648, with Y = 0.27 and solar metal ratios. The cloudy models use
f_sed = 2 and clear between 900 and 1200 K. Tracks start on Hayashi contraction
with deuterium.

![Transition mass and metallicity](transition_mass.png)

Bars span tested cooling and burning masses with metal ionization and
matching physical inputs:

| Initial Z | Clear (solar masses) | Cloudy (solar masses) |
|---|---|---|
| 0.005494 | 0.075–0.07525 | 0.07025–0.0705 |
| 0.009516 | 0.0735–0.074 | 0.068–0.069 |
| 0.01648 | 0.07125–0.0725 | 0.0665–0.06675 |
| 0.02855 | 0.069–0.070 | 0.0635–0.0645 |
| 0.04945 | 0.066–0.0665 | 0.060–0.061 |

Squares have 99–101% nuclear support at 1 Tyr, varying by less than two
percentage points over the final 50 Gyr. Open circles have declined below
0.1%; diamonds are unfinished trials. These intervals measure mass sampling,
not physical uncertainty or permanent main-sequence membership. Shading
connects the five sampled metallicities. The
[neutral-metal EOS comparison](transition_mass_comparison.png) is separate;
different EOS treatments must not be combined into a tighter interval.

The high-Z clear cooling endpoint extends the permitted atmospheric H/He
range slightly. Its composition response remains provisional. Local atmosphere
calculations include corrected electron-continuum terms and independent
flux checks; these checks do not validate a boundary with a different mixture.
The imported Bobcat, Diamondback, SPHINX and MARCS tables were not generated
by this local calculation path. New continuations do not tighten the plotted
mass intervals. Cloud physics, composition response, mixing and dense-matter
thermodynamics remain larger uncertainties than the mass spacing in several
intervals.

Work above three times Z_ref remains experimental and is not included in
the transition-mass diagram. The [five-times-reference-metallicity controls](experimental_highz/Ember_boundary_5Z_experimental.png)
show a 0.06-solar-mass clear-atmosphere model with near-complete nuclear
support through 891.4 Gyr. [Assumptions and data](experimental_highz/README.md).
EOS tables have been assembled through Z = 0.24, but unresolved source
transitions prevent general use. Default table limits are unchanged;
no higher-metallicity hydrogen-burning boundary is established.

The [current cloudy three-times-reference-metallicity comparison](current_cloudy_3z/Ember_boundary_3Z_cloudy_current.png)
shows 0.060, 0.0605 and 0.061 solar masses with metal ionization. The middle
mass uses a wider EOS domain; its dashed late continuation awaits an atmosphere
check. [Assumptions and data](current_cloudy_3z/README.md).

## Neutral-metal comparisons

The following figures use an EOS that omits metal electrons from degenerate
pressure. They illustrate transport sensitivity, rather than the corrected
mass intervals above.

- [Cloudy tracks at three times Z_ref](boundary.png), and the
  [selected curves](boundary_selected.png).
- [Structure of the 0.05835-solar-mass model](interior_05835/Ember_boundary_05835_structure.png).
  A conductive layer separates the core and envelope. Burning fades near
  880 Gyr despite about 98.24% of the initial hydrogen remaining at 897 Gyr.
  The omitted metal pressure and late heat-capacity interpolation make this
  outcome provisional. [Physical interpretation and data](interior_05835/README.md).
- [Mixing comparisons](mixing_metallicity.png) and
  [independent meshes](mesh_comparison.png). Sharp mixing events change with
  timestep and mesh. Coupled heat/composition mixing removes the broad peak
  and dip on both tested meshes; both retain over 99.8% nuclear support at
  1 Tyr, with a 4.945% luminosity difference. The transport prescription
  remains uncertain.
- [Timestep comparison](timestep_comparison.png): the separate
  0.058-solar-mass upturn disappears with shorter steps.

## Reproduce the figures

`plot_transition_mass.py` reproduces the cutoff diagram from
`transition_mass.json`; add `--include-comparison` to regenerate the separate
neutral-metal comparison. Other figures use `plot.py`,
`plot_mixing_metallicity.py`, `plot_mesh_comparison.py`, `plot_interior.py`
and `interior_05835/plot_structure.py`, with NumPy and Matplotlib.
CSV/NPZ files contain plotted states; JSON records endpoints, assumptions,
input fingerprints and numerical checks. These reproduce figures; the full
stellar checkpoints and physics tables are separate.

For the selected curves:
`python plot.py --omit-masses .05815 .05820 .05825 .05830 --output-stem boundary_selected`.
