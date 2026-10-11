# Hydrogen burning near the stellar boundary

Ember tracks compare clear and cloudy atmospheres at five initial metal
abundances, from one third to three times Z_ref=0.01648, with Y=0.27.
The cloudy case uses f_sed=2 and clears between 900 and 1200 K. All tracks
begin on Hayashi contraction with deuterium.

![Transition mass and metallicity](transition_mass.png)

Bars span tested cooling and burning masses. Cooling requires nuclear
support below 0.1%; sustained burning requires 99–101% support at 1 Tyr,
with less than 2 percentage points of variation over the final 50 Gyr.
These are finite-age numerical brackets, not statistical confidence
intervals. Hollow diamonds are unfinished trials. Asterisks flag endpoint
pairs with different EOS coverage or envelope resolution; matched controls
are pending. Atmospheric and mixing uncertainties are additional.

![Cloudy tracks at three times the reference metal abundance](boundary.png)

The detailed figure contains 22 masses at Z=0.04945. Filled endpoints mark
1 Tyr or the fading cutoff; open endpoints are incomplete. The lower panel
shows luminosity colored by effective temperature. Lines connect saved
models without smoothing or extrapolated tails.

## Is the nuclear surplus physical?

![Mixing and metallicity comparisons](mixing_metallicity.png)

At 0.05825 solar masses, the maximum nuclear-to-surface luminosity ratio is
1.027 without semiconvection, 1.009 with weak composition mixing, and 1.004
so far with coupled heat/composition transport. The coupled calculation
assumes layer height H/d=25 and marginal transport between mixing regimes;
its effect depends on that uncertain assumption. Small abrupt excursions
can accompany individual cells joining convection. These results do not
establish a physical surge.

The nearby-metallicity tests hold mass at 0.05815 solar masses and helium
fixed. They stopped at different ages, so their largest saved ratios are
not established final peaks. The 3.03-Z_ref case uses a small declared
atmosphere continuation beyond the native grid.

A separating conductive layer allows composition gradients to develop.
Later, its stability depends on the helium-rich interior beneath the
hydrogen-rich envelope. Such a layer can undergo oscillatory mixing;
conduction does not imply complete isolation.
[Mirouh et al.](https://arxiv.org/abs/1112.4819),
[Moll et al.](https://arxiv.org/abs/1506.07900).
[Interior diagram](Ember_boundary_interior.png) and its 744 saved profiles
show the instantaneous-convection calculation, without semiconvection.

Useful numerical checks: the 0.057 solar-mass burning maximum agrees to
0.0001 in matched 50- and 25-Myr tests. The 0.05825 surplus survives halving
the timestep and resolving more of the envelope; near 129 Gyr its 1.877%
excess agrees with the change in internal plus gravitational energy.
However, the 0.058 upturn disappears with shorter timesteps.
[Timestep comparison](timestep_comparison.png).

Remaining uncertainties include atmospheric H/He response, prescribed cloud
clearing and the additive ideal-metal EOS. Finer mass spacing does not
remove them. Atmosphere-table irregularities also affect some early
shoulders; see [Morley et al.](https://arxiv.org/html/2402.00758v1).

## Reproduce the figures

Run `python plot.py`, `python plot_transition_mass.py`,
`python plot_mixing_metallicity.py` and `python plot_interior.py` with NumPy
and Matplotlib. CSV/NPZ files retain plotted states; JSON records source
fingerprints, endpoints and assumptions. These reproduce figures, not
stellar runs: the full physics tables and checkpoints are not included.
