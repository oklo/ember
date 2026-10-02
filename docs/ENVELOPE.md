# Integrated outer envelope

The lifetime driver can integrate the outer layer between the non-gray atmosphere
at optical depth100 and a fixed enclosed mass. The same boundary calculation is
used during Hayashi contraction, hydrogen burning and cooling. Its EOS, opacity
and mixing-length transport must cover every state queried.

Select the interior EOS throughout the envelope with:

```text
envelope_eos "interior"
envelope_mass_fraction "0.0015"
```

The mesh ends below this layer. Its mass and composition remain in the outer
cell's mixing and energy reservoir; checkpoint version7 records that mass.
The reported radius and effective temperature refer to the photosphere.
Older checkpoints without an envelope retain their original interpretation.

Luminosity and composition are uniform through the integrated layer. Its thermal
energy is represented by the base state, so this approximation needs assessment
when the envelope supplies appreciable luminosity. The driver requires a connected
convective base. Native interior-EOS envelopes currently allow base temperatures
from 10 kK to 2 MK. By default, the absolute reservoir energy change plus pressure
work must stay below 0.1% of luminosity per interval, including below 1 MK.
`envelope_thermal_fraction_limit` can change this assessment bound after a
resolved-envelope or evolution comparison. It cannot exceed the selected time
energy tolerance or 0.5%, and is part of the restart identity. The run records
both the bound and the largest assessed fraction. Failed solves are excluded;
conservation checks remain unchanged. Above 1 MK, base nuclear power times the
layer mass must stay below 0.01%. These estimates are not universal bounds on
the spatial approximation.

`envelope_jacobian_radius "0.001"` optionally reuses nearby boundary derivatives
for Newton iteration. Boundary values are always reintegrated; derivatives refresh
after seven uses or a sufficient state change. The default is zero.

For H/He layers, `envelope_source` can supply a separate thermodynamic table.
The atmosphere's density inversion then uses that same source and its selected
metal approximation. Its tabulated temperature and pressure remain unchanged.
This keeps a cool atmosphere independent of the interior EOS temperature range;
the envelope still rejects unsupported source states. Analytic gray atmospheres
require a full EOS and retain their own thermodynamic calculation.
The native integration accumulates mass inward from the surface to reduce
roundoff in thin layers. Matched Hayashi and main-sequence checks agree within
1.196e-12; a 20 Myr cooling check changes luminosity by 0.003599%.
[Integration checks](results/envelope_density_integration_sept29_v1.json).
`envelope_map` can instead supply values and derivatives tabulated from native
integrations for one mass and envelope depth. Neither option extrapolates through
missing cells or outside composition coverage. These alternatives retain the
same evolution program and require their own physical and interpolation checks.

An optional Gibbs-potential source supplies density and all thermal derivatives
from the same function of temperature and pressure:

```text
envelope_eos "gibbs_hhe"
envelope_gibbs_hydrogen "hydrogen.dat"
envelope_gibbs_helium "helium.dat"
envelope_gibbs_hydrogen_warm "hydrogen_warm.dat"
envelope_gibbs_helium_warm "helium_warm.dat"
envelope_metals "neutral"
envelope_mass_fraction "0.00015"
```

Use this selection without `envelope_source`. Component potentials combine at
common gas pressure, convert atomic to baryon mass and include radiation once.
He3 and trace D use their nuclear number densities with the classical isotope
entropy correction. D must be at most 1e-4; molecular isotope shifts are omitted.
For Z up to 0.04, `neutral` or `ionized` adds ideal GS98 metals in the specified
charge state. The default `reject` permits only Z<=1e-12. These options do not
solve metal ionization, nonideal H/He mixing or element separation.

Five complete early-track envelope comparisons span X=0.1710–0.7000 and He3
up to 0.1013. Neutral versus ionized metals changes the base temperature by at
most 2.501%. A separate consistent-potential sensitivity to the published H/He
interaction correction changes it by up to 4.121%. Neither comparison is a
universal error bound; the interaction correction is not included in this
provider. A fresh evolutionary calculation and comparison with observed
low-mass stars are still required. [Mixture checks](results/gibbs_mixture_oct2_v1.json).

Each `EMBER_GIBBS_GAS_V1` file records the spline degree, knot vectors in ln T
and ln gas pressure, row-major coefficients of G/(R T), the original source
axes, and a Boolean source mask. Each numeric vector begins with its length;
the mask has one entry per source-grid node. Coefficients refer to atomic mass
and CGS units. Queries need a supported source cell and must remain within the
fitted potential. Warm joins must preserve the potential and its first two
derivatives. The four table hashes, isotope treatment and metal selection enter
restart identity; the tables are included in runtime packaging.

`envelope_integration_tolerance` selects checked implicit integration when
positive; its default is 1e-7 for `gibbs_hhe` and zero for the existing methods.
The implicit solver brackets temperature, falls back to bisection, and compares
a full pressure step with two half-steps. Negative thermal expansion uses the
signed buoyancy test; an outward-radiating stable layer remains radiative.
The integration method and tolerance enter restart identity. Their existence
does not validate a physical table or justify changing an evolved star's boundary
without assessing its earlier evolution.

With a positive tolerance, `envelope_integration_method "sdirk2"` selects a
second-order, L-stable pressure integration. Its two implicit stages have
diagonal coefficient `1 - 1/sqrt(2)` and retain the bracketed thermal solves.
Full and half pressure steps control error; the base-mass root uses the same
two-half-step solution as the accepted integration. The default remains
`backward_euler`, preserving existing configurations. This choice affects the
envelope calculation, not the stellar time-integration method.
The radius solve normally targets one tenth of the integration tolerance,
with a minimum relative accuracy of 1e-11.
If adaptive mesh changes prevent an exact radius root, the bracket must still
close within 1e-11 relative and the radius residual must be smaller than
the integration tolerance. Measured Hayashi cases reached machine-scale brackets
with residuals of 5.261e-10 and 2.104e-8; requiring a smaller residual stalled
without improving the physical accuracy.

Across twelve early and cooling envelopes, the second-order integration uses
about five times less CPU. Tightening its tolerance from 1e-7 to 1e-9 changes
base temperature, pressure and density by at most 4.988e-6 relative. Analytic
spherical hydrostatic and radiative-diffusion solutions test convergence
independently of the EOS tables. [Checks](results/envelope_second_order_oct2_v1.json).
Matched fresh Hayashi starts followed to 500 yr take 28.47 versus 132.2 CPU s
and 18.76 versus 74.16 wall s with eight zone threads. Both accept one interval
without rejection; all compared quantities agree within 0.1%. This short
control does not establish a speed ratio for an entire stellar lifetime.

Analytic tests check mass conventions, radiation, caloric energy, entropy
derivatives, missing cells and component joins. Seven complete envelopes from
656 to 4500 K agree with the independently assessed implementation within
4.435e-9. A bracketed secant solve reduces their CPU cost from 58.81 to 12.56 s.
An 8 Myr control with the existing EOS retains four accepted intervals and no
rejections; maximum global differences are below 3.841e-12. These checks validate
the implementation and preservation of the existing calculation, not the full
physical accuracy of the new tables. [Checks](results/gibbs_envelope_oct2_v1.json).

A complete fresh Hayashi-to-cold-white-dwarf calculation remains under validation.
Cool molecular chemistry, nonideal atmospheres and strong quantum effects are
not established by the existence of this boundary solver.

Envelope tests also check mass, energy, restart identity and the layer reservoir.
Source coverage and physical approximation checks remain necessary for each use.

A warm-H source repair follows the measured inconsistent-response region near
100,000 K, fitting the source density and entropy to one stable potential. It
retains the cold join and all previous coverage. In two saved envelopes the
base-temperature change is at most 0.003910%; the failed 4014 K envelope and
hotter 4250 K trial now complete. This does not validate the inconsistent source
derivatives or supply missing atmosphere cells. [Checks](results/gibbs_hydrogen_coverage_oct2.json).
