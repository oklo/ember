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

A complete fresh Hayashi-to-cold-white-dwarf calculation remains under validation.
Cool molecular chemistry, nonideal atmospheres and strong quantum effects are
not established by the existence of this boundary solver.

Validation includes86 test groups with the required local source data, a matched
1Gyr restart comparison with unchanged global quantities, and explicit mass,
energy and checkpoint tests for the envelope reservoir. Source-table coverage
and physical approximation checks are still required for each application.
