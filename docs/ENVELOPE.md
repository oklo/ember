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
from10kK to2MK. The absolute reservoir energy change plus pressure work must stay
below0.1% of luminosity per interval, including when the base is cooler than1MK.
Above1MK, base nuclear power times the layer mass must stay below0.01%. These diagnostics supplement conservation checks and
are not a universal bound on the spatial approximation.

`envelope_jacobian_radius "0.001"` optionally reuses nearby boundary derivatives
for Newton iteration. Boundary values are always reintegrated; derivatives refresh
after seven uses or a sufficient state change. The default is zero.

For H/He layers, `envelope_source` can supply a separate thermodynamic table.
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
