# Very hydrogen-rich gas atmospheres

These sources extend the common atmosphere selection at optical depth 100.
They do not select the optical-depth-25.12 white-dwarf boundary used in the
separate pulse experiments. `sources.json` records each retained input and hash.

`low.dat` contains 27 solved gas atmospheres; `middle.dat` contains 72.
`trace.dat` retains explicit hydrogen and metal dependence, neglecting the
atmospheric isotope correction below a declared helium-3 limit. The nearly
pure-hydrogen sources neglect trace atmospheric helium in the matching
pressure and temperature. Density always uses the actual stellar composition.

`high.dat` contains solved gas atmospheres and explicit interpolated knots.
The complete 4600 K gravity row and the 4625 K, log g = 6.1 column now have
independently checked source solutions. Older warm anchors retain their
recorded provenance. The grid covers 4600–4800 K and log g = 5.9–6.1;
`sources.json` identifies the source calculations and remaining interpolation.
Trace atmospheric helium is still approximated, and helium mass fractions
above the declared limits are rejected.

The latest source replacement preserves all 5932 checked early boundary
queries exactly and supports all 67 saved late structures. On the late
structures, matching temperature changes by at most 0.2367% and pressure by
1.101%. A three-step stellar replay and its final physical checkpoint are
identical before and after the update. The actual continuation also preserves
its initial physical state. These checks establish compatibility over the
sampled structures; future models still require coverage checks. See
[the integration checks](../../../docs/results/highg_atmosphere_integration_sept27_v2.json).

The manifest selects composition and gravity intervals. Each interval queries
covered source endpoints, preserving the original full-composition boundary
outside its declared composition overlap. There is no failure-triggered source
selection or interpolation across missing source cells.

The retained helium approximation was checked at helium mass fraction 0.005
at 4650 K and 4800 K, log g = 5.9: matching temperature changes were below
0.1% and pressure changes below 2.1%. See
[the source comparison](../../../docs/results/hydrogen_envelope_integration_sept27_v1.json).
Those local controls do not establish accuracy throughout the cooler,
higher-gravity region or after a structural instability.
