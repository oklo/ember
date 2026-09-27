# Very hydrogen-rich gas atmospheres

These sources extend the common atmosphere selection at optical depth 100.
They do not select the optical-depth-25.12 white-dwarf boundary used in the
separate pulse experiments. `sources.json` records each retained input and hash.

`low.dat` contains 27 solved gas atmospheres; `middle.dat` contains 72.
`trace.dat` retains explicit hydrogen and metal dependence, neglecting the
atmospheric isotope correction below a declared helium-3 limit. The nearly
pure-hydrogen sources neglect trace atmospheric helium in the matching
pressure and temperature. Density always uses the actual stellar composition.

`high.dat` includes seven solved anchors and **inferred** temperature/gravity
rows, as declared in its header. It is an explicitly bounded continuation,
not a fully computed cool, high-gravity grid. Its coverage ends at 4600 K
and log g = 6.1. New source calculations must test and replace those inferred
rows before the resulting late track can establish atmosphere-independent
helium-3 ignition. Helium mass fractions above the declared limits are rejected.

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
