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
It covers 4600–5400 K and log g = 5.9–6.2. The forward extension uses 12
independently checked columns on the same fine opacity grid. The 4800 K,
log g = 6.2 entry interpolates the solved 4700 K and 5000 K columns; a separate
attempt to calculate it directly reached its CPU limit and was not used.
`sources.json` records the source columns and each interpolated entry.

The extension preserves all 50 original numerical rows and all 5999 saved
atmosphere queries exactly. Another 100 queries test the added domain and its
derivatives. Three actual stellar intervals and their final physical checkpoint
are identical before and after the extension; the production continuation also
preserves its starting physical state. Values are continuous at the former
4800 K and log g = 6.1 boundaries, while the new cells supply different
one-sided slopes. This is the existing piecewise logarithmic interpolation.
[Validation](../../../docs/results/highg_atmosphere_extension_sept27_v1.json).

Trace atmospheric helium remains approximated. The high-gravity metallicity
range is still 0 to 1e-20; separate small-metal source controls have not changed
that domain. Future models outside the declared composition, temperature or
gravity bounds require additional coverage.

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
