# Cloudy models at three times reference metallicity

The two panels show nuclear support and luminosity from the same saved models,
with temperature coloring the lower panel. All tracks include metal ionization
in the interior EOS and begin on Hayashi contraction. Initial Z = 0.04945,
Y = 0.27; clouds use f_sed = 2 and clear between 900 and 1200 K.

The 0.060 and 0.061 solar-mass tracks use the same EOS family. The 0.0605
track uses its wider metallicity extension because settling raises the core
metal fraction beyond the former table range. It is shown as a separate
comparison, not used to tighten the mass interval.

The dashed continuation beyond 1.087 Tyr widens the permitted atmospheric
H/He mismatch without changing boundary values. Its estimated validity remains
provisional: its atmosphere composition response remains under assessment.
Local PICASO calculations now include the corrected electron-dependent
continuum opacity; this alone does not validate the wider composition range.
The plot stops at the last model inside the 1150 K assessment range; the next
accepted state is slightly outside and is omitted.

The earlier 0.058-solar-mass surge figures use a different, neutral-metal EOS.
They should not be read as the current prediction of the transition mass.

Run `python plot.py` with NumPy and Matplotlib. `tracks.csv` contains the plotted
states; `provenance.json` identifies retained histories and their checksums.
