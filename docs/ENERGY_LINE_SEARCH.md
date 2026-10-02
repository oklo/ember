# Energy balance in the Newton line search

The structure solver checks both local residuals and total luminosity balance.
The line search now uses both, scaled by their existing tolerances. Previously,
a correction could improve total energy balance but be rejected because a
smaller local residual had stopped decreasing. No equation or tolerance changes.

A 1 Gyr interval starting at 3.247 Tyr used identical physics and two zone threads
per calculation. CPU time fell from 88.41 to 31.90 s; wall time from 74.43 to
26.21 s. Accepted intervals fell from six to two, rejected trials from three
to zero. The largest compared global difference was 8.294e-7 in luminosity;
the largest pointwise change in ln radius, density or temperature was 8.307e-8.
These measurements apply to this interval, not the whole track.

Both calculations retained the 5e-8 total-energy limit. An ideal-gas stellar
benchmark with a finite-precision boundary reproduces the old failure and
converges with the new line search. Its energy balance is checked independently
by integrating heating and the change in internal energy plus pressure work.
The relaxation, PMS/CN, face-energy and evolution-controller tests pass.

Numerical results: [energy_line_search_oct2.json](results/energy_line_search_oct2.json).
