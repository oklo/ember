# Five-times-reference-metallicity controls

![Nuclear support and luminosity](Ember_boundary_5Z_experimental.png)

A 0.06-solar-mass model with Z = 0.08241 and Y = 0.27 approaches sustained
hydrogen burning in these experimental clear-atmosphere calculations.
Matching the atmosphere at Rosseland optical depth 3 gives nuclear support
near 100% through 891.4 Gyr, at 1204 K. Luminosity and temperature increase
slowly after their minima; no large surge appears. The original and
regularized H/He EOS fits give nearly overlapping curves.

Matching at optical depth 1 produces nonconvective regions and some surface
metal depletion. This arm reaches its fixed-composition atmosphere limit
at 24.80 Gyr, with 98.61% nuclear support. The depth-3 arms reach their
hydrogen-composition allowance near 891.4 Gyr. Repeated step reductions at
that boundary ultimately also trigger the envelope thermal-reservoir check;
neither arm completes 1 Tyr. Open circles mark these limits, not cessation
of hydrogen burning.

These are post-deuterium starts. Initial deuterium heat and the formation
clock are omitted. Atmospheric chemistry and its imposed adiabat use fixed
source mixtures; finite mixing-length convection supplies the deeper
envelope. The EOS refit removes interpolation oscillations but changes a
real source bump. These qualifications preclude adding a new transition-mass
point to the main boundary figure. Cloudy controls and composition-dependent
atmosphere matching remain to be assessed.

Every evolutionary interval retains its time-error and conservation checks.
A 10–15 Gyr comparison of 100-Myr and 1-Gyr step ceilings changes luminosity
by 0.02046% and no convective boundary. The 1-Gyr ceiling reduces elapsed
time from 42.22 to 12.39 seconds in that comparison.

`plot.py` reproduces the figure from `tracks.csv`; `provenance.json` records
source histories, restart checks, endpoints and assumptions. It requires
NumPy and Matplotlib. No uncomputed continuation is drawn.
