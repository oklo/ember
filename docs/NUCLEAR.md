# Nuclear rates and screening for bounded pp evolution

`ember-evolve` now selects `sfii-svh`: direct Solar Fusion II rate integrals
and finite-degeneracy Salpeter–Van Horn screening. The library's default
`PPChains()` retains the old fits and capped classical screening for existing
static benchmarks. Its name now identifies that choice; the earlier attribution
of those fits to Adelberger et al. was insufficient. Every evolution output
records the selected prescription. This update supplies a reproducible nuclear
model for continued initial evolution, not a complete hydrogen-burning network
or precision screening tied to the FreeEOS potential.

## Bare rates and units

The new rates integrate the Maxwell–Gamow kernel using the S-factor expansion
`S(E)=S0+S1*E+S2*E²/2`. Inputs are from Table I and equations 25–26 of
[Adelberger et al. (2011), Solar Fusion II](https://arxiv.org/abs/1004.2318):

| Reaction | S0 (MeV barn) | S1 (barn) | S2 (barn/MeV) |
|---|---:|---:|---:|
| pp | 4.01e-25 | 11.2 × S0 / MeV | 0 (omitted) |
| He3 + He3 | 5.21 | -4.9 | 22 |
| He3 + He4 | 5.60e-4 | -3.60e-4 | 1.51e-4 |

SFII does not supply a recommended pp curvature; setting it to zero is an
explicit truncation, not a claim that it vanishes. Later pp calculations and
rate compilations remain a future update. Kinematic masses are physical atomic
masses minus electron masses; their tiny electronic binding correction is
omitted. Screening is applied separately, with no laboratory screening folded
into these S-factors.

Ember integrates in log energy with two 48-node Gauss–Legendre panels around
the Gamow peak. It returns `N_A <sigma v>` in cm³ mol⁻¹ s⁻¹ and the analytic
thermal response `-3/2 + <E/kT>`. Nuclear abundance coefficients are **moles of
reactions per gram per second**; identical reactants get the factor 1/2 once.
`NA` enters the neutrino energy conversion once. Atomic rest masses determine
all Q values, regardless of the abundance basis. The new common reaction
assembly forms heat and all abundance/thermal derivatives from those same Q
values, avoiding subtraction of nearly equal species mass flows for heating.
Neither abundance renormalization nor energy clipping is used.

The SFII expansion is restricted operationally to T<=2e7 K. The existing
T<1e5 K zero-burning cutoff remains explicit. These choices do not implement
cold dense pycnonuclear burning. The reduced ppII path retains its old .861-MeV
neutrino approximation and assumes rapid completion with available protons;
pep, hep, ppIII and CNO are absent. The profile audit reports the ppII/pp
reaction ratio to measure the importance of this particular reduced branch.
It does not measure errors from the missing channels.

## Screening and electron response

The new default uses

```
Theta_e = (kT/ne) (dne/dmu_e)_T
D^-2 = 4 pi e²/(kT) [sum(ni Zi²) + ne Theta_e]
Gamma_e = e²/(ae kT),  ae = [3/(4 pi ne)]^(1/3)
h_DH = Z1 Z2 e²/(D kT)
h_S = .9 Gamma_e [(Z1+Z2)^(5/3) - Z1^(5/3) - Z2^(5/3)]
h_SVH = h_DH h_S / sqrt(h_DH²+h_S²)
f_screen = exp(h_SVH)
```

These are the finite-degeneracy Debye susceptibility and the Salpeter–Van Horn
interpolation reviewed in equations 10–13 of
[Potekhin & Chabrier (2013)](https://arxiv.org/abs/1310.3162).
The interpolation connects weak and ion-sphere limits smoothly; it is an
approximation in intermediate coupling. Their low-mass-object comparison finds
screening-exponent differences of order tens of percent for this simple
prescription. That is **not** an error bound for ember's pp track.

The relativistic ideal-electron Fermi integrals already used in ember now
provide susceptibility moments too. A checked density inversion determines
eta, `Theta_e`, and its thermal/density derivatives. Composition derivatives
include the electron-density change as well as each ion's charge contribution.
The classical nondegenerate limit is not substituted for a partially degenerate gas.
Fully ionized abundances are assumed. The current forward calculation uses
the shared GS98 elemental inventory for metal ion and electron counts; the
retained legacy inventory uses the nuclear module's carried metal charges.
Ionization effects in negligible-burning outer layers are omitted.
This screening model is not derived from the interior's FreeEOS Coulomb
potential and does not establish thermodynamic consistency between the two
physical approximations.

For the new screening choices, `zeta=3*Gamma_12/tau` must be <=.2. This is an
operational classical-ion restriction, not a fitted accuracy bound;
it rejects quantum-dominated states rather than capping their enhancement.
The 20-Gyr track's actual central and maximum active-layer zeta are recorded
in the audit. The final central zeta is .0503 and the maximum above the
burning cutoff is .1004. Quantum corrections and a full
polarization/mixing free-energy treatment remain unimplemented. See the quantum-regime discussion in
[Potekhin & Chabrier (2012)](https://arxiv.org/abs/1201.2133).

Three screening choices are available. `sfii-svh` is the evolution default;
`sfii-debye` retains the finite-degeneracy weak formula as a sensitivity control;
`sfii-legacy-screening` uses SFII rates with the old classical exp(2) cap.
`legacy` uses both old rates and old screening. A control is not an alternative
validated stellar prediction. At the 20-Gyr model's central state the bare
SFII/legacy rate ratios are .9880, .9334 and 1.0628 for pp, 33 and 34.
The pp screening exponents are .5160 (legacy), .4461 (finite-degeneracy weak)
and .2801 (SVH), with `Theta_e=.4311`. Screening materially changes the star;
the old 10-Gyr numerical reference must not be relabeled as the new physics.

## Independent checks and reproduction

The versioned `tests/data/nuclear_reference.dat` contains 18 rate/thermal-response
queries and 10 independent Fermi/screening queries. Regeneration uses SciPy's
adaptive energy-space quadrature and bracketed electron inversion, independently
of the C++ log-energy panels and Newton solve. The largest rate/response
relative difference is 1.23e-11; the screening comparison is 6.12e-13. These
are numerical checks of the stated formulas, not independent physical data.

```
# Only optional reference regeneration needs SciPy.
python3 scripts/generate_nuclear_reference.py /tmp/nuclear-reference.dat
cmake --build build
ctest --test-dir build --output-on-failure
c++ -std=c++23 -O2 -Iinclude scripts/nuclear_probe.cpp build/src/libember.a -o /tmp/ember-nuclear-probe
python3 scripts/audit_nuclear.py /tmp/ember-nuclear-probe out/evolution-sfii-svh-4096-20gyr.json docs/results/nuclear_m010_20gyr.json
```

Additional tests check derivatives against independent differences for all
species, zero He3, both abundance bases, and mixed pp branches; density and
mass-convention transformations; baryon and nuclear rest-energy conservation;
the dilute classical limit; smooth growth beyond the old cap; and unsupported
input rejection. The coupled test now uses SFII/SVH, checks timestep refinement,
and verifies rollback when burning exceeds actual composition coverage.

A 1024-point `legacy early` reproduction reaches 10 Gyr in the same 23
macrosteps as the old checkpoint; final R, L and He3 agree within 2e-11
relative after the EOS performance changes.

## Extended main-sequence calculation

The atmosphere/opacity/conduction update retains this nuclear prescription.
The [one-trillion-year audit](results/nuclear_m010_1tyr.json) repeats the
screening, rate and branch diagnostics at the final evolved mixture; see
[EXTENDED_EVOLUTION.md](EXTENDED_EVOLUTION.md). He3 builds to about .102
by baryonic mass near 720 Gyr and subsequently declines. The explicit
non-equilibrium He3 abundance remains necessary throughout this calculation.

The omitted pep channel was reconsidered. SFII equation 46 gives a pep/pp
fit with an explicit 10–16 MK validity range. This star's approximately
5 MK core is outside that range and has partially degenerate electrons.
Blind extrapolation would introduce an unsupported physical input. A rate
and electron-capture screening treatment applicable to this regime is still
required; neither the small reduced ppII branch nor the plasmon-only loss
audit establishes an error bound for missing nuclear channels.
