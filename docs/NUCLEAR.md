# Nuclear rates and screening for bounded pp evolution

`ember-evolve` now selects `sfii-svh`: direct Solar Fusion II rate integrals
and finite-degeneracy Salpeter–Van Horn screening. The library's default
`PPChains()` retains the old fits and capped classical screening for existing
static benchmarks. Its name now identifies that choice; the earlier attribution
of those fits to Adelberger et al. was insufficient. Every evolution output
records the selected prescription. This update supplies a reproducible nuclear
model for continued initial evolution, not a complete hydrogen-burning network
or precision screening tied to the FreeEOS potential.

## Solar Fusion III option

`sfiii-svh` selects the pp-chain rates recommended by
[Acharya et al. (2025), Solar Fusion III](https://arxiv.org/abs/2405.06470v3),
with the same Salpeter–Van Horn screening and reduced reaction network.
Equations 8–9 supply the pp normalization, slope and curvature. Section V.C
supplies the helium-3 pair polynomial, and equation 14 supplies the full
helium-3/helium-4 energy dependence. The latter is restricted to its stated
range through 1.6 MeV; its negligible high-energy tail is not extrapolated.
The operational temperature range remains 0.1–20 MK, with zero burning below
the existing low-temperature cutoff and rejection above the upper limit.

The 18 rate and temperature-response comparisons use independently adaptive
integration in energy, including a tail check. Both prescriptions also undergo
the same abundance, thermal-response and nuclear-energy conservation tests.
The exact-response cache distinguishes the rate prescriptions. Changing this
selection changes the nuclear physics and requires a consistent new trajectory;
the currently running star retains `sfii-svh`.
[Bare-rate comparison](results/nuclear_sfiii_bare_rate_comparison_v1.json),
[reference calculation](../scripts/generate_sfiii_reference.py).
All 34 test suites pass. Two 512-point controls complete 10 billion years;
the SFII control retains the preceding executable's entire output exactly.
With SFIII, the initial luminosity increases by 0.7454% and radius by 0.2164%.
This is a short stellar comparison, not a measurement of the full lifetime.
[Validation](results/nuclear_sfiii_validation_v1.json),
[stellar controls](results/nuclear_sfiii_evolution_control_v1.json).

This option does not add pep, hep, ppIII or CNO reactions, or update screening.

## Carbon and nitrogen capture times

A new diagnostic evaluates proton capture by carbon-12, carbon-13 and nitrogen-14
at three saved central states using the low-energy expansions in Solar Fusion II
Table XII and the same approximate screening as Ember. At 3.685 trillion years,
the instantaneous capture times are 23.80 billion years for carbon-12,
5.020 billion years for carbon-13 and 27.37 trillion years for nitrogen-14.
The central temperature is 9.290 million K. Independent energy-space and
log-energy quadratures agree within the recorded numerical criterion.
[Calculated times and inputs](results/cno_capture_times_v1.json),
[diagnostic script](../scripts/audit_cno_capture_times.py).

At the plotted 3.848-trillion-year endpoint, the corresponding central times
are 445.6 million years, 93.69 million years and 289.6 billion years. The
central temperature is 11.75 MK. The diagnostic now supports 1–14 MK while
retaining its original quadrature and tail criteria; all sampled integration
energies remain below the 130 keV range stated for the nitrogen polynomial in
Solar Fusion II, equation 50.
[Endpoint comparison](results/cno_capture_times_3848gyr_v1.json),
[rate-range checks](results/cno_rate_range_14mk_v1.json).

These times do not establish a negligible CNO contribution. They describe
fixed conditions, not the star's changing temperature, mixing or depletion.
Carbon conversion to nitrogen needs a changing-composition calculation before
its heat and composition effects can be bounded. The GS98 elemental pattern
used by the current EOS and opacity must also be reconsidered if carbon and
nitrogen change. No CNO reactions have yet been installed, and this diagnostic
does not claim a total CNO luminosity or a lifetime correction.

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
SFII/legacy rate ratios are .9880, .9334 and 1.063 for pp, 33 and 34.
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
