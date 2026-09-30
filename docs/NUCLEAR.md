# Nuclear reactions

The lifetime calculation evolves initial deuterium, the reduced pp chains and
carbon–nitrogen burning, with Solar Fusion III pp rates and finite-degeneracy
screening. Checkpoints record the selected prescriptions. The network omits
pep, hep, ppIII, oxygen branches and pycnonuclear burning.

## Network and energy accounting

Helium-3 evolves explicitly. Primordial deuterium is a separate inventory;
deuterium produced by pp reactions is included in their reduced rates.
The reduced ppII branch assumes rapid completion with available protons and
uses a 0.861 MeV neutrino energy.

The CN network evolves carbon-12, carbon-13 and nitrogen-14 through

- C12(p,γ)N13(β⁺)C13;
- C13(p,γ)N14;
- N14(p,γ)O15(β⁺)N15(p,α)C12.

Short beta decays and nitrogen-15 capture are eliminated from the explicit
network. The catalyst abundances are separate from the fixed GS98 mixture
used by the material tables. Proton capture changes the physical hydrogen,
helium and metal inventories; the nuclear rest-energy correction accompanies
that change. A complete cycle converts four protons to helium-4, but the
intermediate catalyst abundances need not be in equilibrium.

Reaction coefficients are moles of reactions per gram per second. Identical
reactants receive the factor 1/2 once. Atomic rest masses determine Q values,
with escaping neutrino energy subtracted. Heating and abundance derivatives
use the same reaction assembly. Species and nuclear-energy checks do not
renormalize abundances or clip the energy release.

## Rates

`sfiii-svh` uses the pp normalization, slope and curvature from equations 8–9
of [Solar Fusion III](https://arxiv.org/abs/2405.06470v3), the helium-3 pair
polynomial from section V.C, and the helium-3/helium-4 energy dependence from
equation 14. The latter is restricted to its stated range through 1.6 MeV.
The operational temperature range is 0.1–20 MK: burning is set to zero below
the lower cutoff and unsupported hotter states are rejected.

The retained `sfii-svh` option uses
[Solar Fusion II](https://arxiv.org/abs/1004.2318), Table I and equations 25–26,
with `S(E) = S0 + S1 E + S2 E²/2`:

| Reaction | S0 (MeV barn) | S1 (barn) | S2 (barn/MeV) |
|---|---:|---:|---:|
| pp | 4.01e-25 | 11.2 × S0 / MeV | 0 (omitted) |
| He3 + He3 | 5.21 | −4.9 | 22 |
| He3 + He4 | 5.60e-4 | −3.60e-4 | 1.51e-4 |

Integration uses two 48-node Gauss–Legendre panels in log energy around the
Gamow peak. It returns `N_A ⟨σv⟩` in cm³ mol⁻¹ s⁻¹ and the thermal response
`−3/2 + ⟨E/kT⟩`. Laboratory screening is excluded from the bare S-factors.

The pp-only command defaults to `sfii-svh`; the library's `PPChains()` default
retains historical fits for static comparisons. These defaults do not select
the lifetime calculation's physics. `sfii-debye`, `sfii-legacy-screening` and
`legacy` provide comparison prescriptions.

## Screening

The Salpeter–Van Horn prescription combines finite-degeneracy Debye screening
with the ion-sphere limit:

```text
Theta_e = (kT/ne) (dne/dmu_e)_T
D^-2 = 4 pi e²/(kT) [sum(ni Zi²) + ne Theta_e]
Gamma_e = e²/(ae kT),  ae = [3/(4 pi ne)]^(1/3)
h_DH = Z1 Z2 e²/(D kT)
h_S = .9 Gamma_e [(Z1+Z2)^(5/3) - Z1^(5/3) - Z2^(5/3)]
h_SVH = h_DH h_S / sqrt(h_DH²+h_S²)
f_screen = exp(h_SVH)
```

See [Potekhin & Chabrier (2013), equations 10–13](https://arxiv.org/abs/1310.3162).
Relativistic ideal-electron Fermi integrals supply the susceptibility and its
responses. The model assumes fully ionized matter and is approximate at
intermediate coupling. It is not derived from the interior FreeEOS potential.

Classical-ion evaluation requires `zeta = 3 Gamma_12/tau ≤ 0.2`. The optional
`quantum_screening_zeta_max` adds
`h_CDW(Gamma,zeta) − h_CDW(Gamma,0)` to the SVH exponent, using
[Chugunov & DeWitt (2009), equations 23–25 and A4](https://arxiv.org/abs/0905.3844).
This combination retains the SVH electron response and omits the
energy-dependent S-factor shift; it is not the full CDW prescription.

The option defaults to zero (off), is recorded in checkpoints and accepts
limits through 1.6. Evaluation stops beyond the selected limit or
`Gamma_12 = 200`. This does not supply pycnonuclear burning.
[Quantum-screening checks](results/quantum_screening_sept28_v1.json).

For a depleted cooling core, optional `quantum_burning_fuel_limit` instead
omits an out-of-domain reaction when the local sum of H, D and He3 mass
fractions is below that limit. It defaults to zero (strict refusal) and
accepts at most `1e-6`. Supported reactions are unchanged. Each omitted
reaction has zero rate, heat and derivatives; its fuel remains available
for transport and later burning. Fuel-rich cells still refuse evaluation.
This approximation does not extend the screening formula or supply a
pycnonuclear rate.

The driver records affected cells, mass and their current fuel-energy
ceiling, including subsequent burning of daughter protons. This inventory
does not bound future fuel delivery by diffusion. Selecting the option
requires an independent estimate of missing power and a finite cooling
interval; a small fuel fraction alone does not establish negligible heat.
The option and its limit are part of the restart identity.
[Profile checks and a bounded cooling test](results/cold_trace_burning_sept30_v1.json)
record the current evidence and the rate-estimate limitations.

## Verification

Independent adaptive energy integration checks the rates and thermal responses;
separate Fermi integration checks screening. Tests cover abundance derivatives,
zero helium-3, mass conventions, baryon and nuclear-energy conservation,
timestep refinement and rejection outside supported domains.

[SFIII rate checks](results/nuclear_sfiii_bare_rate_comparison_v1.json) ·
[Stellar comparison](results/nuclear_sfiii_evolution_control_v1.json) ·
[Reference generator](../scripts/generate_sfiii_reference.py) ·
[Run configuration](LIFETIME_DRIVER.md)
