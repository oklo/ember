# Dense-matter thermodynamics

The material EOS derives pressure, energy, entropy and composition responses
from a free energy. Crystallization, phase separation and a common EOS for
extreme-cold interiors and envelopes remain unfinished.

Tables can contain separate valid density intervals. Pressure inversion stays
in the interval containing its density guess and cannot cross a masked cell;
without a guess it uses the dilute interval.

## Degenerate electrons

`ElectronGas` and `IdealEos` integrate relativistic Fermi–Dirac occupations.
For degeneracy parameter eta > 100, paired quadrature about the Fermi surface
avoids cancellation in heat capacity and thermal-pressure derivatives.
Occupation entropy is evaluated directly. The cold limits are

```text
cv = pi^2 n_e k_B^2 T / (rho p_F v_F)
dP/dlnT = rho T cv (2+x_F^2) / [3(1+x_F^2)]
x_F = p_F/(m_e c)
```

See [Baiko & Yakovlev (2019)](https://www.ioffe.ru/astro/Stars/Paper/baiko_yakovlev19mn.pdf).
Tests cover thermal derivatives, first-law and Maxwell relations, and the
warm/cold join. At 100 K and densities 1e3–1e6 g/cm³, heat capacity agrees with
the degenerate limit to about 1e-10. At 10,000 K and 1000 g/cm³ the difference
is 8.2e-7, including finite-temperature corrections.
[Electron checks](results/cold_electron_response_v1_audit.json).
The tabulated metal-bearing FreeEOS material is a separate implementation.

## Ion interactions

A comparison at a 2511 K helium-WD model uses the same H/He/GS98 element numbers
in FreeEOS and [Skye](https://arxiv.org/abs/2104.00691), with quantum ions compared
separately. FreeEOS gives 5.440% less integrated Cv T dm over the inner 90% of
mass. The Skye quantum term contributes −0.4648%. Controlled component tests
attribute most of the classical difference to ion interactions. This thermal
scale is not a measured cooling-age error. An improved prescription must enter
the free energy, including its composition derivatives; rescaling Cv is insufficient.
[Comparison](results/dense_core_eos_comparison_sept28_v1.json).

`ion_classical_liquid_jets` supplies the liquid ion-ion free energy used by
Potekhin & Chabrier (2000), with linear mixing and thermal/composition derivatives.
The nuclear screening calculation shares the same fit. Individual ions beyond
Gamma = 200 use a constant-entropy continuation; this does not choose a stable
phase. Electron polarization and nonlinear mixture terms remain to be added.
The component is not selected in the stellar EOS and must not be added on top
of FreeEOS's existing Coulomb contribution.
[Independent reference and derivative checks](results/classical_ion_sept28_v1.json).

## Weakly quantum liquid ions

The optional `LiquidIonQuantumPotential` adds the liquid-ion free-energy term of
[Baiko & Chugunov (2022)](https://doi.org/10.1093/mnras/stab3613).
Pressure, energy, entropy, heat capacity, chemical forces and transported
enthalpies all derive from that potential. A series through the fourteenth
power avoids cancellation at small plasma-temperature ratio.

The assessed range requires Tp,H/T ≤ 1 and number-weighted ionic coupling ≤ 100.
Below 300 kK, Tp,H/T must be ≤ 0.1 except in two H-rich intervals, both requiring
Z ≤ 1e-8 and He3 at most half the helium:

- T ≥ 200 kK, density 50–150 g/cm³ and X ≥ 0.97. Equilibrium-ionization checks
  give an electron deficit below 2.785e-5.
  [Dense-hydrogen checks](results/quantum_dense_hydrogen_sept28_v1.json).
- T ≥ 50 kK, density 10–200 g/cm³, X ≥ 0.98 and Tp,H/T ≤ 0.5. Equilibrium-ion
  comparisons give quantum-Cv differences up to 0.1347% of classical ion Cv
  in the warmer controls and 0.03264% in the added controls. At 2488 K the
  additional envelope layers contain 0.03875% of the resolved mass; their
  ionization uncertainty contributes less than 1e-8 of integrated stellar Cv.
  The full quantum correction reduces that Cv by 0.5597%.
  [Warm controls](results/quantum_cool_hydrogen_sept29_v1.json) ·
  [Cooling-envelope checks](results/quantum_cooling_envelope_sept30_v1.json). A matched 20 Myr
  comparison changes luminosity by -0.815%; this includes the initial EOS
  readjustment and is not a cooling-age measurement.

Outside those intervals, the smaller quantum-ratio limit bounds the ideal-ion
heat correction to 0.05556%. The potential is never tapered at a boundary;
unsupported states are rejected. This remains a perturbative liquid-ion
correction, without crystallization or a complete partially ionized quantum EOS.

Potential, pressure, energy and Cv agree with the independent `LIQUBC` source
within 6.104e-11; finite differences check thermal and composition derivatives.
A 150 Myr cooling comparison changes luminosity by 0.05775% without changing
the convective extent. Its stellar use with the integrated envelope remains
local. [Component and evolution checks](results/quantum_liquid_sept28_v1.json).

## Independent source derivative checks

`build_eip_probe.py --phase-probe` evaluates the pinned
[Ioffe EOS EIP source](https://www.ioffe.ru/astro/EIP/eipintr.html) at specified
density and phase, excluding ideal electrons and radiation.
`audit_eip_derivatives.py` checks free-energy, pressure and energy derivatives.
Forced metastable phases test formula consistency, not physical applicability.

Two corrections are needed in that source: `EOSFI22` uses the quantum temperature
pressure derivative `PDTQL` where `PDRQL` is required, and `LIQUBC` omits part of
the density dependence of its coefficients. For each mode,

```text
y = C_i(r_s) T_p/T;  D_i = dln C_i/dln r_s;  a_i = 1/2 − D_i/3
u_i = y df_i/dy;  D_i' = dD_i/dln r_s
p_i = a_i u_i
p_i + dp_i/dlnT = a_i c_i
p_i + dp_i/dlnrho = (a_i+a_i^2+D_i'/9) u_i − a_i^2 c_i
D_1' = −D_1(1−D_1);  D_2' = 0;  D_3' = D_3(3D_1−2D_3−1)
```

`--repair-quantum-pressure-derivatives` builds an isolated corrected source,
retaining the original, patch and comparison. Ember differentiates the free
energy itself. [Derivative audit](results/eip_phase_derivatives_consistent_v2_audit.json).
