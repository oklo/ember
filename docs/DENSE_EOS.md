# Dense-matter thermodynamics

The material EOS derives pressure, energy, entropy and composition responses
from a free energy. Crystallization, phase separation and a common EOS for
extreme-cold interiors and envelopes remain unfinished.

Tables can contain separate valid density intervals. Pressure inversion stays
in the interval containing its density guess and cannot cross a masked cell;
without a guess it uses the dilute interval.

Missing composition-plane vertices can be repaired with converged source
evaluations and complete derivative stencils, without changing valid table
values. A repair of 305 warm, dense helium-rich vertices passes 2144 native
queries. Newly supported states agree with direct FreeEOS responses within
1.152e-5; a matched 62.89 Myr stellar interval differs by at most 3.247e-7
in global quantities. These are interpolation and evolution checks, separate
from the physical uncertainty of the EOS.
[Source-support checks](results/cooling_eos_support_sept30_v1.json).

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
enthalpies all derive from that potential. The full published expression is
used at larger plasma-temperature ratios; its series avoids cancellation
when the ratio is small.

The assessed range requires number-weighted ionic coupling ≤ 100. Tp,H/T is
limited to 1, except in a helium-rich core with T ≥ 500 kK, density 1e3–1e5
g/cm³, X ≤ 0.05 and Z ≤ 0.16, where the assessed limit is 4.
The lower-density helium extension passes 136 direct equilibrium-ionization
checks: hydrogen and helium remain fully ionized. Partial metal ionization
produces an electron deficit up to 0.156%; that approximation remains.
[Ionization checks](results/helium_quantum_domain_sept30_v1.json).
Nearly pure-helium layers also permit Tp,H/T up to 4 at T ≥ 400 kK,
density 1e3–1e5 g/cm³, X ≤ 0.01 and Z ≤ 0.001. In 285 direct ionization
checks, H and He remain fully ionized and the metal electron deficit is below
`9.792e-6`. Thermal and composition derivative checks pass. This extends the
same liquid potential; the metal-rich core limits and mean-coupling limit
remain fixed. It does not establish phase stability of the central mixture.
[Cool-helium checks](results/cool_helium_quantum_domain_sept30_v1.json).
Common-electron-density linear mixing follows the prescription discussed by
[Baiko (2022)](https://academic.oup.com/mnras/article/517/3/3962/6712718).
At the assessed 1786 K model, the quantum term reduces the core heat-capacity
integral by 2.229%. Higher terms change it by 0.05104% relative to the leading
quantum term, rising to 0.3069% in a cooler fixed-density projection. These
comparisons measure sensitivity, not a proved mixture error bound or cooling-age
error. The liquid treatment does not determine whether concentrated metals
separate or freeze before the helium. The mean-coupling guard is not a mixture
phase diagram. [Source, derivative and core checks](results/quantum_liquid_full_sept30_v1.json).
Below 300 kK, Tp,H/T must be ≤ 0.1 except in the intervals below. Outside the helium
core, Tp,H/T ≤ 1 unless an interval states otherwise. The hydrogen intervals require
Z ≤ 1e-8; He3 is limited to half the helium unless stated otherwise:

- T ≥ 200 kK, density 50–150 g/cm³ and X ≥ 0.97. Equilibrium-ionization checks
  give an electron deficit below 2.785e-5.
  [Dense-hydrogen checks](results/quantum_dense_hydrogen_sept28_v1.json).
- T ≥ 50 kK, density 10–200 g/cm³, X ≥ 0.98 and Tp,H/T ≤ 0.7. Equilibrium-ion
  comparisons give quantum-Cv differences up to 0.1347% of classical ion Cv
  in the warmer controls, 0.03264% through Tp,H/T = 0.5, and 0.4033%
  through Tp,H/T = 0.7. Pressure differences stay below 0.002253%. At 2488 K the
  additional envelope layers contain 0.03875% of the resolved mass; their
  ionization uncertainty contributes less than 1e-8 of integrated stellar Cv.
  The full quantum correction reduces that Cv by 0.5597%.
  [Warm controls](results/quantum_cool_hydrogen_sept29_v1.json) ·
  [Cooling-envelope checks](results/quantum_cooling_envelope_sept30_v1.json). A matched 20 Myr
  comparison changes luminosity by -0.815%; this includes the initial EOS
  readjustment and is not a cooling-age measurement.
  Above 200 kK this H-rich domain extends to 500 g/cm³, with the same quantum
  ratio and composition limits. Direct equilibrium-ionization controls change
  the quantum Cv by less than 3.090e-10 of classical ion Cv, at the numerical
  differencing floor. The colder layers retain the 200 g/cm³ limit.
  [Dense-envelope controls](results/cold_transport_domain_sept30_v1.json).

- T ≥ 35 kK, density 10–500 g/cm³ (up to 900 above 200 kK), X ≥ 0.99,
  Tp,H/T ≤ 1.1 and He3 at most 60% of the helium. Hydrogen stays fully ionized
  in direct source checks down to 35 kK. Direct ionization controls
  support the hydrogen-rich cooling layers. Trace-helium pressure ionization
  makes the fully ionized quantum correction uncertain by up to 1.15% of
  classical ion Cv at X = 0.99; the corresponding stellar-integrated sensitivity
  in the tested cooler profile is 5.806e-7 of total Cv. This is a model comparison,
  not a rigorous uncertainty bound. Denser sampling also finds peaks of 0.78%
  in the X = 0.98 interval above. Neither interval supplies an ionization model
  for the underlying classical EOS.
  [Cold-envelope controls](results/quantum_cold_envelope_sept30_v1.json).

- Dense H/He at any hydrogen fraction (the helium-rich mantle and the H/He
  transition): T ≥ 200 kK, density 300–4000 g/cm³, Z ≤ 1e-10, He3 at most 75% of
  the helium and Tp,H/T ≤ 1.4. Hydrogen and helium are fully ionized in direct
  source checks (quantum-Cv ionization sensitivity ≤ 1e-8 of classical ion Cv).
  The leading Wigner–Kirkwood and Hansen–Vieillefosse terms differ from the fitted
  correction by ≤ 7% and ≤ 2% of the quantum heat capacity; these are model
  comparisons, not bounds.

In the cold helium liquid, trace hydrogen (mass fraction ≤ 1e-6) may reach
Tp,H/T ≤ 6.5 (otherwise 4). There the fitted correction is used outside its
R_S range. Its heat capacity is small but nonzero (about −1e-7 of the local Cv).
The R_S sensitivity of the hydrogen chemical potential is about 0.01 kT
(rising to 0.025 kT in colder projections).
[Source and stellar checks](results/quantum_dense_mixture_sept30_v1.json).

Outside those intervals, the smaller quantum-ratio limit bounds the ideal-ion
heat correction to 0.05556%. The potential is never tapered at a boundary;
unsupported states are rejected. This remains a perturbative liquid-ion
correction, without crystallization or a complete partially ionized quantum EOS.

Potential, pressure, energy and Cv agree with the independent `LIQUBC` source
within 6.104e-11; finite differences check thermal and composition derivatives.
A 150 Myr cooling comparison changes luminosity by 0.05775% without changing
the convective extent. Its stellar use with the integrated envelope remains
local. [Component and evolution checks](results/quantum_liquid_sept28_v1.json).

## Cold liquid mixture

`eos_cold_helium "liquid_mixture"` optionally replaces the dense, helium-rich
material free energy below 800 kK. It combines degenerate electrons through
the Sommerfeld T² term, classical Coulomb ions, electron polarization and
exchange-correlation from [Potekhin & Chabrier (2013), appendices A–C](https://arxiv.org/abs/1212.3405),
and the existing Baiko–Chugunov quantum correction. Radiation and ideal mixing
are included once by the surrounding EOS. This is a liquid approximation;
it does not select crystallization or latent heat.

The temperature join spans 500–800 kK; the hydrogen join spans X = 0.01–0.05.
An alignment term `a(rho, composition) + b(rho, composition)/T` matches the
two source potentials at both temperature anchors. Its density and composition
dependence carries model differences below the join, including differences
in chemical forces. The cold heat capacity is retained. All resulting thermal
and composition derivatives enter transport and energy accounting.

The model requires density ≥ 1000 g/cm³, Z ≤ 0.16, number-weighted coupling
≤ 100, kT/E_F ≤ 0.05 and Tp,H/T ≤ 4. The helium quantum fit requires
500 ≤ R_S ≤ 1.2e5 and Tp,He/T ≤ 30. Trace-proton quantum terms use the same
linear-mixture approximation as the warmer EOS; comparison with leading-order
terms measures sensitivity, not a rigorous error bound. Tests cover source
responses, derivatives, join continuity and unchanged warm states. A replay from an unchanged warm model has cooled smoothly to 1550 K. At the
same age, interpolation of the comparison history gives a luminosity difference
of 0.1799%; this is an initial cooling comparison, not a final cooling-age test.

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

## Electron collision sources

The coupled transport table covers electron degeneracy −2 ≤ eta ≤ 1024
and screening coordinate 0.02–150. At eta > 256 the source generator uses
1536 points for polynomial normalization and electron–ion integration;
explicit quadrature settings remain available. Independent normalization,
source reproduction and native transport checks accompany the
[table extension](results/electron_pair_eta1024_sept30_v1.json).
This extends the same screened Born collision model; it does not add
strong-correlation or relativistic corrections.
