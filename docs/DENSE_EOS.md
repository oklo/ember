# Dense-matter thermodynamics

The material EOS derives pressure, energy, entropy and composition responses
from a free energy. An optional crystallization model keeps the same composition
in the liquid and solid. Element separation and a common EOS for extreme-cold
interiors and envelopes remain unfinished.

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

## Helium isotope mapping

`eos_helium_isotopes "number_density"` evaluates the classical source using
only its zero-helium-3 planes. With baryonic mass fractions and
`f = 1 + Y3/3`, the transformed state is

```text
rho' = f rho; X' = X/f; Z' = Z/f
Y3' = 0; Y4' = (Y4 + 4 Y3/3)/f
```

These transformations preserve every element's number density in FreeEOS's
classical isotope approximation. Specific thermodynamic potentials scale by
`f`; the physical isotope mixing and translational entropy are restored.
The normalization correction to the mixing-subtracted potential `F/T` is
`-R N_nuclei ln(f)`. Thermal and composition derivatives include the density
shift and this correction. Quantum and phase terms use the original physical
composition and density. This is not a new quantum-mixture or phase-separation
model. The default remains `"tabulated"`, and the option enters restart identity.

Eight direct source pairs agree within 5.009e-12. In 48 cold-envelope controls,
the supported thermal queries increase from 18 to 30 and the largest heat
capacity error falls from 95.58% to 0.008794%. Independent source differences
also check the chemical forces and enthalpies carried by diffusion. A cold
trace-metal enthalpy derivative still differs by 3.723%; its whole tested
Z = 1e-6 inventory corresponds to an energy difference of 7.131e-7 Cv T.
Thus this is not a uniform accuracy bound for all composition derivatives.
Missing thermal source cells remain masked.

Two matched 50 Myr intervals from a 744 K remnant change surface luminosity by
1.242e-8, peak helium-3 by 0.01465%, and local mass fractions by at most
3.232e-5, with the same convective mass. The largest internal luminosity change
is 0.1387% of surface luminosity near m/M = 0.9622. Species and energy checks
pass; the small redistribution is retained rather than called identical physics.

Removing unused isotope planes reduces this table from 1.900 GB to 475.0 MB.
All 4288 native queries are unchanged between full and compact mapped tables.
In that standalone probe peak memory falls from 3.816 GB to 0.9590 GB; this is
not a whole-star speed benchmark. See [checks](results/isotope_mapping_oct1_v1.json)
and [packing instructions](LIFETIME_DRIVER.md#compact-helium-isotope-tables).

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
Z ≤ 1e-8. He3 may be any share of the helium in every interval below: each ion
uses its own mass, and He3's quantum ratio is √(2/3) of the bounded hydrogen ratio.
[He3 range checks](results/he3_share_range_sept30_v1.json). Where trace helium is partly
recombined in the source (below about 110 kK), interpolating a helium-3-rich H envelope across
the hydrogen-fraction planes changes the tabulated Cv by up to 2.2% (0.07% for pure He4);
fully ionized layers reproduce the source at every share:

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

- T ≥ 35 kK, density 10–1000 g/cm³, X ≥ 0.95, Tp,H/T ≤ 1.7. Hydrogen stays
  fully ionized in direct pressure-ionized source checks of the cooling envelope
  (323 layer states, 35–200 kK, to Tp,H/T 1.65); the fitted correction and the leading
  Wigner–Kirkwood and Einstein-harmonic forms differ by ≤ 1.4% of classical ion Cv.
  Where trace helium recombines (below about 130 kK) the ionization change of the
  quantum Cv is ≤ 1.7e-3 of ion Cv at the 90th percentile and 2.9e-2 at most, at
  the recombination front. [Envelope controls](results/envelope_quantum_sept30_v1.json).
  These are model comparisons, not error bounds. The correction does not
  replace the ionization treatment in the underlying classical EOS.

- Nearly pure hydrogen, X ≥ 0.999: T ≥ 70 kK, density 10–1000 g/cm³,
  Z ≤ 1e-8 and Tp,H/T ≤ 2.5. The same quantum potential applies; only its
  assessed range is extended. In 92 direct equilibrium-ionization controls,
  hydrogen remains ionized. Allowing the trace helium to recombine changes
  quantum Cv by at most 0.008271% of classical ion Cv in the projected star.
  These controls have proton rs = 253.1–586.2, partly below the fit's lowest
  PIMC density parameter, rs = 500. This is a stated extrapolation, not a new
  simulation of a proton liquid. Leading Wigner–Kirkwood and harmonic model
  comparisons differ by up to 5.761% of classical ion Cv; their largest
  difference weighted over the newly admitted layers is 0.05611% of the
  projected star's heat capacity. These are sensitivities, not error bounds.
  The small bulk ionization effect does not validate trace-ion diffusion.
  [Controls and restart check](results/cold_hydrogen_oct1_v1.json).

- Dense H/He at any hydrogen fraction (the helium-rich mantle and the H/He
  transition): T ≥ 160 kK, density 300–4000 g/cm³, Z ≤ 1e-10 and
  Tp,H/T ≤ 1.85. Hydrogen and helium are fully ionized in direct
  source checks (quantum-Cv ionization sensitivity ≤ 1e-8 of classical ion Cv;
  ≤ 2e-9 in 69 transition-layer states at 162–222 kK and Tp,H/T ≤ 1.84, where the fitted,
  Wigner–Kirkwood and Einstein-harmonic corrections differ by ≤ 0.9% of classical ion Cv).
  The leading Wigner–Kirkwood and Hansen–Vieillefosse terms differ from the fitted
  correction by ≤ 7% and ≤ 2% of the quantum heat capacity; these are model
  comparisons, not bounds.

  The helium-rich source overlap (X ≤ 0.20) extends to 6000 g/cm³,
  Z ≤ 1e-8 and Tp,H/T ≤ 2.5, with the same 160 kK floor. This supplies
  the warm potential wherever the cold composition overlap has nonzero
  weight. In 27 supported source stencils, H and He remain fully ionized;
  the ionization sensitivity is below 6.317e-10 of classical ion Cv.
  Alternative quantum formulas differ by up to 3.505% of that Cv,
  which is model sensitivity rather than a physical error bound.
  [Source and overlap checks](results/dense_hhe_overlap_oct1_v1.json).

In the cold helium liquid, dilute hydrogen (mass fraction ≤ 1e-3) may reach
Tp,H/T ≤ 8 (otherwise 4). There the fitted correction is used outside its
R_S range (R_S,H ≈ 80). Against the leading Wigner–Kirkwood term and an
Einstein-harmonic form with the same classical limit, the quantum hydrogen
chemical potential differs by up to 0.52 and 0.23 kT at Tp,H/T = 7.6. The
native hydrogen force (chemical-potential difference plus the EOS and collision
enthalpy-temperature term) changes by less than 0.0032 kT per H ion across the
tested liquid layers. Relative changes are below 0.002 at the 90th percentile
and 0.0038 at the 99th; the largest ratios occur where the force is nearly zero.
All 352 dilute-H faces are covered down to half the temperature of the 1100 K
surface model, with density and composition fixed. At 48%, 90 central faces
cross the helium-liquid limit and are excluded. These comparisons bound the
observed model spread, not the physical error. The heat-capacity contribution
scales with the hydrogen fraction (at most about 0.001 of local Cv).
[Trace-hydrogen checks](results/trace_hydrogen_quantum_sept30_v1.json).
With the optional mixture phase model, material with X ≤ 1e-6 supports
Tp,H/T ≤ 12. This retains the same free energy and its derivatives. On a
784.5 K model and cooler temperature projections, replacing the proton quantum
term by an Einstein form changes local Cv by at most 1.613e-6 fractionally.
Removing or doubling the proton phase difference changes Cv by at most
8.737e-5 and the largest carried-heat flux by 3.307e-6 of the reference
surface luminosity. The hydrogen flux can change by 5.411%; its small mass
fraction alone does not establish an accurate diffusion force. These are
model sensitivities, not a calibration of quantum impurities or element
separation. Both phase branches receive the same Einstein correction;
the phase-difference variations are separate controls.
[Cold trace-H checks](results/cold_trace_hydrogen_oct1_v1.json).
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

The temperature join spans 500–800 kK; the hydrogen join spans X = 0.005–0.20.
An alignment term `a(rho, composition) + b(rho, composition)/T` matches the
two source potentials at both temperature anchors. Its density and composition
dependence carries model differences below the join, including differences
in chemical forces. The cold heat capacity is retained. All resulting thermal
and composition derivatives enter transport and energy accounting.

The composition overlap is broad enough to avoid artificial separation from
curvature of the joining weight. In 2066 supported H/He3 queries, the chemical
potential matrix at fixed temperature and pressure is positive definite; its
smallest eigenvalue is 5.309 R. This checks local mixture stability and includes
ideal mixing and density adjustment at fixed pressure. Source masks and physical
domain limits remain enforced. Pressure inversions and derivative tests pass.
Transport uses the joined chemical potentials and exchange enthalpies, so this
choice affects diffusion as well as structure. Short stellar controls pass their
conservation audits but encounter a separate convective switching problem;
a complete cooling-age comparison remains pending.
[Composition stability checks](results/cold_composition_join_oct1_v1.json).

The model requires density ≥ 1000 g/cm³, Z ≤ 0.16, helium coupling
Γ_He ≤ 130 (a practical bound below both assessed melting references: classical OCP
about 175, electron-screened helium 141–147 at 1e4–1e5 g/cm³), kT/E_F ≤ 0.05
and Tp,H/T ≤ 4. The number-weighted mean coupling is not a phase criterion here:
it is sensitive to strongly coupled trace metals (Γ_Fe ≈ 3700 at a 409 kK centre where Γ_He = 52).
The liquid fit continued to those metals overstates their Coulomb heat capacity
relative to a harmonic solid; integrated over the core this is +0.43% of ion Cv at
409 kK and +0.76% at 164 kK. Metal precipitation is not modelled: the calibrated
pure-Fe onset at the centre is about 405 kK (±20% against molecular dynamics for Fe
in C/O), and sinking all core Fe releases at most 3.4e43 erg in the current profile.
This is a conditional energy bound: the large He/Fe charge ratio and quantum
correction to separation are unbenchmarked. The estimated energy is at most
3.1% of the thermal release in the sampled cooling intervals; it is not a
measured cooling-age error.
[Core phase bounds](results/core_phase_sept30_v1.json). The helium quantum fit requires
500 ≤ R_S ≤ 1.2e5 and Tp,He/T ≤ 30. Trace-proton quantum terms use the same
linear-mixture approximation as the warmer EOS; comparison with leading-order
terms measures sensitivity, not a rigorous error bound. Tests cover source
responses, derivatives, join continuity and unchanged warm states. A replay from an unchanged warm model has cooled smoothly to 1550 K. At the
same age, interpolation of the comparison history gives a luminosity difference
of 0.1799%; this is an initial cooling comparison, not a final cooling-age test.

## Dense hydrogen–helium transition

`eos_dense_hhe_transition "liquid"` extends the same dense-material potential
to the metal-depleted H/He transition. It requires `eos_cold_helium` and leaves
that helium-core potential unchanged. The temperature overlap is 200–300 kK
and the density overlap is 300–600 g/cm³. Both weights and their derivatives
act on the free energy, including its composition derivatives. The two warm
source potentials supply the alignment described above; source values remain
required wherever their weight is nonzero.

The assessed range requires T ≥ 100 kK, density ≤ 6000 g/cm³, Z ≤ 1e-8 and
Tp,H/T ≤ 2.5, along with the existing degeneracy and helium-ion limits.
For X ≥ 0.999 and density ≤ 1000 g/cm³, the assessed lower temperature is
70 kK, using the hydrogen controls above. The other limits and warm anchors
remain in place. Together with the classical source repair, all 536 zones
of the 720.8 K profile remain supported at 0.8 times their temperatures.
That projection holds density and composition fixed; it is not a cooling age.
Direct source checks at 100–120 kK retain ionized H and He in the sampled
transition layers; the quantum-Cv sensitivity to ionization is below
5.643e-9 of classical ion Cv. Alternative quantum formulas differ by up to
3.648% of that Cv. This is a model spread, not an error bound. Colder
envelope source gaps remain unsupported.
Both warm anchors must also be supported; above 4000 g/cm³ the quantum
source extension requires X ≤ 0.20. Pressure inversions use this same
density range. On 536 stellar states and five cooler projections through
80% of their temperatures, all source and inversion checks pass; previously
supported values are unchanged. A 100 Myr continuation crosses the former
source boundary with five accepted steps and no rejections; the largest
relative first-law residual is 2.745e-8, with unchanged accuracy settings.
Direct source checks keep H and He fully ionized. Quantum prescriptions differ
by at most 2.772% of classical ion Cv over the admitted comparison points;
this is model sensitivity, not an error bound.

At the 790.3 K stellar model, pressure changes by less than 0.1%, local Cv
by up to 6.269%, and the adiabatic gradient by up to 10.92%. The different
Coulomb prescriptions and the potential overlap contribute to these changes.
A tabulated Cv anomaly at a colder projected state is checked against direct
FreeEOS values rather than treated as a physical reference. All 536 current
states and three cooler projections through 88% of their temperature pass
EOS and pressure-inversion checks; colder states retain the physical limits.
A matched 50 Myr comparison changes luminosity by 0.02604%, with conservation
checks passing. This is not a complete cooling-age test.
[Source, derivative and stellar comparisons](results/dense_hhe_transition_oct1_v1.json).

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

The coupled transport table covers electron degeneracy −2 ≤ eta ≤ 2048
and screening coordinate 0.01–150. At eta > 256 the source generator uses
1536 points for polynomial normalization and electron–ion integration,
increasing the count above eta = 1024 to resolve the Fermi surface over the
larger energy interval. Explicit quadrature settings remain available. Independent normalization,
source reproduction and native transport checks accompany the
[table extension](results/electron_pair_eta2048_sept30_v1.json).
The low-order transport response is tested; the highest raw operator moments
are not converged. A matched 100 Myr stellar interval changes luminosity by
5.103e-13 relative.
This extends the same screened Born collision model; it does not add
strong-correlation or relativistic corrections.


## Optional freezing at fixed composition

`eos_cold_helium "mixture_softmin"` adds liquid and bcc-solid ion free energies
from Potekhin & Chabrier (2013) and Baiko & Chugunov (2022). A smooth minimum,
with width 0.005 kT per ion, supplies the latent heat through the same potential
as pressure and composition forces. The composition dependence of that width
is differentiated. The default remains `liquid_mixture`.

The forced-solid fit is excluded below helium coupling Gamma = 90. The solid
branch must be negligible at that boundary, evaluated at the local density
and composition. Its weight may then increase continuously above the boundary.
[Boundary and cooling-step checks](results/mixture_phase_cut_oct1_v1.json).

`eos_metal_liquid_continuation_gamma "200"` selects a thermodynamically
consistent continuation of the classical liquid fit. Beyond the stated coupling,
its internal energy follows the value and slope at the join; integrating it
gives the free energy. Values 175–300 are available for sensitivity tests; zero
retains the original fit. Both settings enter restart identity. This avoids
using the original liquid fit's unphysical extrapolation below the bcc ground
energy, but does not establish a multicomponent phase diagram.

On the 835.7 K surface-temperature profile, the continued liquid changes total
heat capacity by −0.4283%; the phase term is still negligible. The choices
175, 200 and 300 place the central transition at approximately 188, 183 and
171 kK. These are declared scenarios, not measured uncertainty bounds. Actual
freezing may separate elements. Independent metal velocities, crystal diffusion
and precipitation are not selected by this EOS option. Do not extend a
liquid-ion transport prescription through appreciable freezing without assessing
that approximation. Hydrogen and quantum-domain guards remain active.

The selected Ioffe electron-conduction table already treats liquid and crystal
scattering. A random-impurity comparison on the current profile changes the
required temperature drop across its inner 0.9775% by mass from 140.8 to
499.9 K. This is a fixed-profile sensitivity, not a cooling-age error estimate.

A matched 100 Myr comparison begins with 25 Myr steps and passes every time
and conservation check. With the same liquid continuation, the phase term
changes central temperature by 1.667e-7 relative and luminosity by
1.152e-11: freezing has not yet appreciably begun. With the phase option off,
the original physical checkpoint is reproduced. The thermal and composition
chain-rule implementation agrees with direct five-variable differentiation
within 9.115e-12 over 2680 profile queries. Initial larger-step failures are
retained in the [numerical checks](results/mixture_phase_sept30_v1.json).

## Ion motion during freezing

`solid_ion_mobility_fraction` optionally reduces the full ion mobility matrix
by `1 - (1 - f) w`, where `w` is the same solid weight used by the EOS and
`f` is the remaining fraction of liquid mobility in a fully solid region.
It requires `mixture_softmin`; omitting it leaves transport unchanged.
Values zero and 0.01 provide an immobile limit and a finite-mobility comparison.
They are scenarios, not fitted diffusion coefficients for a helium crystal.

The factor and its thermal and composition derivatives enter species fluxes
and their carried heat together. Scaling the whole matrix preserves its mass
and current constraints. Electron conductivity remains in the selected
conduction treatment. The option is recorded in restart identity.

[Hughto et al. (2011)](https://arxiv.org/abs/1104.4822) find diffusion through
both nearly perfect and defective Coulomb crystals. Consequently, freezing
does not establish permanent chemical trapping. Species-dependent crystal
diffusion and precipitation require additional physics; this option cannot
decide whether lead forms a separate core.

On a matched 100 Myr interval beginning at 814.4 K, both mobility scenarios
change global quantities by less than 5e-9 relative to unchanged transport.
Every interval is step-doubled; the largest first-law residual is 1.725e-8.
The final structure residual and Newton correction are each limited to 2e-6,
200 times tighter than the temporal structure target. Profile checks also
cover nearly solid states. [Numerical checks](results/solid_mobility_sept30_v1.json).
