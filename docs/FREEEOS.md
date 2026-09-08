# A consistent fixed-composition Helmholtz EOS

`HelmholtzTableEos` derives pressure, entropy, internal energy, heat capacities
and transport responses from **one C2 material free-energy interpolant**.
It replaces the inconsistent pressure/entropy pair in the CMS19 experiment.
Radiation is added analytically once. The current data source is
[Alan Irwin's FreeEOS 3.0.0](https://sourceforge.net/projects/freeeos/files/freeeos/3.0.0%20Source/),
using its recommended EOS1 options `(3,1,-2)`. This includes molecules,
partial and pressure ionization, electron degeneracy, Coulomb and exchange
terms. It is an alternative physical model, not repaired CMS19 data.

The numerical H/He mixture retains X=.7 and effective Y=.3 for an isolated
comparison with the previous CMS19 run. Metals are still explicitly treated
as helium, and He3 is unsupported. All queries must match the recorded
fixed interior composition `solar_scaled(.7,.02)`. Opt-in to this mixture
approximation is required at construction. Composition evolution, arbitrary
metallicity, and a full nuclear abundance treatment remain future work.

## Potential and derivatives

The material potential is phi = F/T, in natural logarithmic coordinates
x=ln T and y=ln[rho/(T/1e6)^1.5]. Each node supplies the nine derivatives
phi_(i,j), i,j=0,1,2. Tensor-product quintic Hermite interpolation matches
these jets on a rectangle. phi and its first/second derivatives are
continuous across cells; third derivatives can jump.

P, E, S and the source response coefficients give the potential and its
Hessian at source nodes. Fourth-order centered differences of the Hessian
complete the higher mixed derivatives. The importer removes FreeEOS's
radiation using its derived CODATA radiation constant before forming jets.
Ember adds its own analytic radiation term phi_rad=-a*T^3/(3*rho).

After transforming derivatives to t=ln T at fixed rho and r=ln rho:

```
P = rho*T*phi_r
E = -T*phi_t
S = -phi - phi_t
cv = -phi_t - phi_tt
chiT = 1 + phi_tr/phi_r
chiRho = 1 + phi_rr/phi_r
delta = chiT/chiRho
cp = cv + phi_r*chiT*delta
grad_ad = phi_r*delta/cp
dE/dlnrho = -T*phi_tr
```

Third potential derivatives supply analytic derivatives of cp, delta and
grad_ad. The first-law and Maxwell identities follow from the potential,
and are independently tested by differences of the returned E, S and P.
No independently interpolated entropy or heat-capacity field is overwritten
to make an identity pass.

The current source grid has 305 temperatures × 593 density coordinates,
spaced by .0125 dex: log10 T=3.2..7.0, log10 Q=-5..2.4. Omitting two source
rows on every edge leaves 301 × 589 = 177,289 potential nodes. Their
computational range is log10 T=3.225..6.975, log10 Q=-4.975..2.375.
The importer checks convergence, requested density/temperature, and local
first-law identities before accepting source states. 106 potential nodes
are excluded because their complete derivative stencils encounter a source
state with a nonpositive transport response. Runtime uses only the
contiguous branch connected to the dilute edge, rejects extrapolation,
and also rejects any unstable interpolated state. These are computational
support conditions, not a new phase diagram or an accuracy guarantee for
every state in the rectangle.

## What consistency does and does not establish

All 129 initial direct FreeEOS samples along the old star converged;
their local thermodynamic identities agreed to approximately 1e-12.
However, a finer audit exposed **small discontinuities in the source fits**.
In particular, `mod_molecular_hydrogen.f90` switches its Irwin H2/H2+
partition-function polynomials at 9000 K. Along log10 Q=.882, integrating
the source derivative over log10 T=3.95..3.9625 differs from the endpoint
change in F/T by about 1640 erg/g/K. Local differences away from the join
still match the source derivatives. Thus exact local source identities
alone would have missed this issue.

The C2 potential **regularizes these fit joins over the interpolation
cells**. This is a declared approximation, not an exact reproduction of
every source response. Against 912 off-grid direct evaluations (512 old
stellar-profile points plus 400 seeded random points), the .0125-dex
potential differs by at most .0022% in P, .049% in E, .012% in S, 1.42%
in cp, 1.66% in cv, 1.08% in grad_ad, and 1.70% in delta. These are sampled bounds;
they are not uniform error bounds over the table or physical uncertainties.
Refinement does not monotonically reduce errors at the source-fit joins.
The coarser .025-dex and finer .0125-dex 1024-point stellar models differ
by about .0030% in radius and .021% in luminosity. Do not call this
asymptotic EOS-grid convergence. Smoother underlying molecular/other fits
would be preferable for precision calorimetry and oscillation work.

Internal energy is available for the **implemented fixed-composition
potential**, and its time-dependent zone derivatives and central boundary are tested.
This does not enable burning evolution: composition changes are rejected,
and mixing, composition derivatives and time-step control remain absent.
CMS19 retains its independent static-only guard and its original audit.

## Reproduction

Normal builds use the versioned C++ table and need no Fortran or network.
The source numerical evaluations are retained under `data/eos/sources/`.
Rebuild the potential with standard-library Python:

```
python3 scripts/import_freeeos_potential.py \
  data/eos/sources/freeeos300_hhe_x070_raw.json.gz \
  data/eos/freeeos300_hhe_x070_potential.dat
```

To regenerate the source evaluations with gfortran, CMake and Ninja:

```
curl -fL 'https://downloads.sourceforge.net/project/freeeos/freeeos/3.0.0%20Source/free_eos-3.0.0.tar.gz' -o /tmp/ember-freeeos300.tar.gz
python3 scripts/build_freeeos_probe.py /tmp/ember-freeeos300.tar.gz /tmp/ember-freeeos-source-build
python3 scripts/generate_freeeos_grid.py /tmp/ember-freeeos-source-build/probe /tmp/freeeos-raw.json.gz
python3 scripts/import_freeeos_potential.py /tmp/freeeos-raw.json.gz /tmp/freeeos-potential.dat
```

The archive SHA-256 is
`4ab1c15a51385a3eab3b08c6f3f240739c0105d92ec828d635ac95720edefb09`.
The builder makes one compiler-portability adjustment when extended and
quadruple precision have the same kind: it removes duplicate generic
diagnostic overloads. No source EOS formula is changed. A fresh probe
process per isotherm avoids a demonstrated cached-state failure at the
dense-to-dilute row reset. Every retained calculation must return info=0;
there is no alternative-EOS fallback. Across compilers, compare physical
values within numerical tolerances rather than demanding identical last bits.

The independent source-value and molecular-join audit can also be repeated:

```
c++ -std=c++23 -O2 -Iinclude scripts/helmholtz_probe.cpp build/src/libember.a -o /tmp/ember-helmholtz-probe
build/apps/ember-equilibrium 4096 .1 .15 --eos cms19 --hot-opacity tops --seed-index 1.5 > out/equilibrium-m010-reference.json
python3 scripts/audit_freeeos.py /tmp/ember-freeeos-source-build/probe /tmp/ember-helmholtz-probe out/equilibrium-m010-reference.json
```

The script queries each source state in a fresh process, checks convergence
and the local identities, compares the 912 off-grid states, and integrates
the derivative across the 9000-K molecular fit join. The recorded output is
[results/freeeos_audit.json](results/freeeos_audit.json). Its source identity
residual is 5.98e-13 while the finite-interval potential mismatch remains
1639.73 erg/g/K, demonstrating why both checks are needed.

FreeEOS source is GPL-2.0-or-later and remains in the external build
directory. The versioned files are numerical outputs and ember's own
interpolation/import code; FreeEOS is not linked into ember. Cite Irwin
and the FreeEOS release when using these calculations.
