# Experimental CMS19 pressure/entropy EOS

`Cms19Eos` supplies molecular, partially ionized and dense H/He thermodynamics
for **static structure only**. Its actual interpolation derivatives pass
independent numerical checks, but its pressure and entropy do not form a
thermodynamically exact pair. The source audit below is a prerequisite for
future evolution work, not a solved issue.

## Data and mixture

The source and reproduction commands are in [data/eos/README.md](../data/eos/README.md).
The [CMS19 paper](https://arxiv.org/abs/1902.01852) describes additive-volume
mixing at fixed temperature and material pressure, with a constant ideal
H/He mixing entropy that neglects free-electron mixing. Our implementation
uses that convention and adds radiation pressure and entropy once.

For fixed hydrogen mass fraction X:

```
v(T,Pg) = X/rho_H(T,Pg) + (1-X)/rho_He(T,Pg)
S_material = X*S_H + (1-X)*S_He + S_mix
n_H = X/1.0078; n_He = (1-X)/4.0026; n = n_H+n_He
S_mix = -R_gas * sum(n_i*ln(n_i/n))
rho = 1/v
P = Pg + a_rad*T^4/3
S = S_material + 4*a_rad*T^3*v/3
```

By default metals are rejected. The driver explicitly chooses
`Metals::helium_proxy`, so X=.7, Y=.28, Z=.02 becomes X=.7 and effective
Y=.3 **for the EOS only**. The opacity retains the actual X=.7, Z=.02
GS98 mixture. He3 must be zero; no species abundances or free electron
fraction are inferred from the EOS. This precludes composition evolution
with this implementation even apart from the missing caloric state.

## Interpolation and derivatives

We interpolate log rho and log S on the original log T/log Pg axes using
four-node monotone cubic Hermite interpolation along pressure, then along
temperature. The query is always in the middle interval, with one extra
node on each side. Harmonic slope limiters carry both first and second
derivatives through their active branch using `src/jet2.hpp`. Derivatives
are thus derivatives of the implemented values; separately interpolating
the source's response columns would not provide that guarantee. Second
derivatives can jump at cell boundaries and limiter changes.

Let t=ln T, p=ln Pg, r=ln rho, and q=ln P_total. Subscripts denote actual
interpolant derivatives in t,p. After mixing and adding radiation:

```
cp       = S_t - S_p*q_t/q_p
cv       = S_t - S_p*r_t/r_p
delta    = -r_t + r_p*q_t/q_p
grad_ad  = -S_p/(q_p*cp)
chiRho   = q_p/r_p
chiT     = q_t - q_p*r_t/r_p
Gamma1   = chiRho*cp/cv
```

Here cp and cv are T times the corresponding entropy derivatives.
`grad_ad` follows the constant-entropy direction. First derivatives of cp,
delta and grad_ad use the interpolation Hessians and are transformed to
ln T,ln rho for the analytic stellar Jacobian. Fixed-composition mixing
entropy has no t,p derivatives. Finite-difference checks include actual
dense ionization states visited by the converged star.

Pressure inversion brackets the supported material pressure interval and
uses safeguarded Newton/bisection, accepting only a log-density error below
2e-13. Direct (T,P_total) inversion subtracts radiation and rejects unresolved
or nonpositive material pressure. There is no extrapolation or analytic
fallback. States must have positive finite cp, cv, compressibility, delta
and adiabatic gradient.

## Physical masks and atmosphere

The source README warns that its rectangles include unphysical states.
The runtime restricts queries to log T=3.2..7.3 and requires all stencil
nodes to have -8 <= log rho <= 4, finite usable entropy, and to pass the
classical-fluid/ion-quantum masks. It uses the paper's Wigner–Kirkwood bound
eta²/24 < .7 and melting estimates. Neutral helium uses the Simon curve;
above rho=10 g/cm³ or T=1e6 K the implementation uses the ionized OCP
estimate. These choices define conservative computational support, not a
new phase diagram or a proof of accuracy throughout the retained domain.
The H2 melting peak is below the selected temperature floor.

For each temperature interval, the reader retains one contiguous pressure
interval whose whole stencil passes. Mixtures require the intersection of
the participating pure components. Density bounds are derived from that
pressure interval. The grey atmosphere now brackets its starting state
inside both EOS and opacity support. A default tau_top=1e-6 reaches below
the CMS19 density floor for the stellar atmosphere and correctly fails.
The explicit CMS19 driver default is tau_top=.001, matched at tau=2/3.
Doubling it to .002 is tested on an identical converged mass mesh.

## Reproducible consistency defects

For exact equilibrium thermodynamics,

```
grad_ad = P*delta/(rho*T*cp)
D = P*delta/(rho*T*cp*grad_ad) - 1 = 0.
```

This code does **not** force the equality. The equilibrium JSON reports
the maximum |D|, its location, and the mass-weighted RMS. For the 4096-point
0.1 Msun reference, max |D| is about .24 in the outer dense ionization
layers, while mass RMS is .00435. A small mass RMS does not demonstrate
that envelope structure or luminosity is insensitive to the larger local
defect. These errors are separate from the small stellar residual and
virial error.

`scripts/audit_cms19_energy.py` reads original archive columns without any
ember interpolation. At the pure-H node log T=4.45, log P[GPa]=1.35:
log rho=-1.034, log S[MJ/kg/K]=-.9503,
dlnrho/dlnT at P=-1.107, dlnS/dlnP at T=-.1154. These source columns
give D=-.2650. Thus the discrepancy cannot be attributed entirely to
ember's interpolation. Its interpolant contributes additional differences
and must still be audited against any corrected source.

The same script checks original He internal energies at rho=1 g/cm³:

| T (K) | U (erg/g) |
|---:|---:|
| 891300 | 8.886e+13 |
| 1000000 | 8.227e+13 |

The energy secant is **-6.059e+7 erg/g/K**, while the upper source node's
entropy-derived cv is **+9.725e+7 erg/g/K**. These columns cannot jointly
support a reliable energy equation. Exploratory reconstructions from
U-TS or pressure integration also developed negative heat capacities at
joins and were rejected; none entered the library or imported data.

Energy, its density derivative, mu and free-electron abundance are explicitly
unavailable (NaN). Positive-dt zone and central energy calculations reject
the EOS before evaluating energy differences. Enabling evolution requires
an independently validated caloric EOS with consistent pressure, entropy
and energy through dissociation/ionization and source joins. A corrected
source or a thermodynamic-potential implementation must also supply the
composition support needed for He3 production and eventual H exhaustion.
