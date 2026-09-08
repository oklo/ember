# CMS19 H/He tables

`cms19_h_tp.dat` and `cms19_he_tp.dat` contain the original density and
entropy columns from `TABLE_H_TP_v1` and `TABLE_HE_TP_v1` in the authors'
[DirEOS2019.tar.gz distribution](https://perso.ens-lyon.fr/gilles.chabrier/DirEOS/DirEOS2019.tar.gz),
retrieved 2026-09-07. Cite Chabrier, Mazevet & Soubiran (2019), ApJ 872, 51,
[A new equation of state for dense hydrogen-helium mixtures](https://arxiv.org/abs/1902.01852).
The archive root is `DirTABLES-EOS2019/`. Its README identifies the pure
component files as February 2019 v1; the June 2021 archive update added a
mixture table. These external scientific data are not relicensed under
ember's MIT code licence.

The authors' June 2021 interacting-mixture README recommends the 2019
additive-volume tables for low-mass stars. We use pure H and He, never the
2021 effective-hydrogen table as pure H. The optional metals-as-helium
approximation must be explicitly selected in `Cms19Eos`.

## Reproduction

```
curl -fLO https://perso.ens-lyon.fr/gilles.chabrier/DirEOS/DirEOS2019.tar.gz
python3 scripts/import_cms19.py DirEOS2019.tar.gz data/eos
python3 scripts/audit_cms19_energy.py DirEOS2019.tar.gz
```

The importer and audit verify SHA-256
`736de2a0b26c02b897fbd906504bb6aa144f08a6d0f3cff36da62970054d149b`.
Reimport was checked byte for byte against both versioned files. Python's
standard library is sufficient; tables are read directly by C++ at runtime.

## Format and support

Each file has 121 temperatures by 441 pressures, or 53,361 pairs of original
values. The format is:

```
121 441 CMS19 H: <description>     # HE: for helium
<121 log10 T[K] values>
<441 log10 P[GPa] values>
<log10 rho[g/cm3]> <log10 S[MJ/kg/K]>  # pressure varies fastest
...
```

Axes run from log T=2 to 8 and log P[GPa]=-9 to 13, both in steps of .05.
The importer retains the original printed density/entropy strings, including
unphysical rectangular corners. It neither interpolates nor repairs them.
The runtime supports a much smaller fluid subset: log T=3.2..7.3, with
each entire 4×4 interpolation stencil inside density and phase masks.
Consequently an imported rectangle does **not** imply rectangular physical
coverage. See [CMS19 implementation and audit](../../docs/CMS19.md).

## Static calculations only

Energy columns are deliberately absent. At rho=1 g/cm³, the original He
table's internal energy decreases between log T=5.95 and 6.0 despite a
positive entropy-derived heat capacity. The source H table also violates
the pressure/entropy Maxwell identity by 26.5% at one audited dense
ionization point. These are reproducible diagnostics of this particular
source version, not a claim that every CMS19 state has such errors.

`Cms19Eos::has_internal_energy()` is false. Energy and its density
derivative are NaN; time-dependent structure and central boundary calls
reject this EOS. Density, entropy and their derivatives support an
**experimental static calculation**, with consistency errors exposed in
the equilibrium JSON. Do not manufacture internal energy or enforce one
identity by overwriting an independently derived response.
