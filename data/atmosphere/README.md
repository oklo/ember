# AMES-COND non-grey tau=100 boundary

The source `tau100_T.data` and `tau100_Pgas.data` files are from
[MESA](https://github.com/MESAHub/mesa/tree/fd396fd73d3f936da8063ffdf9d92361882eb557/atm/atm_data),
commit `fd396fd73d3f936da8063ffdf9d92361882eb557`, retrieved 2026-09-07.
They supply temperatures in K and material pressures in dyn/cm² at
Rosseland optical depth 100. Their header identifies solar GN93, [M/H]=0.
The selected cool-dwarf region comes from AMES-COND-2000 structures
(Allard et al. 2001), as documented by
[MESA's atmosphere documentation](https://docs.mesastar.org/en/latest/atm/table.html).
These scientific data are not relicensed under ember's MIT code license.

The full MESA rectangle includes grey integrations, extrapolation and
transition smoothing. The importer **only selects Teff=1800..3300 K in
100-K steps and log g=3.5..6.0 in .5-dex steps** (96 nodes). Inspection of
[`create_tau100.f90`](https://github.com/MESAHub/mesa/blob/fd396fd73d3f936da8063ffdf9d92361882eb557/atm/preprocessor/src/create_tau100.f90)
establishes that this region retains COND values: grey integration is at
log g≤2, C&K overwrite begins at 3500 K, and transition smoothing alters
3400/3500 K. None of those cells is imported. The converter adds PHOENIX's
electron pressure to its heavy-particle pressure, so the distributed
`Pgas` is material pressure; ember adds radiation once.

The source's detailed mixture and helium abundance are not matched to
ember's X=.7,Z=.02 interior. The output is explicitly a **solar-mixture
proxy**, never labeled an exact-composition atmosphere. Version 2 records
both that approximation and the allowed interior composition. The reader
rejects it unless constructed with `Mixture::allow_documented_proxy`;
the CLI selects this with `--atmosphere cond-solar-proxy`.

```
python3 scripts/import_cond_atmosphere.py data/atmosphere/sources \
  data/atmosphere/cond_gn93_tau100_solar_proxy.dat
```

Pinned source SHA-256:

- Pgas: `77606d68e2f0c02602624e04e31792c5150c7a15940b53cc8a364ba3682f5900`
- T: `e9551fac9e2d8d036ab94c55ed87cb66811d2a0c54e82d27c846c1aabdea4c3a`
- MESA converter inspected: `ab5f7d105bf10983d674b998b641ecaa6319f0a32bc7484a42d5b83e490ae57f`

No spectrum file was used to infer pressure. These are extracted boundary
states from model structures, not a complete atmosphere transfer solver.
The interface retains the thin-atmosphere approximation: the tau=100
radius is identified with the photospheric radius and the omitted mass is
neglected. Resolving that geometric extension and matching modern mixture
and atmosphere physics remain further improvements.

## Bounded initial evolution

The underlying table and fixed-composition reader remain strict, including
abundance basis. `ember-evolve` explicitly selects
`FrozenCompositionAtmosphere`: the same GN93 T/P relation with density
recomputed from the actual evolving EOS. It rejects |Xsurface-.7|>.005,
He3>.005, changed metals or unsupported Teff/gravity. These operational
limits do not establish an atmosphere-error bound. A composition-dependent
non-grey source grid is still needed for substantial depletion; see
[EVOLUTION.md](../../docs/EVOLUTION.md).
