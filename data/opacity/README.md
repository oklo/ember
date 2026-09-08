# Opacity tables

## `ferguson_gs98_z020.dat`

Rosseland mean opacities from **Ferguson et al. (2005), ApJ 623, 585**,
"Low-Temperature Opacities", for the Grevesse & Sauval (1998) metal mixture at
Z = 0.020, assembled from the distribution at Wichita State
(`webs.wichita.edu/physics/opacity`, files `f05.gs98.tar.gz`).

Eight hydrogen fractions (X = 0.0, 0.1, 0.2, 0.35, 0.5, 0.7, 0.8, 0.9) × 85
temperatures (log T = 2.70 to 4.50, i.e. 501 to 31,623 K) × 19 densities
(log R = -8.0 to +1.0 in steps of 0.5), where

    log R = log rho - 3 (log T - 6).

Values are log10 of the Rosseland mean in cm^2/g. These include molecular
bands and grain condensation, which is what makes them usable below 1700 K
where the older Alexander, Johnson & Rypma (1983) grid had to be
extrapolated — and where that extrapolation returned kappa ~ 1e-9, removing
the Hayashi limit from cool giants entirely.

Cite Ferguson et al. (2005) in any work using these.

## `opal_gs98_z020.dat`

OPAL Type 1 Rosseland means (Iglesias & Rogers 1996, ApJ 464, 943), GS98
mixture, Z=0.020. Original distribution: `GS98hz`, dated 20011028, obtained
2026-09-07 from Arnold Boothroyd's
[CITA opacity distribution](https://www.cita.utoronto.ca/~boothroy/kappa.html),
[GS98hz.gz](https://www.cita.utoronto.ca/~boothroy/data/kappa/GS98hz.gz).
The LLNL endpoint was unavailable during import.

Ten hydrogen fractions (0, 0.1, 0.2, 0.35, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98)
times 52 temperatures (log T=4.00..7.10) times 19 densities (log R=-8..1).
All 9,880 values are original source cells; none were interpolated or filled.
The full OPAL file is ragged, with `9.999` and blanks indicating missing
values. This import selects only the fully populated rectangle. The upper
temperature is 12.589 million K, **not** the original file's 8.70 dex limit.
The high-temperature, high-density corner and temperatures below log T=4
are deliberately absent. This is radiative opacity, without conduction.

Reproduce with:

```
curl -fLO https://www.cita.utoronto.ca/~boothroy/data/kappa/GS98hz.gz
python3 scripts/import_opal.py GS98hz.gz data/opacity/opal_gs98_z020.dat
```

SHA-256 of the downloaded gzip:
`10b526412dfa242b3759643bd6bbfb5aa43b5b70c285e898abc4504fbe973f71`.
SHA-256 of uncompressed `GS98hz`, checked by the importer:
`dbb0ebb9a231e30201f7347f2033c72783ed8acc1f6abcaa57bb89f2df0b2ca3`.

## `aesopus21_gs98_z020.dat`

ÆSOPUS 2.1 Rosseland mean **gas** opacity, GS98, nominal Z=0.020, from
[Marigo et al. (2024), ApJ 976, 39](https://arxiv.org/abs/2409.10905),
with the preceding method papers Marigo & Aringer (2009), A&A 508, 1539,
and Marigo et al. (2022), ApJ 940, 129. Retrieved 2026-09-07 from the
[authors' solar-mixture distribution](https://stev.oapd.inaf.it/aesopus_2.1/tables/SOLAR/),
[03-GS98.zip](https://stev.oapd.inaf.it/aesopus_2.1/tables/SOLAR/03-GS98.zip).
The selected files include molecular gas opacity and pressure broadening,
but **no grain opacity**. The publication also discusses condensate tables;
that does not make this particular gas archive include grains.

Ten X values (the same as OPAL), 211 temperatures (log T=2..4.5), and 71
densities (log R=-8..6 in steps of 0.2): **149,810 unmodified source cells**.
The importer joins each X's `OPALZ_lowT`, `OPALZ_highT`, and `OPALZ_highR`
files without interpolation or filling. Temperature steps are 0.01 up to
log T=3.7 and 0.02 above it. Archive filenames retain `aesopus2.0_gasbroad`,
but every selected header identifies AESOPUS 2.1. The lowT header reports 91
rows while actually containing 90 (2.00..2.89); the importer validates the
complete actual grid. Some source actual-Z fields print 0.02000001 while
Zref is 0.02000000. The importer checks these two allowed values and labels
all planes by nominal Zref; it does not alter their opacity cells.

Reproduce with:

```
curl -fLO https://stev.oapd.inaf.it/aesopus_2.1/tables/SOLAR/03-GS98.zip
python3 scripts/import_aesopus.py 03-GS98.zip data/opacity/aesopus21_gs98_z020.dat
```

Archive SHA-256, checked before import:
`9379a46e89b0e2c134d018986e8499e27cf6b45038de493f57b2f766ffaaa830`.
The local curl certificate-chain check failed during retrieval; the public
archive was downloaded with a request-specific `-k`, without credentials.

## `tops_gs98_x070_z020_{low,high}.dat`

Radiative Rosseland means returned by the [LANL TOPS service](https://aphysics2.lanl.gov/opacity/lanl),
selecting its `new`/ATOMIC library on 2026-09-07. The returned output does
not identify a unique dated atomic-data release; do not infer one from the
date of a service documentation page. The request uses X=.7, Y=.28, Z=.02
and the 19 GS98 metal mass fractions printed in the OPAL source. This import
supports **only this composition**, without conduction or molecules. The
cool envelope uses AESOPUS gas opacity.

The versioned `sources/tops_gs98_x070_z020.request.json` records the full
request; `sources/tops_gs98_x070_z020.txt` preserves the returned numeric
text, including the normalized 21-element mixture and density warnings.
Requested range: T=.002..2 keV on 50 native temperatures, rho=1e-10..1e4
g/cm³ on 71 log-spaced densities, grey means, plasma-frequency correction
enabled. TOPS silently substitutes a maximum available density at some
low-temperature requested points **but lists those substitutions in its
warnings**. Every warned pair is excluded from our tables.

Two complete rectangles retain original, un-clamped cells:

| File | Native T (keV) | rho (g/cm³) | Cells |
|---|---|---|---:|
| low | .002..2 | 1e-10..251.19 | 50 × 63 = 3150 |
| high | .025..2 | 1e-10..1e4 | 36 × 71 = 2556 |

The rectangles overlap; their cell counts are not unique-source counts.
Some valid ragged points are deliberately omitted. The first temperature
is 23209 K, the high rectangle begins at 290113 K, and the upper temperature
is 23.209 million K. Native printed densities are retained, including their
finite display precision. Unit conversion and logarithms are the only
changes to numeric values; no resampling to log R takes place.

`TopsOpacity` blends the rectangles over log T=5.6..5.7 using the same
smooth log-opacity rule as the low/hot blend. Both must cover an overlap
query, so densities above 251.19 are available only above log T=5.7.
The stellar driver blends AESOPUS to TOPS over log T=4.4..4.5.

Reproduce offline with:

```
python3 scripts/import_tops.py data/opacity/sources/tops_gs98_x070_z020.txt data/opacity
```

SHA-256 of the versioned text, checked by the importer:
`6e99ef76f5409926168a7defe7205c73d12e42166d8a388daffe3a4d44542344`.
The original HTML response had SHA-256
`db0a932ed4921175590264c394464cbe1710ecfbe9216c972b3449939e9d9451`.
The numeric text is obtained from its opacity-result preformatted block;
the source warning and composition blocks are retained for validation.
Both imported files were reproduced byte for byte. These external data
are not relicensed under ember's MIT code licence; attribute LANL TOPS/ATOMIC
and retain this provenance when using them.

Service reproduction notes: prepare the form calculation with `/submit`
before fetching `/results`; the latter can return a previously prepared
calculation despite receiving new form parameters. Validate the returned
composition and dimensions, never just HTTP success. `mixname` was left
blank. The `userid` and calculation identifier in the archived request are
the public form defaults, not credentials. Local retrieval needed a
request-specific curl `-k` for a certificate-chain failure.

## Interpolation and blending

All sources use `TabulatedOpacity`: monotone Hermite interpolation in log R
for Ferguson/AESOPUS/OPAL, or native log rho for TOPS,
then log T, with derivatives through active slope limiters, and linear
interpolation in hydrogen fraction. They reject a composition whose Z
differs from the table's Z by more than 1e-10. The metal pattern is fixed by
the source table; matching X and Z does not change GS98 into another mixture.

`BlendedOpacity` accepts explicit join temperatures (defaults 4.0 and 4.5), using
`ln k = (1-w) ln k_low + w ln k_high`, `w=u²(3-2u)` and
`u=(log10 T-start)/(end-start)`. Its temperature derivative includes
`dw/dlnT * (ln k_high-ln k_low)`. Density bounds inside the overlap are the
intersection of both sources' bounds. Outside the overlap only the selected
source is evaluated. Source errors are never caught as permission to switch
tables or extrapolate. The native log-rho mode omits the -3*dlnk/dlnrho
temperature-coordinate correction needed for log R. Ferguson/OPAL retains a low-temperature X range
0..0.9; AESOPUS/OPAL covers X=0..0.98 throughout.

The stellar driver supports AESOPUS/OPAL or AESOPUS/TOPS and an explicit blend over
log T=4.4..4.5 (25,119..31,623 K), retaining dense gas opacity below that
overlap. The older Ferguson/OPAL join remains tested. Grey integrations at
Teff=2500, 3000 and 4000 K, g=2e5 cm/s² now succeed with AESOPUS. Its
extended density domain is real source coverage, not an opacity clamp.

**Remaining coverage limitations:** OPAL still stops at log R=1, including
within the blend. The old ionized-EOS 0.1 Msun trial still fails there. The
CMS19/TOPS 0.1 Msun model now converges within strict source coverage; TOPS
has only one composition and cannot follow hydrogen depletion. See
[`docs/EQUILIBRIUM.md`](../../docs/EQUILIBRIUM.md).

Decimal grid edges may return a few ulps outside after conversion through
(T,rho). Only roundoff within `16*epsilon*max(1,abs(lo),abs(hi))` dex is
snapped to an axis endpoint. This does not admit physical extrapolation;
tests reject even 1e-10 dex outside the new table.

## Format

```
NX NT NR   <title>
<NR log density-coordinate values>
<NT log T values>
for each X:  X  Z
             NT rows of NR values (log10 kappa)
```

The constructor selects `DensityAxis::logR` (default) or `logRho`; the
calling source wrapper must match the file. A single X plane requires an
exact X match. Axes must be finite and strictly increasing; all X planes must have the same
Z. Missing-cell sentinels, truncated input, trailing tokens, invalid opacity
values and inconsistent compositions are errors.
