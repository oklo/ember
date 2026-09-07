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

## Format

```
NX NT NR   <title>
<NR log R values>
<NT log T values>
for each X:  X  Z
             NT rows of NR values (log10 kappa)
```
