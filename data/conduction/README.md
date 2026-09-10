# Ioffe electron-conduction data

**Storage policy:** new generated tables and raw source archives described below
remain local and are not included in a fresh clone. Generator code, source patches,
small specifications and provenance metadata are versioned. See
[the reproduction guide](../../docs/DATA_REPRODUCTION.md). Existing published
data history is retained. References to installed/archived files describe the
development machine unless explicitly stated otherwise.

The three source tables are from [the Ioffe conductivity archive](https://www.ioffe.ru/astro/conduct/condint.html), retrieved 2026-09-08. They supply **log10 thermal conductivity** in erg cm⁻¹ s⁻¹ K⁻¹, rather than opacity. The retained source variants are:

- `condtab21wd`: 2021 weakly damped Blouin correction (evolution default).
- `condtab21_I`: no Blouin correction; comparison control.
- `condtab21nd`: undamped Blouin correction; comparison control.

Original files are archived as deterministic gzip files in `sources/`, with raw and compressed SHA-256 checksums in `sources/manifest.json`. `import_conduction.py` preserves the original cells for charges 1, 2, 3, 4, 6, 8 and 12; their source mass numbers are 1, 4, 7, 9, 12, 16 and 24. Each plane contains 19 temperatures, log10 T=3..9, and 64 densities, log10 rho=-6..9.75. The webpage's stated upper density of 10⁹ is narrower than the actual file; the reader follows the file. No source cells are extrapolated or filled.

```
python3 scripts/import_conduction.py data/conduction/sources data/conduction
python3 scripts/verify_composition_data.py --extended
```

The format starts with `EMBER_CONDUCTIVITY 1 nZ nT nR label`, followed by temperature and density axes. Each charge/mass pair precedes its density-major rows of source log conductivity. These scientific data are not relicensed under ember's MIT code license.

The optional independent comparison uses the authors' [conduct21.f](https://www.ioffe.ru/astro/conduct/conduct21.f), SHA-256 `37cebbc6adb4567adbb65525dd9cd5fbd4dd4a5d6dc957c8aa4afc5f218fe027`. This source evaluates the conductivities directly and supplies the three H/He correction variants. Normal builds need neither it nor Fortran.

```
gfortran -O2 -std=legacy -ffixed-line-length-none /tmp/conduct21.f \
  scripts/conduction_reference_probe.f90 -o /tmp/conduction-reference
python3 scripts/generate_conduction_reference.py /tmp/conduction-reference \
  tests/data/conduction_reference.dat
```

The 30 independent source-code samples differ from ember's table interpolation by at most 0.895%. This tests interpolation and units, not the underlying physical uncertainty. See [CONDUCTION.md](../../docs/CONDUCTION.md) for transport and mixture assumptions. The authors request citation of [Potekhin, Pons & Page (2015)](https://arxiv.org/abs/1507.06186), together with relevant original work.
