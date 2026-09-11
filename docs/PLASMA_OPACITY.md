# Plasma cutoff and radiative transport

The TOPS source averages require an additional check before the dense opacity
tables can be used for stellar evolution. The current stellar calculation
retains its original opacity inputs. No correction described here has been
installed in the evolution code.

## Averaging convention

For a prescribed temperature gradient, the diffusive radiative flux depends on
the integral of the inverse monochromatic opacity. In the convention of
[Colgan et al. (2016), equation 1](https://arxiv.org/abs/1601.01005),

\[
\frac{1}{\kappa_R} =
\frac{\int_0^\infty n_\nu^3\kappa_\nu^{-1}
 (\partial B_\nu/\partial T)\,d\nu}
 {\int_0^\infty (\partial B_\nu/\partial T)\,d\nu}.
\]

Here `n` is the refractive index. Section 2.8 of that paper describes the web
tables as using `n=1` above the plasma frequency and zero below it. The
[TOPS documentation](https://aphysics2.lanl.gov/static/opacdocs/opac-faq.html)
also describes the exclusion of radiation below the plasma frequency.

Independent integration of four returned spectra finds that the web service's
reported Rosseland means instead agree with an average whose numerator and
denominator are both restricted to frequencies above the cutoff. Write
`f_R` for the fraction of the full Rosseland weight above that cutoff. Within
this step approximation, the opacity appropriate to the equation above is
the conditional mean divided by `f_R`.

The cutoff is inferred from the reported Rosseland mean. Using that same
cutoff independently predicts the reported Planck absorption mean to within
0.02867% in the controls below. Switching the service's plasma option leaves
every returned monochromatic value unchanged. Integration with the cutoff
disabled reproduces the reported Rosseland and Planck means to within
0.003386%. These checks distinguish the averaging operation from a change
in the atomic spectra.

## Measured source controls

The core controls use X=0 and Z=0.02 near the final plotted star's central
temperature and density; they are not exact samples of its composition.
The denser controls use X=0.0375 and Z=0.02. Opacities are in cm²/g.

| Control | Temperature (MK) | Density (g/cm³) | Reported conditional Rosseland mean | `f_R` | Rosseland mean with full normalization |
| --- | ---: | ---: | ---: | ---: | ---: |
| Near the core | 11.60 | 3981 | 10.47 | 0.9749 | 10.74 |
| Near the core | 11.60 | 4000 | 10.50 | 0.9747 | 10.77 |
| Dense material | 2.611 | 1.874e5 | 29.07 | 5.084e-13 | 5.718e13 |
| Dense material | 2.611 | 1.875e5 | 28.85 | 4.645e-13 | 6.211e13 |

The very small `f_R` in the dense controls means that almost none of the
thermal radiation can propagate in the step approximation. A small conditional
mean alone would give a misleadingly large radiative heat flux there.

The source cutoff also moves between discrete frequency points. A narrow
density scan resolves jumps up to 0.7637% in the reported conditional mean.
Refining only the density interpolation does not remove this source effect.

## Effect on the plotted star

A fixed-profile estimate bounds the free-electron density by fully ionizing
the actual hydrogen and helium and using Z/A≤1/2 for the GS98 metals. The
classical plasma cutoff then removes at most 2.455% of the Rosseland weight
anywhere in the 3.848-trillion-year structure. Correcting only the normalization
would increase the radiative opacity by at most 2.517%. Including the star's
existing electron conduction reduces the maximum increase in the combined
opacity to 0.7555%.

This is a fixed-structure estimate for a continuous step cutoff. It excludes
source frequency-bin rounding and the effect of a frequency-dependent
refractive index. A separate collisionless-index integration near the core
gives an opacity of 11.50–11.54 cm²/g, showing that the choice of refractive
index also deserves attention. Neither estimate measures a change in stellar
lifetime or replaces a complete dielectric treatment.

## Comparison with the source electron density

The reported mean free-electron count and normalized composition also give an
independent estimate of the classical plasma frequency. Using the standard
helium atomic weight, the inferred cutoff agrees with this estimate within
0.04185% near the core. Across the two dense states, the classical cutoff
changes smoothly from 39.90 to 39.91 in units of photon energy divided by kT,
whereas the cutoff inferred from the reported mean changes from 39.90 to 40.00.
Direct integration at the continuous classical cutoff removes that discrete
frequency jump. This comparison uses rounded source composition and electron
counts; it does not include relativistic or collisional dispersion.

Conduction strongly limits the importance of this ambiguity in the dense
controls. With the selected conduction prescription, removing radiative heat
transport entirely increases their combined opacity by at most 0.0002651%.
Near the core, the step-cutoff normalization changes the combined opacity by
about 0.7421%; the collisionless refractive-index diagnostic changes it by
about 2.700%. These four fixed-state comparisons do not replace the actual
profile calculation or establish a bound for the whole added table domain.
[Core comparison](results/tops_core_classical_cutoff_v1.json),
[dense comparison](results/tops_dense_classical_cutoff_v1.json).

## Interpolation with conduction included

Across the dense opacity source nodes and independent composition/density
queries, the largest combined-opacity interpolation difference is 0.07885%.
For material temperatures from 1 MK to 20 MK it is 0.03639%, compared with a
radiative-only difference of 0.7176%. These comparisons include the selected
conduction prescription at helium-3 fractions 0 and 0.12. Conduction therefore
reduces the effect of the radiative interpolation errors substantially in these
tested states. Further density refinement alone would have limited benefit.

This does not remove the plasma-normalization question. In some of the tested
states, radiative transport still contributes enough that removing it entirely
would change the combined opacity by much more than the interpolation error.
The diagnostic leaves all stellar inputs unchanged.
[Combined transport comparison](results/tops_density_fine_transport_v1.json).

## Reproduction and remaining work

`scripts/audit_tops_spectral_means.py` checks the source receipts, integrates
the spectra interval by interval and compares two quadrature orders. Constant
and power-law opacity controls verify its normalization independently.
`scripts/audit_profile_plasma_normalization.py` re-evaluates the saved profile's
opacity and conduction with its selected physical tables.

The records are [dense spectra](results/tops_dense_spectral_means_v1.json),
[core spectra](results/tops_core_spectral_means_v1.json),
[density scan](results/tops_density_structure_v1.json) and
[fixed stellar profile](results/plasma_normalization_profile_3848gyr_v1.json).
Raw requests and spectra remain local with recorded hashes.

Before selecting dense tables, the radiative transport prescription must
consistently combine the monochromatic opacity, plasma dispersion and spectral
normalization. Missing source coverage and substituted densities remain
excluded. The existing conduction calculation remains a separate transport
channel.
