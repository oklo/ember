# Physics and remaining work

Ember uses one evolution program from Hayashi contraction through hydrogen
burning and helium-white-dwarf cooling. The aim is to extend that calculation
to extreme cold and conditional nucleon-decay scenarios.

## Implemented

- Deuterium, pp and CN burning, screening and plasma-neutrino losses.
- Composition-dependent free-energy EOS, radiative opacity and electron conduction.
- Non-gray gas atmospheres, including molecular collision-induced absorption.
- Hot H/He/metal diffusion and conservative composition and energy transport.
- Ledoux convection with instantaneous or finite-rate mixing.
- Adaptive time control, checked Richardson extrapolation and conservation audits.
- A tested weakly quantum liquid-ion free-energy correction.

The current stellar configuration uses instantaneous mixing within each
connected convective region. Metals share a GS98 mixture response rather than
independent element velocities. The nuclear network is not full CNO.

## Remaining work

1. Apply the integrated envelope and physically supported EOS joins throughout
   a fresh Hayashi-to-cooling calculation. Validate surface composition, residual
   burning and the helium-3 shell against resolution and physical alternatives.
2. Extend molecular and grain physics, cool diffusion, dense-matter thermodynamics
   and conductive transport as the star enters their relevant regimes.
3. Add helium crystallization, latent heat and strongly quantum mixtures.
4. At sufficiently low luminosity, compare isolation with accretion, tidal encounters
   and conditional dark-matter heating. Treat nucleon decay with an explicit assumed
   lifetime and retain stable baryons as an alternative.
5. Extend the mass/composition survey across the stellar–brown-dwarf boundary.

The integrated cooling envelope and some source tables remain local research
components. The final physics has not yet been validated in one uninterrupted
calculation. No extreme-cold cooling age or proton-decay track is claimed.

[Atmospheres](ATMOSPHERE.md) · [Dense EOS](DENSE_EOS.md) ·
[Nuclear network](NUCLEAR.md) · [Grains](GRAINS.md) ·
[Very cold physics](ULTRACOLD_PHYSICS.md) ·
[Approximation choices](APPROXIMATION_REMIT.md) ·
[Configuration](LIFETIME_DRIVER.md) · [Data](DATA_REPRODUCTION.md)
