# Numerical lineage

**Henyey FORTRAN IV → Bodenheimer FORTRAN 77 → Astra/Opus Ember C++**

Ember follows the Henyey–Bodenheimer tradition of solving stellar structure as
a coupled boundary-value problem. The defining numerical idea is to linearize
the equations about a trial star, eliminate outward through the mass mesh,
then substitute inward for the correction. This keeps the solve proportional
to the number of mesh points. Bodenheimer's
[memoir of Henyey](https://www.nationalacademies.org/read/4961/chapter/10)
describes the development and influence of this method.

The local companion project, `../low_mass_stars`, reconstructs the FORTRAN IV
calculation and a later F77 low-mass program. Its README explicitly identifies
the FORTRAN IV deck as a modern implementation of the published method.
The F77 reconstruction also contains substitutions; it is not the archived
program used for Laughlin, Bodenheimer and Adams (1997). The
[F77 comparison](F77_LBA97_COMPARISON.md) records the substantial differences
from that paper's 0.1-solar-mass milestones. These reconstructions make the
method accessible, but do not establish an unbroken source-code ancestry.

Astra/Opus Ember is a C++ implementation developed through this research
project. Its structural solve is in [henyey.cpp](../src/henyey.cpp), with
[algorithm and numerical checks](HENYEY.md) documented separately. Material
physics, composition evolution and time-step control extend that structural
solve. The research record keeps the inputs and numerical choices associated
with each study; later source changes do not retroactively validate older
tracks.

The historical lesson for further work is economy: establish the governing
balances analytically, use interpolation and controlled extrapolation to map
the problem, and spend detailed integrations on comparisons that can change
the scientific conclusion.
