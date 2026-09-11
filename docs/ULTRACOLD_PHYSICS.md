# Physics for very cold helium remnants

This document specifies calculations still needed for Ember's initially
0.1-solar-mass star through cooling past 100 K and disappearance under explicit
proton/nucleon-decay scenarios. It does not claim that a complete cooling trajectory or
an environmental heating model has been implemented. The immediate stellar
constraint is dense interior opacity; the latest computed model is at
3.890 trillion years. [Current calculation](COLD_REMNANT.md).

An effective temperature of 100 K does not imply a core at 100 K. Evaluate
material properties at the temperature, density and composition of each
layer. Compute an isolated reference track and then environmental cases,
including a remnant retained in a galaxy and one ejected from it. Ejection
from a stellar system need not remove the remnant from its dark-matter halo.

## Diffusion, grains and crystallization

Diffusion and element separation have priority because they can change the
hydrogen envelope and nuclear burning before the cold phase. Their treatment
must conserve each species and include the associated energy changes. A
helium core can retain a hydrogen atmosphere. At strong coupling, ions with
equal mass-to-charge ratios can separate through their Coulomb chemical
potentials; ordering elements by mass alone is insufficient.
[Beznogov & Yakovlev (2013)](https://arxiv.org/abs/1307.6060).

Grain formation must couple elemental depletion, gas chemistry, absorption,
scattering and settling. Surface metal abundance must follow diffusion and
accretion. A permanently fixed solar metal mixture is not an adequate
assumption for the cooling atmosphere. Existing optical calculations and
unfinished coupling are described in [GRAINS.md](GRAINS.md).

Determine solid–liquid coexistence from equal pressure, temperature and
chemical potentials, with latent heat and composition partitioning. Include
quantum heat capacities and conductive transport. Fully ionized Coulomb
models require a different treatment when the freezing region reaches
partially ionized or neutral material. Skye provides a useful fully ionized
comparison; quantum-liquid calculations are particularly relevant to helium.
[Jermyn et al. (2021)](https://arxiv.org/abs/2104.00691),
[Baiko & Yakovlev (2019)](https://arxiv.org/abs/1910.06771).

The absence of present-day room-temperature white dwarfs does not remove all
experimental or theoretical constraints on their material physics.
High-pressure helium melting calculations span 425–10000 K and 15 GPa–35 TPa.
Their applicable range must be checked against the actual envelope; this
range alone does not establish coverage at 100 K.
[Preising & Redmer (2019)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.184107).

## Nuclear heating at very low temperature

Pycnuclear reactions arise from quantum motion of dense ions. A thermal
reaction rate multiplied by an arbitrarily large screening factor is not a
controlled continuation into this regime. Evaluate residual hydrogen and
helium-3 reactions as well as possible helium burning, with explicit rate
uncertainties. A carbon–carbon rate cannot stand in for helium: two helium-4
nuclei do not make a stable burning product, and the three-body route matters.
Published treatments identify different burning regimes and sensitivity to
the arrangement of ions in mixtures.
[Gasques et al. (2005)](https://arxiv.org/abs/astro-ph/0506386),
[Yakovlev et al. (2006)](https://arxiv.org/abs/astro-ph/0608488).

For every candidate reaction, compare its integrated heating with photon
luminosity and its composition-change time with the cooling time. Low mass
suggests weaker density-driven burning than in massive white dwarfs, but this
has not yet been quantified for Ember. Black-dwarf supernova calculations
concern much more massive remnants and do not establish an endpoint for this
helium star. [Caplan (2020)](https://arxiv.org/abs/2008.02296).

## Accretion and encounters

Follow gas and solid interstellar objects separately, with evolving abundance,
velocity and size distributions. Material from a surviving planetary system
is another possible source. Accretion changes mass and composition as well as
energy. The gravitational scale is GM times accreted mass divided by radius;
only the portion retained below the radiating surface delays interior cooling.
Present-day X-ray detection of planetary-debris accretion demonstrates this
process, but does not measure a future interstellar-object flux.
[Cunningham et al. (2022)](https://arxiv.org/abs/2202.12903).

Collision rates alone do not determine tidal heating. For a perturber of
mass m, relative speed at infinity v and pericenter r, the point-mass encounter
cross section is

    sigma(<r) = pi r^2 [1 + 2 G (M + m) / (r v^2)].

Integrate the differential encounter rate against energy deposited in the
star. Near misses are more numerous than collisions, while the tidal coupling
declines steeply with separation. Calculate both dependencies. Model rare
events and the subsequent cooling response, since a mean heating rate can
hide long quiet intervals and short reheating episodes.

The specific proposed study concerns pairs that remain unbound after a close
passage. Define the separation explicitly, for example r/(R1 + R2) just above
unity; distinguish this from the impact parameter at infinity and from the
inverse separation convention often denoted beta. Check that the outgoing
orbital energy remains positive. Separate reversible deformation, surviving
oscillations, spin, irreversible heat and energy carried away by ejecta.
Crystallized stars require the appropriate elastic response and dissipation.

A representative WD–WD SPH encounter study explicitly restricts every system
to initially bound orbits after capture. It supports the distinction raised
here, but does not prove that a systematic unbound heating survey is absent.
That literature search remains open. Likewise, heating in compact helium-WD
binaries supplies useful comparisons without determining heating from an
isolated hyperbolic passage.
[Lorén-Aguilar et al. (2010)](https://academic.oup.com/mnras/article/406/4/2749/1022010),
[Fuller & Lai (2013)](https://arxiv.org/abs/1211.0624).

For numerical encounter calculations, retained heat must exceed the measured
error from isolated-star evolution, particle relaxation and artificial
viscosity. Check orbital plus stellar energy conservation and recover the
weak-tide limit before interpreting heating near contact. Couple resolved
heat deposition back to a stellar cooling calculation rather than identifying
all lost orbital energy with immediate thermal energy.

## Dark-matter capture and heating

The long-term framework begins with
[Adams & Laughlin (1997)](https://arxiv.org/abs/astro-ph/9701131).
Modern capture calculations include realistic stellar profiles, degenerate
electrons, finite temperature, nuclear form factors and multiple scattering.
These effects make capture depend on particle mass, interaction channel and
the actual stellar structure.
[Bell et al. (2021)](https://arxiv.org/abs/2104.14367).
Calculations for heavy dark matter also include phonon emission and absorption
in crystallized cores; their thermalization times can differ greatly from
previous approximations. This progress concerns the conditional capture
problem and does not establish the particle's existence.
[Bell et al. (2024)](https://arxiv.org/abs/2404.16272).

Experimental bounds must be applied to the same mass and interaction used
in a stellar calculation. For example, LZ's published standard
spin-independent result gives a 90% upper limit of 2.2e-48 cm^2 at
40 GeV/c^2. The light-particle search reports no significant dark-matter excess
and evidence for solar-neutrino scattering. These are different signals.
[LZ (2025)](https://arxiv.org/abs/2410.17036),
[LZ light-particle search](https://lz.lbl.gov/wp-content/uploads/sites/6/2025/12/LZ_Paper_Preprint_WS2025_v6.0_20251208.pdf).

The September 2026 extended-recoil search reports one unusual event, with
global significance 2.6 sigma against the background-only hypothesis. This
does not establish WIMP dark matter or a heating rate for a future remnant.
[LZ (2026)](https://arxiv.org/abs/2609.02823).

For self-annihilating captured particles, evolve their number N with capture
rate C, annihilation coefficient A and evaporation coefficient E:

    dN/dt = C - A N^2 - E N.
    L_ann = f_heat m_chi c^2 A N^2.

The often-used L_ann = f_heat m_chi c^2 C requires capture–annihilation
equilibrium and negligible evaporation. The deposited fraction f_heat must
exclude escaping products. Nonannihilating particles do not provide this
rest-mass heating source. Geometric capture is an upper limit, not a default
capture probability. The dark-matter density and velocity distribution must
evolve with the remnant's orbit and halo. These equations specify planned
conditional models; no dark-matter heating is currently installed in Ember.

## Nucleon decay and final disappearance

The continuation through decay must specify a lifetime and branching fractions
for each allowed nucleon process. Experiments constrain partial lifetimes,
not one universal decay time: for example, the Super-Kamiokande limit for
p to positron plus neutral pion is tau/B > 2.4e34 yr at 90% confidence.
The 2026 searches of charged-antilepton plus two-neutral-pion channels also
report results consistent with background. These observations do not establish
a decay process. [Takenaka et al. (2020)](https://arxiv.org/abs/2010.16098),
[Super-Kamiokande (2026)](https://arxiv.org/abs/2604.10975).

Baryogenesis does not require ordinary proton decay. Leptogenesis followed
by high-temperature electroweak baryon-number violation is possible without
grand unification. Keep a stable-proton alternative instead of inferring a
finite proton lifetime from the cosmic baryon asymmetry.
[Fukugita & Yanagida (1986)](https://doi.org/10.1016/0370-2693(86)91126-3).

Gravitationally induced decay through virtual black holes is an additional
scenario. Distinguish it from collective tunneling of a portion of a
degenerate star into a black-hole configuration. The amplitudes, suppression
scales and possible protecting symmetries require a specified quantum-gravity
model; neither the observed baryon asymmetry nor a laboratory lower limit
fixes them. [Adams et al. (1998)](https://arxiv.org/abs/astro-ph/9808250),
[Adams et al. (2001)](https://arxiv.org/abs/hep-ph/0009154).

For the latter paper's Planck-scale interaction, the dimensional estimate is
tau_p ~ m_p^-1 (M_Pl/m_p)^4 ~ 10^45 yr in natural units (equation 2). Treat it
as a conditional comparison, not a known lifetime. Collective stellar
tunneling has a different exponential suppression and cannot be represented
by assigning that same rate without a physical argument.

Follow the nuclear cascade after a decay, including unstable daughter nuclei,
spallation, capture and possible density-driven fusion. Track rest-mass loss,
changes in nuclear binding energy, neutrino and photon escape, and deposition
of charged-particle energy. Begin with coarse composition groups and bracketed
retained fractions. Resolve column-dependent deposition or daughter reactions
when they change the phase order or heating appreciably; the
[approximation remit](APPROXIMATION_REMIT.md) governs this progression.
Heating and chemical evolution can alter its
phase and radius before most nucleons have disappeared.

The continuation must allow loss of degeneracy, neutral or molecular matter,
material supported mainly by interatomic forces, and optical thinning at both
decay-product and thermal wavelengths. These are distinct transitions.
Beyond the continuum regime, follow remaining particles or a probability
distribution. Distinguish a vanishing bound object from the last decay of
formerly associated baryons. No fixed-mass white-dwarf sequence, extrapolated
blackbody atmosphere or artificial mass floor establishes this endpoint.
[Adams & Laughlin (1997), section IV](https://arxiv.org/abs/astro-ph/9701131).

The full timeline must label computed events, literature estimates and
conditional outcomes separately. Parameterize the decay era by the assumed
lifetime; do not assign a finite date to a stable-proton case. A benchmark
with constant independent lifetimes has exponential expected particle number
and an exactly calculable last-event distribution, useful for validating the
discrete continuation. The physical model still needs the nuclear and
environmental effects above.

## Order of implementation and evidence

First resolve the current density coverage, finish atmosphere acceptance and
continue hydrogen exhaustion. Develop diffusion and coupled grain chemistry
before their effects become important along the track. Validate solid and
partially ionized matter in the regions the cooling star actually visits.
Then evaluate nuclear and external heating against its calculated luminosity,
using bounds to identify which detailed calculations are necessary.

For each result, keep numerical error, uncertainty in material properties and
assumptions about the future environment explicit. Report families of cooling
histories and event probabilities when those assumptions permit different
outcomes. A single extrapolated present-day galactic environment cannot supply
a unique temperature floor or age at very late times.
