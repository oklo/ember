# Cluster dynamics and numerical checks

The simulations use [PeTar](https://arxiv.org/abs/2006.16560), FDPS and SDAR
regularized close encounters, SSE/BSE stellar evolution, and galpy's
MWPotential2014 Galactic potential. Gravity is unsoftened. No distant escaper
is removed. Ember supplies low-mass luminosities for the synthetic catalogues.

Each model starts with 3000 components in a virialized Plummer distribution,
truncated at 20 initial half-mass radii. The models omit residual gas and initial
mass segregation. Their sizes and binary populations are controlled assumptions,
not inferred birth conditions. The [17-model comparison](cluster-dynamics.md)
gives the tested initial radii, realizations and outcomes.

The input mass spectrum has dN/dm proportional to m^−0.4 below 0.08 M☉,
m^−1.3 between 0.08 and 0.5 M☉, and m^−2.3 between 0.5 and 50 M☉.
Its lower limit is 0.01 M☉. The substellar normalization is set independently
to the declared birth ratio, usually 0.63, with a Pleiades alternative of 0.55.
Continuity at 0.08 M☉ is not imposed as an additional constraint.

Models include single-star controls and a similar-mass binary prescription:
binary fractions of 15% for brown dwarfs and 30% for stars, no primordial pairs
straddling 0.08 M☉, and log-uniform separations of 0.2–20 AU and 0.2–200 AU,
respectively. Two random-pairing stress tests use a 25% system binary fraction;
they are retained but excluded from the preferred comparisons.

Cluster end ages are 650 Myr for the Hyades and 125 Myr for the Pleiades, with
Z = 0.02307 and 0.01648. Present positions and velocities are integrated backward
in the Galactic potential to obtain starting orbits. The adopted transformation
uses a solar Galactic radius of 8 kpc, circular speed of 220 km/s and the
Schönrich solar motion. Mock sky positions condition the surviving cluster center
on its observed location while retaining simulated relative positions and tails.

The retained controls check orbit round trips, tree steps and encounter grouping.
Two Hyades trajectories continue from saved states after step-size checks;
their costs include both segments. The seed-43 five-Myr comparison changes the
half-mass radius by 1.370e−6 fractionally without changing aperture number ratios.
This local numerical check does not imply identical long-term chaotic trajectories.
Initial and surviving particle identities, masses, endpoint ages and count
accounting are checked in the validation records. Missing initial identities
are assessed against stellar/binary event logs rather than assumed to be escapers.

[Build and numerical records](../provenance/dynamics),
[ensemble validation](../provenance/ensemble/validation_receipt.json),
and [local-data guide](../REPRODUCING.md) identify the retained evidence.
The current figures use all 17 completed models; the separate eight-model pilot
figure has been retired.
