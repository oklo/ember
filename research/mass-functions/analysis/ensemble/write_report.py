"""Write measured results from the validated endpoint tables."""
from pathlib import Path
import json
import pandas as pd

W=Path(__file__).resolve().parent
D=W.parents[1]/'docs/research/independent'
e=pd.read_csv(W/'structural_endpoints.csv')
j=pd.read_csv(W/'joint_checks.csv')
o=pd.read_csv(W/'observed_matched_counts.csv')
v=json.loads((W/'validation_receipt.json').read_text())
f=lambda x:f'{x:.4g}'
pending=v['pending_runs']
new=e[e.workspace==W.name]
status=(f"{len(e)} full-age simulations are included, of which {len(new)} are new. "
        + ('All planned runs have finished.' if not pending else
           f"{len(pending)} additional Hyades realization remains in progress and is excluded from every endpoint comparison."))
obsrows=[]
for cl,sc,label in [('hyades','all_union','Hyades, all associated candidates'),
                   ('hyades','central_matched','Hyades, within 10 pc'),
                   ('hyades','ir_matched','Hyades, matched infrared region'),
                   ('pleiades','all_union','Pleiades, all associated candidates'),
                   ('pleiades','ir_matched','Pleiades, central DANCe sample')]:
    r=o[(o.cluster==cl)&(o.scenario==sc)].iloc[0]
    obsrows.append(f"| {label} | {f(r.N_BD)} | {f(r.N_star)} | {f(r.ratio)} | {f(r.ratio_lo)}–{f(r.ratio_hi)} |")
dynrows=[]
for r in new.itertuples():
    dynrows.append(f"| {r.cluster.capitalize()} | {f(r.initial_rh_pc)} | {r.seed} | {f(r.expected_birth_ratio)} | {f(r.aperture_mass_excluding_BH_NS_Msun)} | {f(r.aperture_halfmass_pc)} |")
sel=j[(j.structural_screen_pass)&~j.random_pair_stress&(j.faint_extension!='continued_slope')]
pr=[]
for birth in [.63,.55]:
    x=sel[(sel.cluster=='pleiades')&(sel.scenario=='ir_matched')&(sel.expected_birth_ratio==birth)]
    pr.append((f(x.model_ratio.min()),f(x.model_ratio.max())))
h=sel[(sel.cluster=='hyades')&sel.scenario.str.startswith('wide')]
hr=(f(h.model_ratio.min()),f(h.model_ratio.max()))
report=f'''# Cluster dissolution and the low-mass population

**The additional calculations strengthen the case that unequal survey reach
inflates the apparent brown-dwarf deficit. They still do not establish that
dissolution and incompleteness explain the entire difference, or that the
clusters formed different mass functions.** {status}

All associated members and tidal-tail candidates remain the scientific target.
Smaller apertures below are diagnostic comparisons, not a redefinition of the
population. No escaper is deleted from the simulations.

![Structure, spatial distributions and number ratios](CLUSTER_DISSOLUTION_FOLLOWUP_20261009.png)

The panels compare fixed-aperture cluster mass and size, reference-star spatial
distributions, and brown-dwarf/reference-star ratios. Boxes are broad structural
checks, not confidence regions. Each model span varies the stated atmosphere
and survey choices; it is not a completeness posterior. Grey shading shows the
field component-ratio sensitivity; simulated catalogues blend unresolved pairs.
Random-pairing stress tests are preserved in the tables but omitted from the plot.

## The observational comparison

Define R = N(0.03–0.08 solar masses)/N(0.10–0.20 solar masses). These are fixed
comparison intervals, not classifications of hydrogen burning for every object.
The counts use the existing Ember luminosity-to-mass inference, including its
age and atmosphere choices. Intervals add approximate counting uncertainty;
they do not include unknown contamination, completeness or binary distributions.

| Sample | Brown dwarfs | Reference stars | R | Approximate 68% interval |
|---|---:|---:|---:|---:|
{chr(10).join(obsrows)}

The DANCe comparison uses membership probability above 0.84, angular separation
below 3 degrees, and Ks < 17. This probability-only cut is slightly stricter than
the paper's membership-plus-uncertainty candidate definition. Its published
classifier tests give 90% recovery and 4.3% contamination at that threshold.
Their factor 0.9404 is used only as an absolute-count sensitivity; assuming it
acts equally on both mass intervals leaves R unchanged. It is not a calibration
of brown-dwarf recovery in unobserved tails.
[Olivares et al. 2018](https://arxiv.org/abs/1804.08360),
[Bouy et al. 2015](https://arxiv.org/abs/1502.03728).

Missing Hyades J/K photometry was filled using exact Gaia identifiers and the
Gaia–2MASS crossmatch. Missing values no longer silently fail the brightness cut.
Bolometric luminosities and masses were unchanged. The additional measurements
also improve the empirical infrared magnitude–luminosity relation used to observe
the simulations. Photometry coverage and query checksums are saved separately.
The Hyades infrared region uses the published sky bounds, an 18-pc radius and
J < 17.5; its heterogeneous member list is not a complete infrared census.
[Schneider et al. 2022](https://arxiv.org/abs/2203.11090),
[Schneider et al. 2024](https://arxiv.org/abs/2408.10112).

## What the new simulations establish

The simulations use the same unsoftened PeTar dynamics, stellar mass loss,
Galactic potential and similar-mass binary prescription as the
[first calculation](CLUSTER_DISSOLUTION_RESULTS_20261009.md). Each starts with
3000 components. Seeds 41–43 vary the initial realization. Initial half-mass
radii are 2.5 pc for the Hyades and 3 pc for the principal Pleiades ensemble.
Three additional Pleiades calculations at 4 pc use R = 0.55, the lower end of
the independently obtained field sensitivity range; this is not a fitted birth
ratio, nor a corresponding Hyades test. The realized birth ratios fluctuate
around their expected values and are recorded individually.
The field inference and its assumptions are described in the
[field comparison](FIELD_MASS_FUNCTION_20261009.md).

| Cluster | Initial half-mass radius (pc) | Seed | Expected birth R | Final aperture mass (solar masses) | Final half-mass radius (pc) |
|---|---:|---:|---:|---:|---:|
{chr(10).join(dynrows)}

Mass and half-mass radius refer to the same fixed aperture: 9 pc for the Hyades,
11.6 pc for the Pleiades. White dwarfs are included; black holes and neutron stars
are excluded. Broad mass ranges of 240–360 and 600–1000 solar masses, respectively,
were specified before examining the new endpoints. A model just outside a box
is not physically ruled out. Definitions and literature mass estimates differ.
[Hyades census](https://arxiv.org/abs/1105.6093),
[Pleiades mass estimates](https://academic.oup.com/mnras/article/525/2/2315/7233111).

**Pleiades:** structurally plausible models predict central DANCe ratios
{pr[0][0]}–{pr[0][1]} for expected birth R = 0.63 and {pr[1][0]}–{pr[1][1]}
for R = 0.55. The latter approach the observed 0.5019, but their reference-star
counts remain too large in the center. After the aggregate classifier factor,
they predict 317–359 reference stars versus about 275 observed. Other models
can come closer in absolute counts while missing the brown-dwarf ratio. The
wide reference-star distribution also extends farther than most simulations.
There is no satisfactory joint fit yet.

**Hyades:** a new realization has 333.3 solar masses inside 9 pc and half-mass
radius 3.336 pc, satisfying the broad structural check. Yet it predicts a wide
catalogue ratio of 0.3049–0.3514, versus 0.1346 observed within 200 pc. It also
contains too few wide reference stars while overpredicting their central count.
That cannot be repaired simply by declaring additional missed observed stars.
Across all completed structurally screened Hyades models the wide ratio range
is {hr[0]}–{hr[1]}. The realization with 232.5 solar masses is close to, but below,
the broad mass range. The stochastic structural spread matters.

**Multiplicity:** using the published DANCe near-equal-binary probabilities and
two assumed component flux ratios raises the central component estimate from
about 0.50–0.52 to 0.53–0.55. This is a modest effect. It does not recover unseen
low-mass companions or tails, and the component-corrected result must not be
compared directly with a model whose binaries remain blended.

![Normalized mass distributions in matched regions](CLUSTER_MASS_SHAPES_20261009.png)

The second figure checks distribution shape after normalizing to reference
stars. Absolute counts are checked separately above. Its upper plotted mass
is 0.20 solar masses to avoid interpreting the synthetic-light grid's 0.25-solar-mass
edge as a physical fall in the mass function. Model shading spans the completed
realizations and atmosphere choices, not a fitted probability interval.
The stellar distribution shapes also disagree: fitting the integrated reference
count alone would overlook this. The imposed stellar birth-mass slope has not
been inferred from these clusters.

A direct magnitude-count check includes DANCe candidates even when they lack
an assigned luminosity or mass. In the central sample, 13 of 334 candidates
with 13.5 <= Ks < 17 lack luminosities. Including them leaves a model excess in
most magnitude bins. For example, the lower-field models predict 67–75 sources
at Ks = 14.5–15 after the classifier adjustment, against 42 observed. The
disagreement therefore is not just loss of these unconverted candidates.
The empirical conversion from model luminosity to magnitude remains an
assumption. Counts and the full missing-luminosity list are saved.

## What remains decisive

The all-associated Pleiades sample need not be the debris of one bound cluster.
The Greater Pleiades Complex includes multiple coeval groups over about 600 pc
whose common origin is a giant molecular cloud. A single initially virialized
Plummer cluster is consequently a restricted hypothesis for that entire union.
[Boyle, Bouma & Mann 2025](https://arxiv.org/abs/2511.07533).

The next useful calculation is a joint fit to stellar counts and spatial
distributions with a calibrated survey selection, allowing an extended birth
population, gas loss and initial mass segregation where required. The brown-dwarf
counts should then test that fit. Unconstrained extra N-body runs or freely fitted
completeness factors would not resolve the ambiguity. Current results make an
intrinsic birth-mass-function claim premature.

## Validation and reproduction

All included endpoints reach 650 Myr or 125 Myr. Unique initial/surviving IDs,
total mass, and the agreement of spatial and mass-bin counts were checked.
The Hyades seed-43 continuation used a doubled global tree step after 200 Myr;
a 100–105 Myr control changed half-mass radius by 1.370e-6 fractionally, with
unchanged aperture number ratios. Binary grouping, unsoftened forces and close
encounter treatment were retained. This short test is not a proof of identical
long-term chaotic trajectories.

Completed simulations consumed {f(v['cpu_hours_completed'])} CPU-hours in total;
the new completed subset consumed {f(v['cpu_hours_new_completed'])}. These totals
include discarded work before saved-state continuations, exclude installation
and short probes, and {'exclude the still-running realization' if pending else 'include every planned full-age run'}.
They are CPU time, not elapsed campaign time.

Workspace: `out/cluster-dissolution-oct9-v2/`. Its README contains execution order.
`cluster_followup.sqlite` stores structural outcomes, counts, initial components,
final descendants and cost. `validation_receipt.json` identifies pending runs and
hashes every included endpoint. `joint_checks.csv`, `normalized_mass_shapes.csv`,
`dance_binary_sensitivity.csv` and the full draw tables preserve the comparisons.
No shared Ember source or stellar table was changed.
'''
(D/'CLUSTER_DISSOLUTION_FOLLOWUP_20261009.md').write_text(report)
print(status)
