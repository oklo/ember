"""Compare observed cluster mass-distribution shapes on a common normalization.

Normalize each saved sensitivity realization by its own count at 0.10--0.20
solar masses. This avoids letting the unequal faint-survey depths set the
normalization, and avoids the upper model boundary at 0.25 solar masses.
It does not establish completeness in the reference interval or elsewhere.
"""
from pathlib import Path
import csv
import hashlib
import json
import shutil

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

W = Path(__file__).resolve().parent
DOC = W.parents[1] / 'docs/research/independent'
REFERENCE = (0.10, 0.20)
# Explicit normalization boundaries permit an exact integral check. Narrower
# bins at the upper end retain the information in the earlier figure.
BINS = np.array([.008, .010, .014, .020, .028, .040, .056, .080,
                 .10, np.sqrt(.10 * .20), .20, .25])
WIDTH = np.diff(np.log10(BINS))
REF_BINS = (BINS[:-1] >= REFERENCE[0]) & (BINS[1:] <= REFERENCE[1])
source = json.loads((W / 'observed_mass_function_receipt.json').read_text())
receipt = dict(
    normalization_interval_Msun=list(REFERENCE),
    definition='(dN/dlog10M) / N(0.10 <= M < 0.20 Msun), separately for every saved realization',
    reason='Common stellar interval, above uncertain substellar coverage and below the upper model boundary. No claim of completeness.',
    population='All selected candidates with a model-supported mass, including identified tails and associations; no core-only cut.',
    uncertainty='16th--84th percentiles of the existing age/cloud/initial-state/luminosity sensitivity draws, including normalization variation. Not counting errors or a complete posterior.',
    missing_corrections=['detection completeness', 'membership contamination', 'unresolved binaries'],
    source_receipt='observed_mass_function_receipt.json',
    bins_Msun=BINS.tolist(), clusters={})

plt.rcParams.update({'font.size': 12, 'axes.spines.top': False,
                     'axes.spines.right': False, 'savefig.dpi': 220})
fig, ax = plt.subplots(figsize=(9.2, 6.0))
ax.axvspan(*REFERENCE, color='#5d6065', alpha=.055, zorder=0)
rows = []
reference_counts = {}
for cluster, color, label in (
    ('hyades', '#267783', r'Hyades: [Fe/H] = +0.15'),
    ('pleiades', '#b96535', r'Pleiades: [Fe/H] = 0')):
    path = W / f'{cluster}_mass_sensitivity_draws.npz'
    with np.load(path, allow_pickle=False) as data:
        mass = data['mass_Msun']
    nref = ((mass >= REFERENCE[0]) & (mass < REFERENCE[1])).sum(axis=1)
    assert np.all(nref > 0), 'Every realization needs a populated reference interval.'
    counts = np.array([np.histogram(m, BINS)[0] for m in mass])
    normalized = counts / nref[:, None] / WIDTH
    # The same exact reference boundaries are used for counts and histograms.
    area = (normalized[:, REF_BINS] * WIDTH[REF_BINS]).sum(axis=1)
    assert np.allclose(area, 1., rtol=0, atol=1e-12)
    low, median, high = np.quantile(normalized, [.16, .50, .84], axis=0)
    nquant = np.quantile(nref, [.16, .50, .84])
    reference_counts[cluster] = nref
    receipt['clusters'][cluster] = dict(
        normalization_count_quantiles=nquant.tolist(),
        maximum_normalization_error=float(np.max(abs(area - 1))),
        minimum_model_mass_Msun=source['clusters'][cluster]['minimum_mass_model'],
        input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        age_scenarios_years=source['age_scenarios_years'][cluster],
        relative_density_quantiles=[low.tolist(), median.tolist(), high.tolist()])
    # Do not draw unsupported bins or pretend a zero median is measurable on
    # a log axis. Positive upper sensitivity bounds remain visible as shading.
    minimum = source['clusters'][cluster]['minimum_mass_model']
    supported = BINS[:-1] >= minimum * (1 - 1e-10)
    valid = supported & (high > 0)
    lower = np.where(valid, np.maximum(low, 1e-5), np.nan)
    upper = np.where(valid, high, np.nan)
    ax.fill_between(BINS, np.r_[lower, lower[-1]], np.r_[upper, upper[-1]],
                    step='post', color=color, alpha=.18, linewidth=0, zorder=1)
    ax.stairs(np.where(supported & (median > 0), median, np.nan), BINS,
              baseline=None, color=color, linewidth=2.3, label=label, zorder=3)
    for i, (a, b) in enumerate(zip(BINS[:-1], BINS[1:])):
        rows.append(dict(cluster=cluster, mass_lower_Msun=a, mass_upper_Msun=b,
                         model_covered=bool(supported[i]),
                         relative_density_low=low[i], relative_density_median=median[i],
                         relative_density_high=high[i]))

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(.008, .25)
ax.set_ylim(.005, 12)
ax.xaxis.set_major_locator(FixedLocator([.01, .02, .04, .08, .10, .15, .20, .25]))
ax.xaxis.set_major_formatter(FuncFormatter(lambda value, pos: f'{value:g}'))
ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xlabel(r'Ember mass ($M_\odot$)', labelpad=8)
ax.set_ylabel('Relative number per dex in mass', labelpad=8)
ax.grid(axis='y', which='major', alpha=.13)
ax.legend(loc='upper left', frameon=False, fontsize=12,
          handlelength=2.6, borderaxespad=.7, labelspacing=.65)
ax.text(.035, .785, 'Includes identified tails and associated members',
        transform=ax.transAxes, fontsize=10.5, color='#444444')
ax.text(np.sqrt(.10 * .20), 10.5, 'Equal counts over\n' + r'0.10–0.20 $M_\odot$',
        ha='center', va='top', fontsize=10, color='#45484b')
fig.text(.11, .063, 'Shading: age, clouds, initial conditions and luminosity sensitivity.',
         fontsize=10, color='#333333')
fig.text(.11, .030, 'Observed candidates; incompleteness, contamination and unresolved binaries remain uncorrected.',
         fontsize=9.7, color='#333333')
fig.subplots_adjust(left=.115, right=.975, bottom=.195, top=.955)
stem = 'normalized_observed_mass_functions'
for extension in ('png', 'pdf'):
    output = W / f'{stem}.{extension}'
    fig.savefig(output)
    shutil.copy2(output, DOC / f'CLUSTER_MASS_FUNCTIONS_NORMALIZED_20261009.{extension}')
with (W / f'{stem}_bins.csv').open('w') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
(W / f'{stem}_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps({cl: item['normalization_count_quantiles']
                  for cl, item in receipt['clusters'].items()}, indent=2))
