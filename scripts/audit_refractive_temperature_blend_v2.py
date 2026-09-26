#!/usr/bin/env python3
"""Compare the actual opacity transition to independent native-temperature sources.

Keep lower-temperature discrepancies distinct from changes introduced by the
new high-temperature family. No new source temperatures or spectra are implied.
"""
import argparse
import json
import math
from pathlib import Path

import audit_refractive_opacity_family_v2 as runtime_tools
from audit_refractive_opacity_family_v2 import run, stable_change
from audit_tops_electron_dispersion import KEV, KB
from reduce_tops_group_factors import add_inputs, verify
from fetch_tops_composition import digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('runtime', 'family', 'old_family', 'probe', 'scratch', 'output'):
        p.add_argument('--'+name.replace('_', '-'), type=Path, required=True)
    p.add_argument('--reference', type=Path, nargs='+', required=True)
    p.add_argument('--prior-scratch', type=Path, required=True)
    p.add_argument('--prior-controller', type=Path, required=True)
    a = p.parse_args()
    if a.scratch.exists() or a.output.exists():
        raise FileExistsError('preserve completed diagnostics')
    previous = json.loads(a.runtime.read_text())
    if not all(previous[k] for k in ('runtime_numerical_checks_passed',
                                    'independent_high_transport_comparison_passed',
                                    'full_spectrum_comparison_passed')):
        raise ValueError('complete hot-family comparison must pass first')
    inputs = dict(previous['input_sha256'])
    references = []
    for path in a.reference:
        ref = json.loads(path.read_text())
        add_inputs(inputs, ref['input_sha256'])
        add_inputs(inputs, ref['output_sha256'])
        add_inputs(inputs, {str(path.resolve()): digest(path)})
        for row in ref['records']:
            if row['temperature_keV'] <= .06:
                references.append({'X': ref['X'], 'Z': ref['Z'], **row})
    for path in (a.runtime, a.probe, Path(__file__),
                 Path(__file__).with_name('audit_refractive_opacity_family_v2.py')):
        add_inputs(inputs, {str(path.resolve()): digest(path)})
    prior = json.loads(a.prior_controller.read_text())
    if prior['status'] != 'released' or prior['jobs'][0]['exit_code'] != 1:
        raise ValueError('expected terminal diagnostic with unsupported test locations')
    add_inputs(inputs, prior['jobs'][0]['input_sha256'])
    add_inputs(inputs, {str(a.prior_controller.resolve()): digest(a.prior_controller)})
    runtime_tools.REUSE_DIRECTORY = a.prior_scratch
    runtime_tools.REUSE_HASHES = {str(p.resolve()): digest(p) for p in a.prior_scratch.glob('*.json.gz')}
    add_inputs(inputs, runtime_tools.REUSE_HASHES)
    verify(inputs)
    a.scratch.mkdir()
    points = [(r['X'], r['Z'], 0., r['temperature_keV']*KEV/KB,
               r['density_atomic_g_cm3']) for r in references]
    new = run(a.probe, 'atomic-blend', a.family, points, a.scratch, 'reference-new')
    old = run(a.probe, 'atomic-blend', a.old_family, points, a.scratch, 'reference-old')
    comparisons, omissions = [], []
    for ref, n, o in zip(references, new, old, strict=True):
        if not n['covered'] or not o['covered']:
            omissions.append({'query': n['query'], 'new': n, 'old': o})
            continue
        direct = ref['native_uncut_times_ratio_atomic_cm2_g']
        dnew, dold = n['kappa']/direct-1, o['kappa']/direct-1
        cn = stable_change(dnew, direct, n['conductive_opacity']) if n.get('conduction_covered') else dnew
        co = stable_change(dold, direct, o['conductive_opacity']) if o.get('conduction_covered') else dold
        lt = math.log10(n['query'][3])
        region = 'lower_temperature' if lt <= 5.6 else 'transition' if lt < 5.7 else 'high_temperature'
        comparisons.append({'query': n['query'], 'region': region, 'source': ref['source'],
                            'reference': direct, 'new_opacity': n['kappa'], 'old_opacity': o['kappa'],
                            'new_radiative_relative_error': dnew, 'old_radiative_relative_error': dold,
                            'new_combined_error_or_radiative_bound': cn,
                            'old_combined_error_or_radiative_bound': co})
    # Check both transition endpoints and its interior through actual production
    # classes. The step is in ln T, and finite differences cross each endpoint.
    centers = [(x, z, 0., 10**lt, rho) for x in (.1, .175, .3, .7)
               for z in (.015, .02, .025) for rho in (.003, .03, .3, 3., 30., 300., 3000.)
               for lt in (5.6, 5.625, 5.65, 5.675, 5.7)]
    h = 1e-6
    query = []
    for q in centers:
        query.extend([q, (*q[:3], q[3]*math.exp(-h), q[4]),
                      (*q[:3], q[3]*math.exp(h), q[4])])
    values = run(a.probe, 'atomic-blend', a.family, query, a.scratch, 'continuity')
    checks, unsupported_centers = [], []
    for i, q in enumerate(centers):
        c, lo, hi = values[3*i:3*i+3]
        if not all(r['covered'] for r in (c, lo, hi)):
            unsupported_centers.append({'query': q, 'responses': [c, lo, hi]})
            continue
        finite = (math.log(hi['kappa'])-math.log(lo['kappa']))/(2*h)
        error = abs(finite-c['dlnk_dlnT'])/max(1., abs(finite), abs(c['dlnk_dlnT']))
        checks.append({'query': q, 'scaled_derivative_error': error,
                       'left_value_relative_change': lo['kappa']/c['kappa']-1,
                       'right_value_relative_change': hi['kappa']/c['kappa']-1})
    summary = {}
    for region in ('lower_temperature', 'transition', 'high_temperature'):
        rows = [r for r in comparisons if r['region'] == region]
        if not rows:
            raise ValueError('missing independent comparison region')
        summary[region] = {'comparisons': len(rows),
                          'maximum_new_combined_error_or_bound': max(abs(r['new_combined_error_or_radiative_bound']) for r in rows),
                          'maximum_old_combined_error_or_bound': max(abs(r['old_combined_error_or_radiative_bound']) for r in rows),
                          'new_failures_above_criterion': sum(abs(r['new_combined_error_or_radiative_bound']) > .005 for r in rows),
                          'maximum_opacity_change': max(abs(r['new_opacity']/r['old_opacity']-1) for r in rows)}
    derivative_max = max(r['scaled_derivative_error'] for r in checks)
    passed = (derivative_max < 3e-5 and summary['transition']['new_failures_above_criterion'] == 0
              and summary['high_temperature']['new_failures_above_criterion'] == 0
              and summary['lower_temperature']['maximum_opacity_change'] == 0)
    verify(inputs)
    report = {'scope': __doc__, 'transition_checks_passed': passed,
              'accepted_for_stellar_opacity': False, 'relative_comparison_criterion': .005,
              'summary': summary, 'maximum_scaled_derivative_error': derivative_max,
              'continuity_centers': len(checks), 'comparisons': comparisons,
              'unsupported_continuity_centers': unsupported_centers,
              'reused_archives': runtime_tools.REUSED_ARCHIVES,
              'omissions': omissions, 'continuity_checks': checks,
              'limitations': ['The direct comparison uses native source temperatures, including 0.04 keV inside the transition.',
                              'Group transport is normalized to the native uncut mean; earlier group/native recovery failures remain recorded.',
                              'No full spectrum lies inside this transition in the current independent controls.',
                              'Existing lower-temperature errors are reported without changing their source data or conduction activation.',
                              'Finite differences check interpolation consistency, not independent thermal accuracy.'],
              'input_sha256': inputs,
              'scratch_sha256': {p.name: digest(p) for p in a.scratch.iterdir() if p.is_file()}}
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'passed': passed, 'summary': summary, 'derivative_error': derivative_max}), flush=True)
    if not passed:
        raise ValueError('one or more transition checks fail; retain all diagnostics')


if __name__ == '__main__':
    main()
