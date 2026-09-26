#!/usr/bin/env python3
"""Compare archived MESA OPLIB means with saved, uncut TOPS means.

This is an offline source comparison, not a production table importer.
The MESA mixture has 25 elements; Ember's saved TOPS mixture has 21.
Two- and four-node tensor interpolation measure sensitivity, not accuracy.
Only source states inside the printed MESA grid are evaluated.
"""
import bisect
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import re

from import_tops_composition import read as read_tops

ROOT = Path(__file__).resolve().parents[1]
SUBSET = Path('/tmp/ember-mesa-oplib-subset-v1')
TASKS = Path('/tmp/ember-tops-cool-family-run-v1/tasks.json')
OUTPUT = ROOT / 'docs/results/mesa_oplib_mean_comparison_v2.json'
KELVIN_PER_KEV = 1e3 * 1.602176634e-12 / 1.380649e-16


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_table(entry):
    path = Path(entry['path'])
    assert digest(path) == entry['sha256']
    lines = path.read_text().splitlines()
    header = lines[0]
    assert '25 Elements' in header and 'Grevesse & Sauval (1998)' in header
    match = re.search(r'logT_max\s+(.+?)logT\s+logR', header)
    assert match is not None
    meta = list(map(float, match.group(1).split()))
    assert len(meta) == 10 and meta[:2] == [1, 37]
    _, _, x, z, nr, rmin, rmax, nt, tmin, tmax = meta
    assert (nr, rmin, rmax, nt, tmin, tmax) == (39, -8, 1.5, 213, 3.75, 9.05)
    name = re.fullmatch(r'oplib_gs98_z(.+)_x(.+)\.data', path.name)
    assert name is not None and (z, x) == tuple(map(float, name.groups()))
    rr = list(map(float, header.split('logRho - 3*logT + 18')[1].split()))
    assert rr == [-8 + .25 * i for i in range(39)]
    rows = [list(map(float, line.split())) for line in lines[1:] if line.strip()]
    assert len(rows) == 213 and all(len(row) == 40 for row in rows)
    assert all(math.isfinite(v) for row in rows for v in row)
    tt = [row[0] for row in rows]
    assert all(abs(v - (3.75 + .025 * i)) < 2e-14 for i, v in enumerate(tt))
    values = [row[1:] for row in rows]
    # Reject values outside this inspected subset's range. This is not a
    # general MESA missing-value decoder or proof of native source support.
    assert all(-8 < v < 8 for row in values for v in row)
    return dict(X=x, Z=z, logT=tt, logR=rr, logk=values, source=entry,
                logk_min=min(map(min, values)), logk_max=max(map(max, values)))


def stencil(grid, q, order):
    if not math.isfinite(q) or q < grid[0] or q > grid[-1]:
        raise ValueError('outside printed table domain')
    if order not in (2, 4) or len(grid) < order:
        raise ValueError('invalid interpolation order')
    lower = min(len(grid) - 2, max(0, bisect.bisect_right(grid, q) - 1))
    first = min(len(grid) - order, max(0, lower - (order // 2 - 1)))
    ii = list(range(first, first + order))
    ww = [math.prod((q - grid[j]) / (grid[i] - grid[j])
                    for j in ii if j != i) for i in ii]
    return ii, ww


def interpolate(table, t, r, order):
    ti, tw = stencil(table['logT'], t, order)
    ri, rw = stencil(table['logR'], r, order)
    return math.fsum(a * b * table['logk'][i][j]
                     for i, a in zip(ti, tw) for j, b in zip(ri, rw))


def controls():
    tt = [0, .3, .8, 1.4, 2.1]
    rr = [-2, -.8, -.3, .5, 1.1]
    polynomial = lambda t, r: 1 + .3*t - .2*r + .4*t*r + .03*t**3*r**2
    table = dict(logT=tt, logR=rr, logk=[[polynomial(t, r) for r in rr] for t in tt])
    checks = []
    for t, r in [(0, -2), (.12, -.41), (.75, .02), (1.81, .91), (2.1, 1.1)]:
        error = abs(interpolate(table, t, r, 4) - polynomial(t, r))
        assert error < 1e-13
        checks.append(dict(kind='independent_polynomial', logT=t, logR=r, absolute_error=error))
    bilinear = lambda t, r: 2 - .2*t + .7*r + .1*t*r
    table['logk'] = [[bilinear(t, r) for r in rr] for t in tt]
    for order in (2, 4):
        error = abs(interpolate(table, .66, -.53, order) - bilinear(.66, -.53))
        assert error < 1e-13
        checks.append(dict(kind='bilinear', order=order, absolute_error=error))
    for t, r in [(-.001, 0), (2.101, 0), (.5, -2.001), (.5, 1.101)]:
        try:
            interpolate(table, t, r, 4)
        except ValueError:
            checks.append(dict(kind='outside_domain_rejected', logT=t, logR=r))
        else:
            raise AssertionError('extrapolation accepted')
    return checks


def comparison(arg):
    table, tasks = arg
    x, z = table['X'], table['Z']
    rows, sources, skipped, exclusions = [], [], [], []
    plans = {}
    for task in tasks:
        planpath = ROOT / task['plan']
        if task['plan'] not in plans:
            plans[task['plan']] = json.loads(planpath.read_text())
        plan = plans[task['plan']]
        request_plan = plan['requests'][task['request_index']]
        t = request_plan['temperatures_keV'][0]
        assert len(request_plan['temperatures_keV']) == 1
        logt = math.log10(t * KELVIN_PER_KEV)
        densities = request_plan['densities_g_cm3']
        minr = math.log10(min(densities)) - 3 * logt + 18
        maxr = math.log10(max(densities)) - 3 * logt + 18
        if logt < table['logT'][0] or logt > table['logT'][-1] or minr > 1.5 or maxr < -8:
            skipped.append(dict(work=task['work'], reason='entire request outside printed MESA grid'))
            continue
        work = Path(task['work'])
        receipt_path = work / 'receipt.json'
        receipt = json.loads(receipt_path.read_text())
        assert (receipt['X'], receipt['Z']) == (x, z)
        assert receipt['plasma_cutoff'] == 'off'
        assert receipt['validation_role'] == 'official_form_uncut_groups'
        request_path = work / receipt['request']
        assert digest(request_path) == receipt['request_sha256']
        request = json.loads(request_path.read_text())
        assert request['plasnu'] == 'off' and request['lib'] == 'new'
        source_path = work / receipt['file']
        ts, rs, cells, excluded = read_tops(source_path, receipt, dimensions=tuple(receipt['dimensions']))
        assert len(ts) == 1 and abs(ts[0] / t - 1) < 1e-5
        source_index = len(sources)
        sources.append(dict(work=str(work), source_sha256=receipt['sha256'],
                            receipt_sha256=digest(receipt_path), request_sha256=receipt['request_sha256'],
                            plan=task['plan'], plan_sha256=digest(planpath),
                            request_index=task['request_index'], excluded_substituted_count=len(excluded),
                            source_bytes=source_path.stat().st_size))
        for (t, rho), kappa in sorted(cells.items()):
            if (t, rho) in excluded:
                exclusions.append(dict(source_index=source_index, T_keV=t, rho=rho, reason='TOPS substituted density'))
                continue
            logt = math.log10(t * KELVIN_PER_KEV)
            logr = math.log10(rho) - 3 * logt + 18
            if not (-8 <= logr <= 1.5 and table['logT'][0] <= logt <= table['logT'][-1]):
                exclusions.append(dict(source_index=source_index, T_keV=t, rho=rho, reason='outside printed MESA grid'))
                continue
            linear = 10 ** interpolate(table, logt, logr, 2)
            cubic = 10 ** interpolate(table, logt, logr, 4)
            rows.append(dict(source_index=source_index, T_keV=t, rho=rho, logT=logt, logR=logr,
                             tops_uncut=kappa, mesa_two_nodes=linear, mesa_four_nodes=cubic,
                             linear_relative_to_tops=linear/kappa-1, cubic_relative_to_tops=cubic/kappa-1,
                             interpolation_relative_change=cubic/linear-1))
    assert rows and len({(row['T_keV'], row['rho']) for row in rows}) == len(rows)
    def stats(key):
        vv = sorted(abs(row[key]) for row in rows)
        return dict(max_absolute=vv[-1], median_absolute=vv[len(vv)//2],
                    worst= max(rows, key=lambda row:abs(row[key])))
    return dict(X=x, Z=z, row_count=len(rows), sources=sources, skipped_requests=skipped,
                excluded_cells=exclusions, rows=rows,
                summaries={key:stats(key) for key in ['linear_relative_to_tops','cubic_relative_to_tops',
                                                      'interpolation_relative_change']})


def main():
    assert not OUTPUT.exists()
    analytic = controls()
    manifest = json.loads((SUBSET / 'manifest.json').read_text())
    tables = [read_table(entry) for entry in manifest['files']]
    tasks = json.loads(TASKS.read_text())
    args = []
    for table in tables:
        if table['Z'] not in [.01, .03] or table['X'] not in [0, .05]:
            continue
        component = f"/x{round(table['X']*1000):03d}-z{round(table['Z']*1000):03d}/"
        selected = [task for task in tasks if task['kind'] == 'groups' and component in task['work']]
        assert len(selected) == 76
        args.append((table, selected))
    assert len(args) == 4
    with ProcessPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(comparison, args))
    table_info = [{k:v for k,v in table.items() if k != 'logk'} for table in tables]
    report = dict(outcome='completed_source_comparison', accepted_for_stellar_opacity=False,
                  input_sha256={str(p):digest(p) for p in [Path(__file__),ROOT/'scripts/import_tops_composition.py',
                                 SUBSET/'manifest.json',TASKS]},
                  method='Log-opacity tensor Lagrange interpolation in logT and logR, two versus four nodes; no composition interpolation or domain extrapolation.',
                  analytic_controls=analytic, tables=table_info, comparisons=results,
                  assumptions=['MESA has 25 GS98 elements; saved TOPS requests have 21; these are not identical mixtures.',
                               'F, Sc, V and Co are present in MESA and absent from the TOPS mixture.',
                               'Compare the printed atomic mass density in both archives; no stellar baryonic-density conversion is attempted.',
                               'MESA tables are already interpolated from 74 to 213 temperatures; printed-grid support does not prove native ATOMIC support.',
                               'No explicit missing-value flag occurs in the extracted numeric subset; no general MESA missing-value convention is assumed.',
                               'TOPS means explicitly disable the plasma cutoff. MESA plasma conventions have not been independently established.',
                               'Differences combine abundance, source interpolation and physical conventions; they are not errors of the new interpolator.',
                               'No mean has been multiplied by an OP plasma factor or selected for stellar evolution.'],
                  references=['https://arxiv.org/html/2406.02845v1','https://docs.mesastar.org/en/latest/kap/overview.html'])
    with OUTPUT.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(report=str(OUTPUT),sha256=digest(OUTPUT),rows=sum(r['row_count'] for r in results),
                          comparisons=[{k:r[k] for k in ['X','Z','row_count','summaries']} for r in results]),indent=2))


if __name__ == '__main__':
    main()
