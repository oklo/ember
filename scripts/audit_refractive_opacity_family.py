#!/usr/bin/env python3
"""Check assembled refractive opacities through Ember's production classes.

Test every composition/temperature row, independent source means and spectra,
material derivatives, isotope mapping, domain guards, and the fixed current
profile. Retain failures of the unblended cold radiative component separately
from the high-temperature combined transport comparison. No stellar evolution
or complete physical acceptance follows from this material-state audit alone.
"""
import argparse
import bisect
import csv
import gzip
import json
import math
from pathlib import Path
import subprocess

from audit_tops_electron_dispersion import KEV, KB
from fetch_tops_composition import digest
from reduce_tops_group_factors import add_inputs, verify


def run(probe, mode, family, points, scratch, label):
    request = ''.join(' '.join(format(x, '.17g') for x in q)+'\n' for q in points)
    command = [str(probe), mode, str(family), 'data/conduction/condtab21wd_metals.dat']
    response = subprocess.run(command, input=request, text=True, capture_output=True, timeout=60)
    raw = {'command': command, 'request': request, 'returncode': response.returncode,
           'stdout': response.stdout, 'stderr': response.stderr}
    archive = scratch / (label+'.json.gz')
    archive.write_bytes(gzip.compress(json.dumps(raw, separators=(',', ':')).encode(), mtime=0))
    response.check_returncode()
    rows = [json.loads(line) for line in response.stdout.splitlines()]
    if len(rows) != len(points) or any(row['query'] != list(q) for row, q in zip(rows, points, strict=True)):
        raise ValueError('runtime changed query order or coordinates')
    for row in rows:
        if row['covered'] and any(not math.isfinite(row[k]) for k in
                                   ('kappa', 'dlnk_dlnT', 'dlnk_dlnrho')):
            raise ValueError('nonfinite covered runtime value')
    return rows


def read_table(path):
    lines = path.read_text().splitlines()
    head = lines[0].split()
    if head[:2] != ['EMBER_OPACITY_TABLE', '2']:
        raise ValueError('expected density-prefix opacity format')
    nx, nt, nr = map(int, head[2:5])
    rr, tt = [list(map(float, row.split())) for row in lines[1:3]]
    if len(rr) != nr or len(tt) != nt:
        raise ValueError('invalid table axes')
    planes = {}
    offset = 3
    for _ in range(nx):
        x, z = map(float, lines[offset].split()); offset += 1
        values = []
        for _ in range(nt):
            row = lines[offset].split(); offset += 1
            n = int(row[0]); data = list(map(float, row[1:]))
            if n != len(data) or not 4 <= n <= nr:
                raise ValueError('invalid table density prefix')
            values.append(data)
        if x in planes:
            raise ValueError('duplicate table composition')
        planes[x] = values
    if offset != len(lines):
        raise ValueError('trailing table rows')
    return {'X': sorted(planes), 'Z': z, 'logT': tt, 'logrho': rr, 'planes': planes}


def stable_change(delta, rad, cond):
    return delta if cond is None else delta/(1+(1+delta)*rad/cond)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('assembly', 'family', 'old_family', 'probe', 'probe_build', 'mass_probe',
                 'profile', 'dense_spectra', 'composition_spectra', 'scratch', 'output'):
        p.add_argument('--'+name.replace('_', '-'), type=Path, required=True)
    p.add_argument('--reference', type=Path, nargs='+', required=True)
    a = p.parse_args()
    if a.scratch.exists() or a.output.exists():
        raise FileExistsError('preserve previous runtime comparisons')
    assembly, build, dense, mixed = [json.loads(path.read_text()) for path in
                                     (a.assembly, a.probe_build, a.dense_spectra, a.composition_spectra)]
    inputs = {}
    for report in (assembly, build, dense, mixed):
        add_inputs(inputs, report['input_sha256'])
        add_inputs(inputs, report.get('output_sha256', {}))
    if digest(a.probe) != build['executable_sha256'] or build['exit_code'] != 0:
        raise ValueError('runtime probe differs from completed build')
    references = []
    for path in a.reference:
        r = json.loads(path.read_text()); references.append(r)
        add_inputs(inputs, r['input_sha256']); add_inputs(inputs, r['output_sha256'])
        add_inputs(inputs, {str(path.resolve()): digest(path)})
    for path in (a.assembly, a.probe, a.probe_build, a.mass_probe, a.profile,
                 a.dense_spectra, a.composition_spectra, Path(__file__),
                 Path('/tmp/ember-opacity-mass-basis-v1.cpp'),
                 Path('data/conduction/condtab21wd_metals.dat')):
        add_inputs(inputs, {str(path.resolve()): digest(path)})
    for path in a.old_family.glob('*.dat'):
        add_inputs(inputs, {str(path.resolve()): digest(path)})
    verify(inputs)
    a.scratch.mkdir(parents=True)
    tables = [read_table(a.family/f'tops_gs98_mixture_z{round(z*1000):03d}_high.dat')
              for z in assembly['Z']]
    nodes, truth = [], []
    for table in tables:
        xs, tt, rr = table['X'], table['logT'], table['logrho']
        for x in xs:
            ix = min(max(0, bisect.bisect_right(xs, x)-1), len(xs)-2)
            for it, logt in enumerate(tt):
                # Exact knots can round toward either temperature interval.
                count = min(len(table['planes'][xx][j]) for xx in xs[ix:ix+2]
                            for j in range(max(0,it-2), min(len(tt),it+3)))
                for ir in sorted({0, count//2, count-1}):
                    nodes.append((x, table['Z'], 0., 10**logt, 10**rr[ir]))
                    truth.append(10**table['planes'][x][it][ir])
    node_values = run(a.probe, 'high', a.family, nodes, a.scratch, 'nodes')
    if any(not r['covered'] for r in node_values):
        raise ValueError('a conservatively supported table node was rejected')
    node_error = max(abs(r['kappa']/value-1) for r, value in zip(node_values, truth, strict=True))
    independent, points = [], []
    for ref in references:
        for row in ref['records']:
            points.append((ref['X'], ref['Z'], 0., row['temperature_keV']*KEV/KB, row['density_atomic_g_cm3']))
            independent.append({'X': ref['X'], 'Z': ref['Z'], **row})
    values = run(a.probe, 'high', a.family, points, a.scratch, 'independent-groups')
    compared, omitted = [], []
    for source, row in zip(independent, values, strict=True):
        if not row['covered']:
            omitted.append({'source': source['source'], 'query': row['query'], 'reason': row['error']})
            continue
        direct = source['native_uncut_times_ratio_atomic_cm2_g']
        delta = row['kappa']/direct-1
        transport = stable_change(delta, direct, row['conductive_opacity']) if row.get('conduction_covered') else None
        compared.append({'X': source['X'], 'Z': source['Z'],
                         'temperature_keV': source['temperature_keV'],
                         'density_atomic_g_cm3': source['density_atomic_g_cm3'],
                         'source': source['source'], 'radiative_reference': direct,
                         'radiative_interpolated': row['kappa'], 'radiative_relative_error': delta,
                         'combined_relative_error': transport,
                         'combined_absolute_error_bound': abs(delta) if transport is None else abs(transport)})
    spectra = [{'X': 0., 'Z': .02, 'temperature_keV': r['temperature_keV'],
                'density_atomic_g_cm3': r['density_atomic_g_cm3'], 'reference': r['reference_opacity']}
               for r in dense['spectral_comparisons']]
    spectra += [{'X': r['X'], 'Z': r['Z'], 'temperature_keV': r['temperature_keV'],
                 'density_atomic_g_cm3': r['density_atomic_g_cm3'],
                 'reference': r['reference_rosseland_atomic_cm2_g']}
                for r in mixed['records']]
    points = [(r['X'],r['Z'],0.,r['temperature_keV']*KEV/KB,r['density_atomic_g_cm3']) for r in spectra]
    values = run(a.probe,'high',a.family,points,a.scratch,'full-spectra')
    spectrum_checks, spectrum_omissions = [], []
    for ref, row in zip(spectra,values,strict=True):
        if row['covered']:
            delta = row['kappa']/ref['reference']-1
            spectrum_checks.append({**ref, 'actual':row['kappa'],'relative_error':delta})
        else:
            spectrum_omissions.append({**ref,'reason':row['error']})
    # Scalar derivatives are checked at interior composition/thermal points.
    centers = [(x,z,0.,T,rho) for x in (.0125,.0875,.1375,.1875,.35,.65,.72)
               for z in (.015,.025) for T,rho in ((4.2e5,.3),(4.8e5,3.),(2e6,100.),(1.3e7,5000.))]
    centers += [(x,.02,.05,1.3e7,12000.) for x in (.0125,.0875,.1375)]
    derivative_records=[]
    for mode, selected in [('high',[q for q in centers if q[2]==0]),
                           ('atomic-blend',[q for q in centers if q[2]==0 and q[3]<1e6]),
                           ('stellar',[q for q in centers if q[2]>0])]:
        queries=[]; axes=[]
        for point in selected:
            queries.append(point); local=[]
            for index,field,h,logarithmic in [(3,'dlnk_dlnT',2e-5,True),(4,'dlnk_dlnrho',2e-5,True),
                                            (0,'dlnk_dX',1e-6,False),(1,'dlnk_dZ',1e-6,False)]:
                if mode=='stellar' and index==1:
                    continue
                lo,hi=list(point),list(point)
                lo[index]=point[index]*math.exp(-h) if logarithmic else point[index]-h
                hi[index]=point[index]*math.exp(h) if logarithmic else point[index]+h
                local.append((field,h,len(queries),len(queries)+1));queries.extend([tuple(lo),tuple(hi)])
            if mode=='stellar':
                lo,hi=list(point),list(point);lo[2]-=1e-6;hi[2]+=1e-6
                local.append(('dlnk_dY3',1e-6,len(queries),len(queries)+1));queries.extend([tuple(lo),tuple(hi)])
            axes.append((point,len(queries)-1-2*len(local),local))
        result=run(a.probe,mode,a.family,queries,a.scratch,'derivatives-'+mode)
        if any(not row['covered'] for row in result):
            raise ValueError('an interior derivative query is unsupported')
        for point,i,local in axes:
            for field,h,j,k in local:
                finite=(math.log(result[k]['kappa'])-math.log(result[j]['kappa']))/(2*h)
                actual=result[i][field]
                if actual is None:raise ValueError('missing required material derivative')
                error=abs(finite-actual)/max(1.,abs(finite),abs(actual))
                derivative_records.append({'mode':mode,'query':point,'field':field,'finite_difference':finite,
                                           'runtime':actual,'scaled_error':error})
    hs,hes,zs=map(float,subprocess.check_output([str(a.mass_probe)],text=True).split())
    isotope_points=[(x,.02,y3,T,rho) for x in (.0125,.0875,.1375,.35,.65)
                    for y3 in (.001,.05) for T,rho in ((2e6,100.),(1.3e7,5000.))]
    atomic_points=[];scales=[]
    for x,z,y3,T,rho in isotope_points:
        scale=x*hs+y3*hes*4/3+(1-x-z-y3)*hes+z*zs
        scales.append(scale);atomic_points.append((x*hs/scale,z*zs/scale,0.,T,rho*scale))
    mapped=run(a.probe,'stellar',a.family,isotope_points,a.scratch,'mapped-isotopes')
    atomic=run(a.probe,'high',a.family,atomic_points,a.scratch,'atomic-isotopes')
    if any(not r['covered'] for r in mapped+atomic):raise ValueError('unsupported isotope comparison')
    isotope_error=max(abs(b['kappa']/(s*r['kappa'])-1) for b,r,s in zip(mapped,atomic,scales,strict=True))
    profile=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(a.profile.open())]
    points=[(r['X'],.02,r['Y3'],r['temperature_K'],r['density_g_cm3']) for r in profile]
    old=run(a.probe,'stellar',a.old_family,points,a.scratch,'profile-old')
    new=run(a.probe,'stellar',a.family,points,a.scratch,'profile-new')
    if any(not r['covered'] or not r.get('conduction_covered') for r in old+new):
        raise ValueError('current stellar profile has unsupported transport')
    profile_checks=[]
    for i,(point,before,after) in enumerate(zip(points,old,new,strict=True)):
        profile_checks.append({'zone':i,'query':point,
                               'radiative_relative_change':after['kappa']/before['kappa']-1,
                               'combined_relative_change':after['combined_opacity']/before['combined_opacity']-1,
                               'density_limit_before':before['density_range'][1],
                               'density_limit_after':after['density_range'][1]})
    unchanged=[i for i,q in enumerate(points) if q[3]<10**5.6]
    lower_exact=all(all(old[i][k]==new[i][k] for k in
                      ('kappa','dlnk_dlnT','dlnk_dlnrho','dlnk_dX','dlnk_dY3')) for i in unchanged)
    bad=[(.01,.02,0.,1.3e7,100100.),(.3,.02,0.,1.3e7,10100.),(.1,.02,0.,4e5,10100.)]
    rejected=run(a.probe,'high',a.family,bad,a.scratch,'unsupported')
    guards=all(not row['covered'] for row in rejected)
    high=[r for r in compared if r['temperature_keV']*KEV/KB>=10**5.7]
    if not high or not spectrum_checks:raise ValueError('missing independent comparisons')
    derivative_max=max(r['scaled_error'] for r in derivative_records)
    spectral_max=max(abs(r['relative_error']) for r in spectrum_checks)
    combined_max=max(r['combined_absolute_error_bound'] for r in high)
    numerical=node_error<1e-11 and derivative_max<3e-5 and isotope_error<1e-11 and lower_exact and guards
    verify(inputs)
    result={'scope':__doc__,'accepted_for_stellar_opacity':False,
            'source_row_node_queries':len(nodes),'maximum_node_relative_error':node_error,
            'independent_group_comparisons':len(compared),'independent_group_omissions':omitted,
            'maximum_unblended_radiative_relative_error':max(abs(r['radiative_relative_error']) for r in compared),
            'fully_selected_high_group_comparisons':len(high),'high_combined_absolute_error_upper_bound':combined_max,
            'independent_high_transport_comparison_passed':combined_max<=.005,
            'full_spectrum_comparisons':len(spectrum_checks),'full_spectrum_omissions':spectrum_omissions,
            'maximum_full_spectrum_relative_error':spectral_max,'full_spectrum_comparison_passed':spectral_max<=.005,
            'relative_comparison_criterion':.005,'derivative_checks':derivative_records,
            'maximum_scaled_derivative_error':derivative_max,'isotope_mapping_queries':len(isotope_points),
            'maximum_isotope_mapping_relative_error':isotope_error,'unsupported_queries_rejected':guards,
            'profile_zones':len(points),'profile_lower_temperature_zones':len(unchanged),
            'profile_lower_temperature_values_and_derivatives_exact':lower_exact,
            'maximum_profile_radiative_relative_change':max(abs(r['radiative_relative_change']) for r in profile_checks),
            'maximum_profile_combined_relative_change':max(abs(r['combined_relative_change']) for r in profile_checks),
            'profile_checks':profile_checks,'runtime_numerical_checks_passed':numerical,
            'group_comparisons':compared,'spectrum_checks':spectrum_checks,
            'limitations':['The lower-temperature opacity prescription is held fixed, not independently revalidated by its unchanged bytes.',
                           'The actual temperature blend has derivative and profile checks here, but no independent full-spectrum blend comparison.',
                           'The full-spectrum comparisons retain their limited sampled mixtures and native temperatures.',
                           'Global cold radiative errors and earlier strict group/native recovery failures remain explicit.',
                           'A changed-physics stellar calculation and its evolution checks remain separate.'],
            'input_sha256':inputs,'scratch_sha256':{path.name:digest(path) for path in a.scratch.iterdir() if path.is_file()}}
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in
                     ('scope','derivative_checks','profile_checks','group_comparisons','spectrum_checks','input_sha256',
                      'scratch_sha256','limitations','independent_group_omissions','full_spectrum_omissions')}),flush=True)
    if not numerical or not result['independent_high_transport_comparison_passed'] or not result['full_spectrum_comparison_passed']:
        raise ValueError('assembled family fails one or more recorded numerical or independent comparisons')


if __name__=='__main__':
    main()
