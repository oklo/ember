#!/usr/bin/env python3
"""Build a local OP comparison family without re-running atomic calculations.

Native element spectra are mixed on their electron-density grid. OP mixture
corrections and finite-frequency plasma transport are integrated locally.
Four missing trace elements retain mass but supply no charge or opacity.
This is a conditional source family, not a selected stellar input.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import re
import time

import numpy as np
from scipy.interpolate import PchipInterpolator

from op_native_spectra_v2 import read_spectra, read_mesh, native_rosseland
from op_mixture_corrections import corrections
from op_plasma_transport import finite_means
from audit_tops_electron_dispersion import electron_moments, HBAR, ME, KB

ROOT = Path(__file__).resolve().parents[1]
WORK = Path('/tmp/ember-op-warm-family-v1')
PREVIOUS = ROOT/'docs/results/op_warm_envelope_v1.json'
OUTPUT = ROOT/'docs/results/op_warm_family_v1.json'
MESH = Path('/tmp/ember-op-selected-spectra-v1/OP4STARS_1.3/mono/m01.mesh')
XGRID = [0.,.001,.01,.025,.05,.075,.1,.125,.15,.175,.2,.25,.3,.4,.5,.6,.7,.75]
ZGRID = [.01,.02,.03]
TEMPERATURES = list(range(172,208,2))
LOGR_MIN, LOGR_MAX = -3., 2.25
MU = 1.66053906660e-24
A0 = 5.29177210903e-9
QE = 4.80320471257e-10


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic_json(path, value):
    temp = path.with_suffix('.json.tmp')
    temp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    temp.replace(path)


@lru_cache(maxsize=1)
def configuration():
    return json.loads((WORK/'configuration.json').read_text())


@lru_cache(maxsize=1)
def native_elements(it):
    config = configuration()
    return {int(z):read_spectra(path) for z,path in config['source_paths'][str(it)].items()}


@lru_cache(maxsize=1)
def mesh():
    return read_mesh(MESH)


@lru_cache(maxsize=2048)
def electron_state(it, jn):
    t = 10**(.025*it); ne = 10**(.25*jn)
    e = electron_moments(t,ne)
    cutoff = HBAR*math.sqrt(4*math.pi*QE*QE*ne/ME)/(KB*t)*math.sqrt(e['plasma_frequency_squared_ratio'])
    return cutoff,e['vstar_squared']


def integrate_state(it, jn, sigma, ions, electrons):
    dv,u = mesh(); t = 10**(.025*it); ne = 10**(.25*jn)
    conversion = A0*A0/MU*(-np.expm1(-u))
    means = []
    for order in (32,64):
        c = corrections(t,ne,ions,electrons,u,order)
        total = sigma+c['scattering_change']+c['freefree_change']
        if np.any(total <= 0):
            raise ValueError('nonpositive corrected spectrum')
        means.append(float(A0*A0/MU*native_rosseland(total,dv)))
    thermal_change = means[-1]/means[-2]-1
    if abs(thermal_change) > 1e-6:
        raise ValueError('free-free thermal quadrature is unresolved')
    remainder = total-c['scattering']
    if np.any(remainder <= 0):
        raise ValueError('known scattering exceeds total opacity')
    cutoff,vstar = electron_state(it,jn)
    trials = []
    for order in (8,16,32):
        f = finite_means(u,remainder*conversion,c['scattering']*conversion,cutoff,vstar,order)
        trials.append(f)
        if len(trials)>1:
            change = max(abs(f[k]/trials[-2][k]-1) for k in ['uncut_finite_mean','scattering_extreme_mean','absorption_extreme_mean'])
            if change < 1e-8:
                break
    if change >= 1e-8:
        raise ValueError('frequency quadrature is unresolved')
    return dict(ordinary=f['uncut_finite_mean'],scattering=f['scattering_extreme_mean'],
                absorption=f['absorption_extreme_mean'],thermal_quadrature_change=thermal_change,
                frequency_quadrature_change=change,continuous_relative_to_native=f['uncut_finite_mean']/means[-1]-1,
                electron_eta=c['electron_response']['eta'],cutoff_u=cutoff)


def plane_task(task):
    it,x,z = task
    config = configuration(); elements = native_elements(it)
    weights = {1:x/1.00794,2:(1-x-z)/4.002602}
    weights.update({int(e[0]):z*e[3]/e[2] for e in config['metals'] if str(int(e[0])) in config['source_paths'][str(it)]})
    weights = {k:w for k,w in weights.items() if w>0}
    jns = np.array([s['electron_index'] for s in elements[1]['states']])
    n = len(jns); sigma = np.zeros((n,10000)); electrons = np.zeros(n)
    ions = np.zeros((n,28)); valid = np.ones(n,dtype=bool)
    for element,w in weights.items():
        states = elements[element]['states']
        if [s['electron_index'] for s in states] != jns.tolist():
            raise ValueError('unaligned native electron-density grids')
        for j,s in enumerate(states):
            sigma[j] += w*s['cross_section_atomic']
            electrons[j] += w*s['electrons_per_atom']
            valid[j] &= s['positive_cross_section']
            for i,fraction in zip(s['ion_indices'],s['ion_fractions']):
                q = element-1-int(i)
                if q>0:
                    ions[j,q-1] += w*fraction
    if np.any(electrons<=0):
        raise ValueError('nonpositive source electron count')
    rho = MU*10**(.25*jns)/electrons
    logr = np.log10(rho)-3*.025*it+18
    if np.any(np.diff(logr)<=0):
        raise ValueError('nonmonotonic native mixture density')
    # Use the connected original positive segment containing the requested
    # dilute edge. Never bridge a masked native spectrum.
    segments = np.split(np.arange(n),np.flatnonzero(np.diff(valid.astype(int)))+1)
    eligible = [s for s in segments if valid[s[0]] and len(s)>=4 and logr[s[0]]<=LOGR_MIN<=logr[s[-1]]]
    result = dict(temperature_index=it,logT=.025*it,X=x,Z=z,atom_weights=weights,
                  native_states=[],failures=[],source_masks=[int(jns[j]) for j in np.flatnonzero(~valid)])
    if len(eligible)!=1:
        result['status']='no_connected_dilute_support';return result
    segment = eligible[0]
    first = max(int(segment[0]),int(np.searchsorted(logr,LOGR_MIN))-2)
    last = min(int(segment[-1]),int(np.searchsorted(logr,LOGR_MAX))+1)
    for j in range(first,last+1):
        try:
            record = integrate_state(it,int(jns[j]),sigma[j],ions[j],float(electrons[j]))
        except (ValueError,RuntimeError) as error:
            result['failures'].append(dict(electron_index=int(jns[j]),logR=float(logr[j]),reason=str(error)))
            # Keep all outcomes; the assembler later chooses a connected
            # positive segment and never silently joins across a failure.
            continue
        result['native_states'].append(dict(electron_index=int(jns[j]),rho=float(rho[j]),logR=float(logr[j]),**record))
    result['status']='evaluated'
    return result


def connected_states(row):
    all_states = row['native_states']
    if len(all_states)<4:
        raise ValueError('fewer than four corrected native density states')
    # Source electron indices differ by two on the archived OP grid.
    segments = np.split(np.arange(len(all_states)),np.flatnonzero(np.diff([s['electron_index'] for s in all_states])!=2)+1)
    segments = [s for s in segments if len(s)>=4 and all_states[s[0]]['logR']<=LOGR_MIN<=all_states[s[-1]]['logR']]
    if len(segments)!=1:
        raise ValueError('no connected corrected segment at the dilute boundary')
    return [all_states[i] for i in segments[0]]


def write_tables(rows, step, label):
    folder = WORK/label;folder.mkdir()
    grid = np.arange(round((LOGR_MAX-LOGR_MIN)/step)+1)*step+LOGR_MIN
    lookup = {(r['temperature_index'],r['X'],r['Z']):r for r in rows}
    files = []
    for component in ['ordinary','absorption','scattering']:
        for z in ZGRID:
            lines = [f'EMBER_OPACITY_TABLE 2 {len(XGRID)} {len(TEMPERATURES)} {len(grid)} OP_GS98_17_elements_finite_frequency_{component}',
                     ' '.join(f'{v:.17g}' for v in grid),
                     ' '.join(f'{.025*it:.17g}' for it in TEMPERATURES)]
            for x in XGRID:
                lines.append(f'{x:.17g} {z:.17g}')
                for it in TEMPERATURES:
                    states = connected_states(lookup[it,x,z])
                    rr = np.array([s['logR'] for s in states])
                    count = int(np.searchsorted(grid,rr[-1],side='right'))
                    if count<4:
                        raise ValueError('insufficient dense coverage for '+str((it,x,z)))
                    values = PchipInterpolator(rr,np.log10([s[component] for s in states]),extrapolate=False)(grid[:count])
                    if not np.isfinite(values).all():
                        raise ValueError('density regridding leaves native support')
                    lines.append(str(count)+' '+' '.join(f'{v:.17g}' for v in values))
            path = folder/f'op_gs98_z{round(z*1000):03d}_{component}.dat'
            with path.open('x') as f:
                f.write('\n'.join(lines)+'\n')
            files.append(dict(path=str(path),sha256=digest(path),bytes=path.stat().st_size,Z=z,component=component))
        path=folder/f'op_gs98_{component}.dat'
        lines=[f'EMBER_OPACITY_MIXTURE 1 3 logR OP_GS98_17_elements_finite_frequency_{component}']
        lines += [f'{z:.17g} "op_gs98_z{round(z*1000):03d}_{component}.dat"' for z in ZGRID]
        path.write_text('\n'.join(lines)+'\n')
        files.append(dict(path=str(path),sha256=digest(path),bytes=path.stat().st_size,component=component))
    return files


def main():
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    WORK.mkdir(exist_ok=True); (WORK/'planes').mkdir(exist_ok=True)
    previous=json.loads(PREVIOUS.read_text())
    inputs={str(p):digest(p) for p in [Path(__file__),PREVIOUS,MESH,ROOT/'include/ember/gs98_mixture.hpp',
        ROOT/'scripts/op_native_spectra_v2.py',ROOT/'scripts/op_mixture_corrections.py',
        ROOT/'scripts/op_plasma_transport.py',ROOT/'scripts/audit_tops_electron_dispersion.py']}
    paths={str(it):{} for it in TEMPERATURES}
    for name,h in previous['input_sha256'].items():
        match=re.fullmatch(r'm(\d+)\.(\d+)',Path(name).name)
        if match and int(match[2]) in TEMPERATURES:
            assert digest(name)==h
            paths[str(int(match[2]))][str(int(match[1]))]=name;inputs[name]=h
    assert all(len(p)==17 for p in paths.values())
    metals=[tuple(map(float,r.split(','))) for r in re.findall(r'^\s*\{([^{}]+)\}, //',(ROOT/'include/ember/gs98_mixture.hpp').read_text(),re.M)]
    config=dict(source_paths=paths,metals=metals,X=XGRID,Z=ZGRID,temperatures=TEMPERATURES,
                input_sha256=inputs,logR_min=LOGR_MIN,logR_max=LOGR_MAX)
    config_path=WORK/'configuration.json'
    if config_path.exists():
        assert json.loads(config_path.read_text())==config, 'resumption inputs changed'
    else:
        atomic_json(config_path,config)
    tasks=[(it,x,z) for it in TEMPERATURES for z in ZGRID for x in XGRID]
    rows=[];pending=[];reused=0
    def output_path(task):
        it,x,z=task;return WORK/'planes'/f't{it:03d}_x{round(x*1000000):06d}_z{round(z*1000000):06d}.json'
    for task in tasks:
        p=output_path(task)
        if p.exists():
            row=json.loads(p.read_text());assert (row['temperature_index'],row['X'],row['Z'])==task
            rows.append(row);reused+=1
        else:
            pending.append(task)
    start=time.monotonic();done=0
    print(f'OP warm family: {len(tasks)} planes, {len(pending)} new, {reused} reused',flush=True)
    with ProcessPoolExecutor(max_workers=8) as pool:
        futures={pool.submit(plane_task,task):task for task in pending}
        for future in as_completed(futures):
            task=futures[future];row=future.result()
            atomic_json(output_path(task),row);rows.append(row);done+=1
            if done%18==0 or done==len(pending):
                status=dict(completed=len(rows),total=len(tasks),elapsed_seconds=time.monotonic()-start,
                            native_states=sum(len(r['native_states']) for r in rows),failed_native_states=sum(len(r['failures']) for r in rows))
                atomic_json(WORK/'progress.json',status)
                print(json.dumps(status),flush=True)
    rows.sort(key=lambda r:(r['temperature_index'],r['Z'],r['X']))
    unsupported=[]
    for row in rows:
        try: connected_states(row)
        except ValueError as error: unsupported.append(dict(T=row['temperature_index'],X=row['X'],Z=row['Z'],reason=str(error)))
    files=[]
    if not unsupported:
        files=write_tables(rows,.05,'density005')+write_tables(rows,.025,'density0025')
    for p,h in inputs.items():assert digest(p)==h,p
    report=dict(outcome='assembled_conditional_warm_family' if not unsupported else 'incomplete_corrected_source_support',
        accepted_for_stellar_opacity=False,input_sha256=inputs,configuration=str(config_path),
        configuration_sha256=digest(config_path),planes=len(rows),reused_planes=reused,
        native_states=sum(len(r['native_states']) for r in rows),
        failed_native_states=sum(len(r['failures']) for r in rows),unsupported_planes=unsupported,
        plane_files=[dict(path=str(output_path((r['temperature_index'],r['X'],r['Z']))),sha256=digest(output_path((r['temperature_index'],r['X'],r['Z'])))) for r in rows],
        table_files=files,
        maximum_thermal_quadrature_change=max(abs(s['thermal_quadrature_change']) for r in rows for s in r['native_states']),
        maximum_frequency_quadrature_change=max(s['frequency_quadrature_change'] for r in rows for s in r['native_states']),
        limitations=['17 OP elements; P, Cl, K and Ti supply neither opacity nor charge while retaining reference mass.',
            'The three tables are ordinary transport and fixed absorption/scattering endpoint assumptions over the finite supplied frequency interval.',
            'Interpolating endpoint models does not establish bounds under nonlinear interpolation; no such bound is claimed for these tables.',
            'Native-density PCHIP followed by Ember temperature/density interpolation must be checked for regridding and composition errors.',
            'No atmosphere, full stellar step, or selected opacity input is changed.'])
    atomic_json(OUTPUT,report)
    print(json.dumps(dict(outcome=report['outcome'],report=str(OUTPUT),sha256=digest(OUTPUT),native_states=report['native_states'],unsupported_planes=len(unsupported))),flush=True)


if __name__=='__main__':
    main()
