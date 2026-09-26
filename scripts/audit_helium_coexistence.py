#!/usr/bin/env python3
"""Solve pure-He liquid/solid coexistence at equal total pressure and Gibbs energy.

This is an offline source comparison, not an installed remnant EOS. It combines
EOSFI22 ideal-ion/interacting terms with Ember's ideal relativistic electrons.
No radiation, composition separation, partially ionized join or cooling age is
included. Both phase densities are solved; Gamma=175 only sets an initial bracket.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

RGAS=1.380649e-16/(4*1.66053906660e-24)
FIELDS=['F','U','S','P','cv','dP_dlnT','dP_dlnrho']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Probe:
    def __init__(self, executable, columns, work, name):
        self.errors=(work/(name+'.stderr')).open('w')
        self.raw=(work/(name+'.jsonl')).open('w')
        self.process=subprocess.Popen([str(Path(executable).resolve())],stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,stderr=self.errors,text=True,bufsize=1)
        self.columns=columns

    def query(self, values):
        self.process.stdin.write(' '.join(format(x,'.17g') for x in values)+'\n')
        self.process.stdin.flush()
        line=self.process.stdout.readline()
        result=list(map(float,line.split()))
        if len(result)!=len(self.columns) or not all(math.isfinite(v) for v in result):
            raise ValueError('invalid or incomplete source response: '+repr(line))
        self.raw.write(json.dumps({'input':values,'output':result},allow_nan=False)+'\n')
        return dict(zip(self.columns,result,strict=True))

    def close(self):
        self.process.stdin.close()
        code=self.process.wait(timeout=10)
        self.errors.close();self.raw.close()
        if code:raise RuntimeError('source probe failed')


class Matter:
    def __init__(self, ions, electrons):
        self.ions=ions;self.electrons=electrons
        self.cache={};self.electron_cache={}

    def state(self, phase, rho, T):
        key=(phase,rho,T)
        if key in self.cache:return self.cache[key]
        ion=self.ions.query([phase,2,4,rho,T])
        if (rho,T) not in self.electron_cache:
            self.electron_cache[rho,T]=self.electrons.query([rho,T])
        electron=self.electron_cache[rho,T]
        state={k:ion[k]+electron[k] for k in FIELDS}
        state.update({k:ion[k] for k in ['electron_rs','Gamma','Tp_over_T']})
        state.update(rho=rho,T=T,phase=phase,ion_rs=3*(ion['Gamma']/ion['Tp_over_T'])**2)
        state['G']=state['F']+state['P']/rho
        if state['P']<=0 or state['dP_dlnrho']<=0:
            raise ValueError('nonpositive total pressure or compressibility')
        self.cache[key]=state
        return state

    def pressure_pair(self, rho, T):
        liquid=self.state(0,rho,T);solid_rho=rho
        for iteration in range(30):
            solid=self.state(1,solid_rho,T)
            residual=math.log(solid['P']/liquid['P'])
            if abs(residual)<1e-13:
                return liquid,solid,(solid['G']-liquid['G'])/(RGAS*T)
            step=residual*solid['P']/solid['dP_dlnrho']
            solid_rho*=math.exp(-max(-.05,min(.05,step)))
        raise ValueError('equal-pressure density solve failed')

    def melting(self, rho):
        initial=self.state(0,rho,1e5)
        guess=initial['Gamma']*1e5/175
        lo,hi=map(math.log,[.6*guess,1.4*guess])
        if not self.pressure_pair(rho,math.exp(lo))[2]<0<self.pressure_pair(rho,math.exp(hi))[2]:
            raise ValueError('melting root is not bracketed')
        for iteration in range(70):
            mid=(lo+hi)/2;T=math.exp(mid)
            liquid,solid,error=self.pressure_pair(rho,T)
            if abs(error)<1e-10 and hi-lo<1e-8:break
            if error>0:hi=mid
            else:lo=mid
        else:raise ValueError('equal-Gibbs-energy solve failed')
        if not all(500<=s['ion_rs']<=120000 and s['Tp_over_T']<=30 and s['cv']>0
                   for s in [liquid,solid]):
            raise ValueError('coexistence outside the declared test range or nonpositive heat capacity')
        latent=T*(liquid['S']-solid['S'])
        dh=liquid['U']+liquid['P']/rho-solid['U']-solid['P']/solid['rho']
        dv=1/rho-1/solid['rho']
        if latent<=0 or dv<=0:raise ValueError('unexpected latent-heat or density-jump sign')
        return {'liquid':liquid,'solid':solid,'T_K':T,'pressure':liquid['P'],
            'relative_pressure_residual':solid['P']/liquid['P']-1,
            'Gibbs_difference_in_ideal_ion_units':error,
            'fractional_density_jump':solid['rho']/rho-1,
            'latent_heat_erg_g':latent,'latent_heat_in_ion_kT':latent/(RGAS*T),
            'enthalpy_identity_residual_in_ion_kT':(dh-latent)/(RGAS*T),
            'clapeyron_dP_dT':(liquid['S']-solid['S'])/dv}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['source_manifest','electron_probe','work','output']:
        parser.add_argument(name,type=Path)
    args=parser.parse_args()
    manifest=json.loads(args.source_manifest.read_text())
    if manifest['probe_kind']!='forced-phase interaction terms':raise ValueError('wrong source wrapper')
    if sha(manifest['executable'])!=manifest['executable_sha256']:raise ValueError('source executable changed')
    if args.work.exists():raise FileExistsError('use a new work directory')
    args.work.mkdir(parents=True)
    inputs={str(args.source_manifest):sha(args.source_manifest),str(args.electron_probe):sha(args.electron_probe),
            manifest['executable']:sha(manifest['executable'])}
    ions=Probe(manifest['executable'],FIELDS+['electron_rs','Gamma','Tp_over_T'],args.work,'ions')
    electrons=Probe(args.electron_probe,FIELDS,args.work,'electrons')
    matter=Matter(ions,electrons);rows=[]
    try:
        for rho in [1e3,3e3,1e4,3e4,1e5,3e5,1e6]:
            row=matter.melting(rho);checks=[]
            for h in [1e-3,3e-4]:
                before=matter.melting(rho*math.exp(-h));after=matter.melting(rho*math.exp(h))
                slope=(after['pressure']-before['pressure'])/(after['T_K']-before['T_K'])
                checks.append({'density_log_step':h,'finite_difference_dP_dT':slope,
                               'relative_difference':slope/row['clapeyron_dP_dT']-1})
            row['clapeyron_checks']=checks;rows.append(row)
            print(f"rho={rho:.4g}: Tm={row['T_K']:.4g} K, Gamma={row['liquid']['Gamma']:.4g}, density jump={row['fractional_density_jump']:.4g}",flush=True)
    finally:
        ions.close();electrons.close()
    if any(sha(p)!=value for p,value in inputs.items()):raise ValueError('input changed during calculation')
    report={'scope':__doc__,'source_manifest':manifest,'input_sha256':inputs,
        'script_sha256':sha(__file__),'electron_wrapper_sha256':sha(Path(__file__).with_name('helium_electron_probe.cpp')),
        'domain':'Pure fully ionized He4, fixed baryonic mass convention, 500<=ion_rs<=120000 and Tp/T<=30 at coexistence. Native interaction fits add physics beyond the rigid-background ion-only calculation in Baiko and Chugunov2022; this is not a reproduction of their ion-only melting fit.',
        'reference':'https://arxiv.org/abs/2112.04822','states':rows,
        'max_abs_clapeyron_relative_difference':max(abs(row['clapeyron_checks'][-1]['relative_difference']) for row in rows),
        'raw_sha256':{p.name:sha(p) for p in args.work.iterdir()}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
