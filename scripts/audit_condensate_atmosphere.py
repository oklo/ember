#!/usr/bin/env python3
"""Independent chemistry/convection audit of one completed depleted atmosphere.

Uses numpy and the separately built FastChem source probe. The settled-grain
experiment omits grain enthalpy; reject condensation in heat-carrying
convective layers instead of silently applying an incomplete adiabat there.
"""
import argparse
import gzip
import json
from pathlib import Path
import subprocess
import numpy as np
from audit_condensates import pack,sha
from generate_nongrey_grid import composition
from prepare_nongrey_sources import SOURCES
from condensate_chemistry_profile import audit_profile,HOT_JOIN_TEMPERATURE


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('model',type=Path)
    p.add_argument('--diagnostic-control',action='store_true',
                   help='audit an explicit heat-capacity control; report unsupported grain transport without accepting a production model')
    a=p.parse_args()
    validated=json.loads((a.model/('initializer.json' if a.diagnostic_control else 'validated.json')).read_text());provenance=validated['provenance']
    if a.diagnostic_control and provenance['prepared'].get('initialization_only',{}).get('method')!='outer-layer heat-capacity sensitivity control; not an accepted atmosphere':
        raise ValueError('explicit diagnostic thermal-control source required')
    if provenance['mode']!='equilibrium':raise ValueError('expected a depleted-gas atmosphere')
    receipt=json.loads((a.model/'completed.json').read_text())
    for name,digest in receipt['outputs'].items():
        if sha((a.model/name).read_bytes())!=digest:raise ValueError('completed source changed')
    fc=provenance['prepared']['depletion']['fastchem']
    if sha(Path(fc['executable']).read_bytes())!=fc['executable_sha256']:raise ValueError('chemistry probe changed')
    profile=[]
    for line in (a.model/'run.log').read_text().rsplit('FINAL MODEL ATMOSPHERE',1)[1].splitlines():
        t=line.replace('D','E').split()
        if len(t)==11 and t[0].isdigit():profile.append(list(map(float,t)))
    profile.reverse()
    inputs=''.join(f'{r[3]:.17g} {r[6]/1e6:.17g}\n' for r in profile if r[3]<=HOT_JOIN_TEMPERATURE)
    result=subprocess.run([fc['executable'],fc['source'],str((a.model/'condensate-abundances.dat').resolve()),'equilibrium'],
        input=inputs,text=True,capture_output=True)
    (a.model/'chemistry-audit-input.dat').write_text(inputs)
    (a.model/'chemistry-audit-source.json.gz').write_bytes(gzip.compress(result.stdout.encode(),mtime=0))
    (a.model/'chemistry-audit.log').write_text(result.stderr)
    result.check_returncode()
    chemistry=json.loads(result.stdout)
    hot=[r for r in profile if r[3]>HOT_JOIN_TEMPERATURE];hot_join=None
    if hot:
        join_inputs=''.join(f'{HOT_JOIN_TEMPERATURE:.17g} {r[6]/1e6:.17g}\n' for r in hot)
        join_result=subprocess.run([fc['executable'],fc['source'],str((a.model/'condensate-abundances.dat').resolve()),'equilibrium'],
                                  input=join_inputs,text=True,capture_output=True)
        (a.model/'chemistry-hot-join-input.dat').write_text(join_inputs)
        (a.model/'chemistry-hot-join-source.json.gz').write_bytes(gzip.compress(join_result.stdout.encode(),mtime=0))
        (a.model/'chemistry-hot-join.log').write_text(join_result.stderr)
        join_result.check_returncode();hot_join=json.loads(join_result.stdout)
    abundance,weights=composition(provenance['XH'],provenance['X3'],provenance['specification']['metals'])
    metadata=json.loads((SOURCES/'synple-elements.json').read_text());numbers={};masses={}
    for s,n,w in zip(metadata['symbol'],abundance,weights,strict=True):
        if n>1e-90:numbers[s.capitalize()]=n;masses[s.capitalize()]=w*1.67333e-24/1.66053906660e-24
    audit,condensed=audit_profile(chemistry,profile,numbers,masses,hot_join)
    # 1e-14 condensate particles per nucleus is the existing chemistry
    # audit's reporting floor. A 1e-8 heat-flux fraction excludes roundoff.
    flux=np.array([r[9] for r in profile])
    audit['maximum_convective_flux_fraction_in_condensing_layers']=float(np.abs(flux[condensed]).max()) if condensed.any() else 0.
    audit['grain_enthalpy_supported']=not bool(np.any(condensed & (np.abs(flux)>1e-8)))
    audit['provenance']={'atmosphere_receipt_sha256':sha((a.model/'completed.json').read_bytes()),
        'chemistry_source_sha256':sha((a.model/'chemistry-audit-source.json.gz').read_bytes()),'probe_sha256':fc['executable_sha256']}
    if hot:audit['provenance']['hot_join_source_sha256']=sha((a.model/'chemistry-hot-join-source.json.gz').read_bytes())
    if a.diagnostic_control:audit['scope']='Diagnostic heat-capacity perturbation; not an accepted atmosphere or a grain EOS.'
    (a.model/('chemistry-diagnostic.json' if a.diagnostic_control else 'chemistry-audit.json')).write_bytes(pack(audit))
    print(json.dumps(audit,indent=2))
    if not a.diagnostic_control and not audit['grain_enthalpy_supported']:raise ValueError('condensing convective layers require grain enthalpy physics')


if __name__=='__main__':main()
