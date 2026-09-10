#!/usr/bin/env python3
"""Recheck original depleted-atmosphere outputs, inputs and chemistry.

This offline validator needs numpy for independent species conservation.
It never accepts a numerical initializer as a completed atmosphere.
"""
import json
import hashlib
import math
from pathlib import Path
import numpy as np
from generate_nongrey_grid import (composition, input_fingerprint, temperatures, sequence,
                                   continuation_structure,resample_initial_structure)
from import_nongrey_grid import source_inputs, source_state, read_text
from prepare_nongrey_sources import digest, SOURCES
from condensate_chemistry_profile import audit_profile,HOT_JOIN_TEMPERATURE

PHYSICAL_KEYS = ['metals','log_temperature','temperature_K','log_density',
    'opacity_frequencies','opacity_method','wavelength_A','line_threshold',
    'synthesis_spacing_A','microturbulence_km_s','atmosphere_frequencies',
    'depths','alpha','tau_top','tau_bottom','tau','convection_derivative_step']


def original(directory, name):
    path=Path(directory)/name
    return path if path.exists() else path.with_name(name+'.gz')


def validate(directory, specification=None, prepared=None, coordinates=None, opacity_sha256=None,
             *, diagnostic_grain_transport=False):
    """Recheck a model; diagnostic mode reports an omitted-enthalpy failure.

    The diagnostic option is for controlled sensitivity experiments only.
    Production import, collection, and interpolation retain the default guard.
    All source, conservation and interior matching checks remain mandatory.
    """
    directory=Path(directory)
    provenance=json.loads(read_text(original(directory,'provenance.json')))
    spec=provenance['specification'];source=provenance['prepared']
    if provenance['mode']!='equilibrium' or 'depletion' not in source or 'initialization_only' in source:
        raise ValueError('canonical equilibrium-depleted source required')
    if specification and any(spec.get(k)!=specification.get(k) for k in PHYSICAL_KEYS):
        raise ValueError('condensate family physical inputs differ')
    if prepared and source!=prepared:raise ValueError('condensate family source identity differs')
    key=tuple(provenance[k] for k in ['XH','X3','teff_K','log_g'])
    if coordinates and key!=tuple(coordinates):raise ValueError('condensate coordinate mismatch')
    if opacity_sha256 and provenance['opacity_sha256']!=opacity_sha256:
        raise ValueError('condensate model used a different opacity plane')
    receipt=json.loads(read_text(original(directory,'completed.json')))
    for name,expected in receipt['outputs'].items():
        if hashlib.sha256(read_text(original(directory,name)).encode()).hexdigest()!=expected:
            raise ValueError('condensate source output changed')
    if not {'fort.7','fort.9','run.log'}.issubset(receipt['outputs']):
        raise ValueError('incomplete condensate source receipt')
    if input_fingerprint(source['executables']['tlusty'],directory,archived=True)!=receipt['input_sha256']:
        raise ValueError('condensate input fingerprint changed')
    if read_text(original(directory,'opacity.sha256')).strip()!=provenance['opacity_sha256']:
        raise ValueError('loaded condensate opacity hash differs')
    cfg=read_text(original(directory,'ember-condensates.cfg')).splitlines()
    if len(cfg)!=3 or cfg[0]!='equilibrium' or cfg[1]!=source['depletion']['fastchem']['source']:
        raise ValueError('condensate source mode or chemistry path differs')
    hashes=json.loads(read_text(original(directory,'condensates.sha256')))
    if hashes!=source['depletion']['fastchem']['chemistry_data_sha256']:
        raise ValueError('condensate chemistry identity differs')
    abundance,weights=composition(key[0],key[1],spec['metals'])
    metadata=json.loads((SOURCES/'synple-elements.json').read_text())
    numbers={};masses={}
    for s,n,w in zip(metadata['symbol'],abundance,weights,strict=True):
        if n>1e-90:numbers[s.capitalize()]=n;masses[s.capitalize()]=w*1.67333e-24/1.66053906660e-24
    declared={}
    for line in read_text(original(directory,'condensate-abundances.dat')).splitlines():
        if not line.strip() or line.startswith('#'):continue
        s,n=line.split()
        if s!='e-':declared[s]=10**(float(n)-12)
    if set(declared)!=set(numbers) or any(not math.isclose(declared[s],v,rel_tol=1e-12) for s,v in numbers.items()):
        raise ValueError('condensate abundance input differs from family mixture')
    inputs={k:read_text(original(directory,n)) for k,n in [
        ('atmosphere_input','fort.5'),('element_masses','ember-masses.dat'),
        ('parameters','tas'),('initial_structure','fort.8')]}
    if 'initialization_provenance' in provenance:
        proof=provenance['initialization_provenance'];initial=directory/'initializer'
        if not initial.exists():initial=Path(provenance['initial'])
        init_source=proof['prepared']
        if init_source.get('initialization_only',{}).get('canonical_prepared')!=source:
            raise ValueError('initializer belongs to a different canonical source')
        for name,expected_hash in [('completed.json',proof['source_receipt_sha256']),
                                   ('provenance.json',proof['provenance_sha256'])]:
            if hashlib.sha256(read_text(original(initial,name)).encode()).hexdigest()!=expected_hash:
                raise ValueError('initializer provenance changed')
        saved_init=json.loads(read_text(original(initial,'completed.json')))
        if input_fingerprint(init_source['executables']['tlusty'],initial,archived=True)!=saved_init['input_sha256']:
            raise ValueError('initializer input fingerprint changed')
        for name,expected_hash in saved_init['outputs'].items():
            if hashlib.sha256(read_text(original(initial,name)).encode()).hexdigest()!=expected_hash:
                raise ValueError('initializer output changed')
        old_t,old_g=map(float,read_text(original(initial,'fort.5')).splitlines()[0].split())
        guess=continuation_structure(read_text(original(initial,'fort.7')),old_t,old_g,key[2],key[3])
        if resample_initial_structure(guess,spec['depths'])!=inputs['initial_structure']:
            raise ValueError('canonical initial structure does not reproduce its recorded initializer')
    log=read_text(original(directory,'run.log'))
    source_inputs(inputs,spec,*key,log)
    state=source_state(log,read_text(original(directory,'fort.9')),key[2],key[3],
        {'temperature_K':temperatures(spec),'density_g_cm3':sequence(spec['log_density'])},spec['tau'])
    if state['depths']!=spec['depths']:raise ValueError('condensate depth count differs')
    profile=[]
    for line in log.rsplit('FINAL MODEL ATMOSPHERE',1)[1].splitlines():
        words=line.replace('D','E').split()
        if len(words)==11 and words[0].isdigit():profile.append(list(map(float,words)))
    profile.reverse()
    chemical=original(directory,'chemistry-audit-source.json')
    chemistry=json.loads(read_text(chemical))
    if chemistry['mode']!='equilibrium':raise ValueError('equilibrium chemistry audit required')
    hot_source=original(directory,'chemistry-hot-join-source.json')
    hot_join=json.loads(read_text(hot_source)) if any(r[3]>HOT_JOIN_TEMPERATURE for r in profile) else None
    audit,condensed=audit_profile(chemistry,profile,numbers,masses,hot_join)
    flux=np.array([r[9] for r in profile])
    audit['maximum_convective_flux_fraction_in_condensing_layers']=float(np.abs(flux[condensed]).max()) if condensed.any() else 0.
    supported=not bool(np.any(condensed & (np.abs(flux)>1e-8)))
    if not supported and not diagnostic_grain_transport:
        raise ValueError('condensing convective layers require grain enthalpy physics')
    # The interior boundary uses the undepleted bulk metal composition.
    # Both matching nodes and all deeper layers must therefore be free of
    # appreciable condensed particles, including in a radiative layer.
    optical_depth=np.array([r[2] for r in profile])
    if np.any(condensed & (optical_depth>=state['tau_bracket'][0])):
        raise ValueError('condensation at the matching depth requires a depleted interior composition')
    audit['grain_enthalpy_supported']=supported
    saved=json.loads(read_text(original(directory,'chemistry-audit.json')))
    if any(saved[k]!=v for k,v in audit.items()):raise ValueError('independent condensate audit differs')
    if saved['provenance']['chemistry_source_sha256']!=digest(chemical):
        raise ValueError('independent chemistry output checksum differs')
    raw_receipt=read_text(original(directory,'completed.json')).encode()
    if saved['provenance']['atmosphere_receipt_sha256']!=hashlib.sha256(raw_receipt).hexdigest():
        raise ValueError('chemistry audit belongs to another source atmosphere')
    if saved['provenance']['probe_sha256']!=source['depletion']['fastchem']['executable_sha256']:
        raise ValueError('chemistry probe identity differs')
    if hot_join is not None and saved['provenance'].get('hot_join_source_sha256')!=digest(hot_source):
        raise ValueError('vaporized-join chemistry checksum differs')
    return {'coordinates':key,'diagnostics':state,'chemistry':audit,'provenance':provenance}


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('model',type=Path)
    print(json.dumps(validate(p.parse_args().model),indent=2))
