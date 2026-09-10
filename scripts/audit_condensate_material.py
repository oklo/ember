#!/usr/bin/env python3
"""Isolate opacity-temperature refinement with identical atmosphere equations."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

from generate_nongrey_grid import temperatures,sequence,input_fingerprint
from import_nongrey_grid import read_text,source_inputs,source_state,source_diagnostics_match
from nongrey_opacity import read_table
from prepare_nongrey_sources import digest
from validate_condensate_model import validate,original,PHYSICAL_KEYS


def thermal_sensitivity(path,reference,directory):
    """Check archived thermal counterfactuals for this exact fine model.

    This qualifies an interpolation diagnostic, never a production model or
    a bound on the missing physical EOS. Ten ppm is the requested diagnostic
    precision, substantially below the material interpolation criterion.
    """
    report=json.loads(Path(path).read_text());base=reference['provenance']
    if report['reference_source_sha256']!=hashlib.sha256(
            read_text(original(directory,'provenance.json')).encode()).hexdigest():
        raise ValueError('thermal controls belong to a different reference model')
    results=[]
    if sorted(r['extra_cp_erg_g_K'] for r in report['controls'])!=[1e6,1e7]:
        raise ValueError('both specified heat-capacity controls are required')
    for record in report['controls']:
        directory=Path(record['archive'])
        saved=json.loads((directory/'archive.json').read_text())
        if saved!=record:raise ValueError('thermal archive receipt differs')
        for name,expected in record['files_sha256'].items():
            if digest(directory/name)!=expected:raise ValueError('thermal archive bytes changed')
        source=json.loads(read_text(directory/'source-prepared.json.gz'))
        physics=source['initialization_only']
        if (physics['method']!='outer-layer heat-capacity sensitivity control; not an accepted atmosphere'
                or physics['canonical_prepared']!=base['prepared']
                or physics['extra_cp_erg_g_K']!=record['extra_cp_erg_g_K']
                or physics['temperature_taper_K']!=[2000,2400]):
            raise ValueError('thermal control uses different source physics')
        provenance=json.loads(read_text(directory/'provenance.json.gz'))
        if provenance['prepared']!=source or any(provenance[k]!=base[k] for k in
                ['specification','mode','XH','X3','teff_K','log_g','opacity_sha256']):
            raise ValueError('thermal control physical inputs differ')
        receipt=json.loads(read_text(directory/'completed.json.gz'))
        if input_fingerprint(source['executables']['tlusty'],directory,archived=True)!=receipt['input_sha256']:
            raise ValueError('thermal input fingerprint differs')
        for name,expected in receipt['outputs'].items():
            if hashlib.sha256(read_text(directory/(name+'.gz')).encode()).hexdigest()!=expected:
                raise ValueError('thermal source output differs')
        spec=base['specification'];coordinates=reference['coordinates']
        inputs={k:read_text(directory/(n+'.gz')) for k,n in [
            ('atmosphere_input','fort.5'),('element_masses','ember-masses.dat'),
            ('parameters','tas'),('initial_structure','fort.8')]}
        log=read_text(directory/'run.log.gz')
        source_inputs(inputs,spec,*coordinates,log)
        state=source_state(log,read_text(directory/'fort.9.gz'),coordinates[2],coordinates[3],
            {'temperature_K':temperatures(spec),'density_g_cm3':sequence(spec['log_density'])},spec['tau'])
        if not source_diagnostics_match(record['diagnostics'],state):
            raise ValueError('thermal diagnostics do not reproduce source outputs')
        relative={k:state[k]/reference['diagnostics'][k]-1 for k in ['T','Pgas','source_density']}
        if relative!=record['relative_boundary_difference']:raise ValueError('thermal comparison differs')
        if any(not np.isfinite(v) or abs(v)>1e-5 for v in relative.values()):
            raise ValueError('thermal sensitivity exceeds the diagnostic precision')
        results.append({'extra_cp_erg_g_K':record['extra_cp_erg_g_K'],
                        'relative_boundary_difference':relative,'archive_receipt_sha256':digest(directory/'archive.json')})
    if reference['chemistry']['maximum_convective_flux_fraction_in_condensing_layers']>1e-5:
        raise ValueError('material probe requires further grain transport investigation')
    return {'controls':results,'maximum_relative_boundary_difference':1e-5,
            'scope':'Thermal sensitivity of this interpolation diagnostic only; not a missing-physics bound or production grain-enthalpy acceptance.'}


def audit(coarse,fine,output,maximum_boundary_error=.005,*,thermal_controls=None):
    coarse=Path(coarse);fine=Path(fine)
    old=validate(coarse);new=validate(fine,diagnostic_grain_transport=thermal_controls is not None)
    thermal=thermal_sensitivity(thermal_controls,new,fine) if thermal_controls else None
    a,b=old['provenance'],new['provenance']
    if old['coordinates']!=new['coordinates'] or a['prepared']!=b['prepared']:
        raise ValueError('material refinement coordinates or equation source differ')
    if any(a['specification'].get(k)!=b['specification'].get(k) for k in PHYSICAL_KEYS
           if k not in ['log_temperature','temperature_K']):
        raise ValueError('material refinement changes other physical inputs')
    if temperatures(a['specification'])!=temperatures(b['specification'])[::2]:
        raise ValueError('refinement does not retain the original temperature nodes exactly')
    paths=[(d/'opacity.bin').resolve(strict=True) for d in [coarse,fine]]
    for path,provenance in zip(paths,[a,b],strict=True):
        if digest(path)!=provenance['opacity_sha256']:raise ValueError('opacity table changed')
    tables=[read_table(p) for p in paths]
    for key in ['abundance','ifmol','tmolim','flags','log_density','frequency']:
        if tables[0][key]!=tables[1][key]:raise ValueError('opacity inputs differ beyond temperature refinement')
    nf,nr,nt=tables[0]['shape'];new_nf,new_nr,new_nt=tables[1]['shape']
    if (new_nf,new_nr,new_nt)!=(nf,nr,2*nt-1):raise ValueError('unexpected refinement dimensions')
    original_values=np.asarray(tables[0]['log_opacity']).reshape(nf,nr,nt)
    retained_values=np.asarray(tables[1]['log_opacity']).reshape(nf,nr,new_nt)[:,:,::2]
    original_electrons=np.asarray(tables[0]['log_electron_density']).reshape(nt,nr)
    retained_electrons=np.asarray(tables[1]['log_electron_density']).reshape(new_nt,nr)[::2,:]
    if not np.array_equal(original_values,retained_values) or not np.array_equal(original_electrons,retained_electrons):
        raise ValueError('refinement changed original source values')
    relative={k:new['diagnostics'][k]/old['diagnostics'][k]-1 for k in ['T','Pgas','source_density']}
    error=max(abs(relative[k]) for k in ['T','Pgas'])
    report={'description':'Nested material-temperature refinement; original source values retained bit for bit',
            'coordinates':old['coordinates'],'coarse':str(coarse),'fine':str(fine),
            'temperature_rows':[nt,new_nt],'opacity_sha256':[a['opacity_sha256'],b['opacity_sha256']],
            'same_equation_executable_sha256':a['prepared']['executables']['tlusty'],
            'original_source_values_unchanged':True,'relative_boundary_difference':relative,
            'maximum_allowed_relative_boundary_error':maximum_boundary_error,'passed':error<=maximum_boundary_error,
            'note':'This tests temperature interpolation, not density/frequency resolution or a global physical uncertainty.'}
    if thermal:
        report['fine_model_production_accepted']=False
        report['fine_model_grain_enthalpy_supported']=new['chemistry']['grain_enthalpy_supported']
        report['thermal_sensitivity']=thermal
        report['note']+=' The fine model is a diagnostic only; production grid cells must separately pass the unchanged grain-enthalpy guard.'
    Path(output).write_text(json.dumps(report,indent=2)+'\n')
    if not report['passed']:raise ValueError('material temperature resolution requires review; report retained')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['coarse','fine','output']:p.add_argument(name,type=Path)
    p.add_argument('--maximum-boundary-error',type=float,default=.005)
    p.add_argument('--thermal-controls',type=Path,
                   help='qualify a diagnostic fine model using archived heat-capacity experiments; does not accept it into a production grid')
    a=p.parse_args()
    if not 0<a.maximum_boundary_error<1:raise ValueError('invalid comparison criterion')
    print(json.dumps(audit(a.coarse,a.fine,a.output,a.maximum_boundary_error,
                           thermal_controls=a.thermal_controls)['relative_boundary_difference']))
