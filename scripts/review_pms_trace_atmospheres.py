"""Revalidate PMS atmosphere controls and measure their trace-composition response.

Can report completed members of a running batch; uncompleted members remain
explicitly missing. This never accepts a frozen-composition stellar boundary.
"""
import datetime
import json
import math
from pathlib import Path
import sys

from generate_nongrey_grid import (
    completed_initial_structure, input_fingerprint, resample_initial_structure,
    sequence, temperatures,
)
from import_nongrey_grid import source_inputs, source_state
from prepare_nongrey_sources import digest


def review(work):
    work=Path(work).resolve()
    plan=json.loads((work/'plan.json').read_text())
    assert digest(work/'run.py')==plan['script_sha256']
    reference_path=Path(plan['reference_review'])
    assert digest(reference_path)==plan['reference_review_sha256']
    reference=json.loads(reference_path.read_text())
    assert reference['all_checks_pass']
    for path,expected in reference['input_sha256'].items():
        assert digest(path)==expected,path
    old=max(reference['records'],key=lambda r:r['state']['optical_depth_range'][1])
    donor=Path(old['work']).resolve()
    old_plan=json.loads((donor.parent/'plan.json').read_text())
    old_case=next(c for c in old_plan['cases']if c['name']==donor.name)
    prepared=plan['prepared']
    assert digest(prepared['tlusty'])==prepared['executables']['tlusty']
    initial=completed_initial_structure(donor,old_case['specification'],
        old['Teff_K'],old['log_g'],300)
    initial=resample_initial_structure(initial,300)
    hashes={str(p):digest(p)for p in [work/'plan.json',work/'run.py',reference_path,Path(__file__)]}
    records=[];missing=[]
    for case in plan['cases']:
        p=work/case['name']
        if not(p/'validated.json').exists():
            missing.append(case['name']);continue
        spec=case['specification'];T,g=case['teff_K'],case['log_g']
        assert (T,g)==(old['Teff_K'],old['log_g'])
        assert Path(case['donor']).resolve()==donor
        assert digest(donor/'completed.json')==case['donor_completed_sha256']
        assert digest(donor/'validated.json')==case['donor_validation_sha256']
        assert spec['depths']==300 and spec['atmosphere_frequencies']==20000
        for key,value in old_case['specification'].items():
            if key not in ('hydrogen','helium3','source'):assert spec[key]==value,key
        assert spec['hydrogen']==[case['XH']] and spec['helium3']==[case['X3']]
        assert digest(case['opacity'])==case['opacity_sha256']
        assert (p/'fort.8').read_text()==initial
        assert digest(p/'fort.8')==case['initial_sha256']
        completed=json.loads((p/'completed.json').read_text())
        assert input_fingerprint(prepared['executables']['tlusty'],p)==completed['input_sha256']
        for filename,expected in completed['outputs'].items():assert digest(p/filename)==expected
        validation=json.loads((p/'validated.json').read_text())
        assert validation['plan_sha256']==digest(work/'plan.json')
        log=(p/'run.log').read_text()
        source_inputs({k:(p/f).read_text()for k,f in [
            ('atmosphere_input','fort.5'),('parameters','tas'),
            ('element_masses','ember-masses.dat'),('initial_structure','fort.8')]},
            spec,case['XH'],case['X3'],T,g,log)
        state=source_state(log,(p/'fort.9').read_text(),T,g,
            dict(temperature_K=temperatures(spec),density_g_cm3=sequence(spec['log_density'])),100.)
        assert state==validation['diagnostics']
        records.append(dict(case=case['name'],XH=case['XH'],XHe3=case['X3'],
            Teff_K=T,log_g=g,state=state,
            relative_changes={k:state[k]/old['state'][k]-1 for k in ('T','Pgas','source_density')},
            source_cpu_seconds=validation['source_cpu_seconds']))
        for name in ['validated.json','completed.json','run.log','fort.5','fort.7','fort.8',
                     'fort.9','tas','ember-masses.dat','opacity.sha256','physics.json']:
            hashes[str(p/name)]=digest(p/name)
        hashes[case['opacity']]=case['opacity_sha256']
    result=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        all_completed_source_checks_pass=bool(records),complete_rectangle=not missing,
        missing_cases=missing,records=records,reference_state=old,
        selected_for_physical_PMS_track=False,input_sha256=hashes,
        limitations=[
            'Measures composition sensitivity at one temperature and gravity, not throughout the PMS atmosphere grid.',
            'Each comparison uses the same initial depth column; composition-dependent depth sensitivity is not independently measured.',
            'Fixed-composition boundary approximation is not selected by this measurement.',
            'The source approximates helium isotopes by their mean elemental mass; D chemistry and isotope band shifts are not resolved.'])
    if not missing:
        by={(r['XH'],r['XHe3']):r for r in records}
        result['response']={}
        for field in ('T','Pgas'):
            f0=old['state'][field]
            dH=math.log(by[(.6999,0.)]['state'][field]/f0)
            d3=math.log(by[(.7,.00012)]['state'][field]/f0)
            mixed=math.log(by[(.6999,.00012)]['state'][field]/f0)-dH-d3
            result['response'][field]=dict(
                dlnvalue_dXH=dH/(-.0001),dlnvalue_dXHe3=d3/.00012,
                mixed_log_change=mixed,
                note='Finite composition differences at this source column; mixed change measures departure from an additive linear response.')
    return result


if __name__=='__main__':
    result=review(sys.argv[1]);destination=Path(sys.argv[2])
    assert not destination.exists(),'use a fresh review path; do not replace a partial receipt'
    destination.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ['complete_rectangle','missing_cases','records']},indent=2))
