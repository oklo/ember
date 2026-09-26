"""Reconstruct and compare full-resolution PMS atmosphere depth pairs."""
import datetime,json,sys
from pathlib import Path
from generate_nongrey_grid import (completed_initial_structure,truncate_initial_structure,
    scale_initial_column,resample_initial_structure,input_fingerprint,temperatures,sequence,checkpoint_initial_structure)
from import_nongrey_grid import source_inputs,source_state
from prepare_nongrey_sources import digest


def review(work,prefix,other=None):
    candidates=[];hashes={}
    for root in [Path(work)]+([Path(other)]if other else []):
        plan=json.loads((root/'plan.json').read_text())
        assert digest(root/'run.py')==plan['script_sha256']
        hashes[str(root/'plan.json')]=digest(root/'plan.json');hashes[str(root/'run.py')]=digest(root/'run.py')
        assert digest(plan['opacity'])==plan['opacity_sha256']
        assert digest(plan['prepared']['tlusty'])==plan['prepared']['executables']['tlusty']
        candidates.extend((root,plan,c)for c in plan['cases']if c['name'].startswith(prefix+'_') and (root/c['name']/'validated.json').exists())
    assert len(candidates)==2, 'both completed source members required'
    assert candidates[0][2]['specification']==candidates[1][2]['specification']
    assert candidates[0][1]['prepared']==candidates[1][1]['prepared']
    assert candidates[0][1]['opacity_sha256']==candidates[1][1]['opacity_sha256']
    records=[]
    for work,plan,c in candidates:
        d=work/c['name'];spec=c['specification'];T,g=c['teff_K'],c['log_g']
        assert spec['depths']==300 and spec['atmosphere_frequencies']==20000
        completed=json.loads((d/'completed.json').read_text())
        assert input_fingerprint(plan['prepared']['executables']['tlusty'],d)==completed['input_sha256']
        for name,expected in completed['outputs'].items():assert digest(d/name)==expected
        assert json.loads((d/'validated.json').read_text())['plan_sha256']==digest(work/'plan.json')
        donor=Path(c['donor']);old=json.loads((donor.parent/'plan.json').read_text())
        assert digest(donor.parent/'run.py')==old['script_sha256']
        donor_receipts=[donor.parent/'plan.json']
        if c.get('initial_kind')=='attested_interrupted':
            assert digest(donor/'checkpoint.json')==c['donor_checkpoint_sha256']
            initial=checkpoint_initial_structure(donor,.7,0,T,g)
            donor_receipts.append(donor/'checkpoint.json')
        else:
            assert digest(donor/'completed.json')==c['donor_completed_sha256']
            assert digest(donor/'validated.json')==c['donor_validation_sha256']
            previous=next(x for x in old['cases']if x['name']==donor.name)
            initial=completed_initial_structure(donor,previous['specification'],T,g,300)
            if c['truncate_tau'] is not None:
                initial=truncate_initial_structure(initial,(donor/'run.log').read_text(),c['truncate_tau'])
            initial=resample_initial_structure(scale_initial_column(initial,c['initial_column_factor']),300)
            donor_receipts.extend([donor/'completed.json',donor/'validated.json'])
        # atmosphere_inputs serializes the saved structure even when its
        # depth count is unchanged. Reconstruct that exact input formatting.
        initial=resample_initial_structure(initial,spec['depths'])
        assert initial==(d/'fort.8').read_text()
        assert digest(d/'fort.8')==c['initial_sha256']
        log=(d/'run.log').read_text()
        source_inputs({key:(d/name).read_text()for key,name in
            [('atmosphere_input','fort.5'),('parameters','tas'),('element_masses','ember-masses.dat'),
             ('initial_structure','fort.8')]},spec,.7,0,T,g,log)
        state=source_state(log,(d/'fort.9').read_text(),T,g,
            dict(temperature_K=temperatures(spec),density_g_cm3=sequence(spec['log_density'])),100)
        records.append(dict(work=str(d),state=state,Teff_K=T,log_g=g))
        for path in [d/name for name in ['completed.json','validated.json','fort.5','fort.7','fort.8','fort.9','run.log','tas','ember-masses.dat','physics.json']]+donor_receipts:
            hashes[str(path)]=digest(path)
    shallow,deep=sorted(records,key=lambda r:r['state']['optical_depth_range'][1])
    differences={key:deep['state'][key]/shallow['state'][key]-1 for key in ['T','Pgas','source_density']}
    checks=dict(source_checks_pass=True,independent_bottoms=deep['state']['optical_depth_range'][1]/shallow['state']['optical_depth_range'][1]>1.5,
        temperature=abs(differences['T'])<1e-4,pressure=abs(differences['Pgas'])<1e-4)
    return dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        all_checks_pass=all(checks.values()),checks=checks,records=records,relative_depth_difference=differences,
        depth_comparison_criterion=1e-4,criterion_reference='scripts/audit_nongrey_lower_boundary.py default',
        accepted_as_gas_source=all(checks.values()),selected_for_physical_PMS_track=False,
        input_sha256=hashes,reviewer_sha256=digest(__file__),limitations=[
            'Gas source at one fixed composition. Actual-profile condensation and evolving isotope coverage remain separate checks.',
            'Both starting columns are reconstructed from the same converged coarse source; each solves afresh at full resolution.',
            'Independent bottom-depth sensitivity is tested, not complete opacity or wavelength accuracy.'])


if __name__=='__main__':
    result=review(sys.argv[1],sys.argv[2],sys.argv[4]if len(sys.argv)==5 else None);Path(sys.argv[3]).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ['all_checks_pass','checks','relative_depth_difference']},indent=2))
    raise SystemExit(0 if result['all_checks_pass']else 1)
