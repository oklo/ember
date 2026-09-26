"""Recheck two low-gravity gas columns, including their depth sensitivity."""
import datetime
import json
from pathlib import Path
import sys
from generate_nongrey_grid import resample_initial_structure, truncate_initial_structure, temperatures, sequence
from import_nongrey_grid import source_inputs, source_state
from prepare_nongrey_sources import digest


def review(paths):
    records=[];hashes={};identity=None
    for path in paths:
        path=Path(path);plan=json.loads((path/'plan.json').read_text());work=path/'source'
        spec=plan['specification'];current=(spec,plan['prepared'],plan['opacity_sha256'])
        if identity is None:identity=current
        assert identity==current, 'physical or numerical inputs differ'
        assert digest(path/'run.py')==plan['script_sha256']
        assert digest(work/'opacity.bin')==plan['opacity_sha256']
        completed=json.loads((work/'completed.json').read_text())
        for name,sha in completed['outputs'].items():assert digest(work/name)==sha
        donor=Path(plan['initial_donor'])
        assert digest(donor/'fort.7')==plan['donor_checkpoint_sha256']
        assert digest(donor/'run.log')==plan['donor_log_sha256']
        initial=truncate_initial_structure((donor/'fort.7').read_text(),
            (donor/'run.log').read_text(),plan['initial_tau_bottom'])
        initial=resample_initial_structure(initial,spec['depths'])
        assert initial==(work/'fort.8').read_text(), 'initial column is not the declared transformation'
        log=(work/'run.log').read_text()
        source_inputs({key:(work/name).read_text() for key,name in
            [('atmosphere_input','fort.5'),('parameters','tas'),('element_masses','ember-masses.dat'),
             ('initial_structure','fort.8')]},spec,.7,0,3000,3.5,log)
        state=source_state(log,(work/'fort.9').read_text(),3000,3.5,
            dict(temperature_K=temperatures(spec),density_g_cm3=sequence(spec['log_density'])),100)
        records.append(dict(work=str(path),state=state,initial_tau_bottom=plan['initial_tau_bottom']))
        for name in ['plan.json','run.py','source/completed.json','source/fort.7','source/fort.9','source/run.log','source/fort.8']:
            hashes[str(path/name)]=digest(path/name)
    a,b=sorted(records,key=lambda r:r['state']['optical_depth_range'][1])
    difference={k:b['state'][k]/a['state'][k]-1 for k in ['T','Pgas','source_density']}
    checks=dict(source_checks_pass=True,
        independent_bottoms=b['state']['optical_depth_range'][1]/a['state']['optical_depth_range'][1]>1.5,
        temperature=abs(difference['T'])<1e-4,pressure=abs(difference['Pgas'])<1e-4)
    return dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        all_checks_pass=all(checks.values()),checks=checks,records=records,relative_depth_difference=difference,
        depth_comparison_criterion=1e-4,criterion_reference='scripts/audit_nongrey_lower_boundary.py default',
        selected_for_physical_PMS_track=False,accepted_as_gas_source=all(checks.values()),input_sha256=hashes,
        reviewer_sha256=digest(__file__),limitations=[
            'Gas-only, one effective temperature and gravity. Condensate removal, grains and a useful boundary grid remain to be assessed.',
            'Both full-resolution solves start from separately truncated versions of the same approximate coarse solution. Its unsupported deep tail was never accepted as a boundary.',
            'This comparison measures lower-boundary sensitivity at this cell, not complete wavelength or material-grid convergence.',
        ])


if __name__=='__main__':
    result=review(sys.argv[1:3]);Path(sys.argv[3]).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_checks_pass','checks','relative_depth_difference']},indent=2))
    raise SystemExit(0 if result['all_checks_pass'] else 1)
