"""Review bounded MESA controls without selecting a production approximation."""
from pathlib import Path
import argparse,hashlib,json,re
from datetime import datetime,timezone
import numpy as np
from write_scientific_result import write_result

def history(p):
    ls=p.read_text().splitlines();keys=ls[5].split();rows=[]
    for line in ls[6:]:
        parts=line.split()
        if len(parts)!=len(keys):continue
        row=dict(zip(keys,map(float,parts)))
        while rows and rows[-1]['model_number']>=row['model_number']:rows.pop()
        rows.append(row)
    return rows

def review(d):
    plan=json.loads((d/'plan.json').read_text())
    done=json.loads((d/'completion.json').read_text())
    assert done['returncode']==0 and not done['budget_timeout']
    for p,h in plan['input_sha256'].items(): assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
    log=(d/'run.log').read_text();assert 'termination code: max_model_number' in log
    rs=history(d/'LOGS/history.data')
    assert len(rs)==plan['last_model']-plan['start_model']
    assert rs[0]['model_number']==plan['start_model']+1 and rs[-1]['model_number']==plan['last_model']
    assert all(b['star_age']>=a['star_age'] for a,b in zip(rs,rs[1:]))
    # MESA checks the administrative model-number stop before the hard
    # temperature-step check. Preserve and report such terminal rows, but
    # exclude them from the checked continuation statistics.
    raw_count=len(rs);excluded=[]
    bad=[i+1 for i,d in enumerate(np.diff([r['log_Teff'] for r in rs]))
         if abs(d)>.02+1e-12]
    if bad:
        assert bad==[len(rs)-1],bad
        excluded=[rs[-1]];rs=rs[:-1]
    t=np.array([10**r['log_Teff'] for r in rs]);dl=np.diff([r['log_Teff'] for r in rs]); assert max(abs(dl))<=.02+1e-12
    profiles=[]
    for p in sorted((d/'LOGS').glob('profile*.data')):
        ls=p.read_text().splitlines();head=dict(zip(ls[1].split(),ls[2].split()));keys=ls[5].split();x=np.loadtxt(p,skiprows=6)
        assert np.isfinite(x).all()
        if int(head['model_number'])>rs[-1]['model_number']:continue
        get=lambda k:x[:,keys.index(k)]
        dm=-np.diff(np.r_[get('mass'),0])*1.9884e33
        profiles.append({'model_number':int(head['model_number']),'surface_X':get('x_mass_fraction_H')[0],'surface_Z':get('z_mass_fraction_metals')[0],'diffusion_heat_to_surface_luminosity':float(np.sum(dm*get('eps_diffusion'))/(3.828e33*get('luminosity')[0])),'unmixed_mass_fraction':float(sum(dm[get('log_D_mix')<=5])/sum(dm))})
    cpu=[float(x) for x in re.findall(r'(?m)^(?:user|sys)\s+([\d.]+)$',log)]
    return {'case':d.name,'rows':len(rs),'raw_saved_rows':raw_count,'terminal_rows_excluded_from_statistics':excluded,'elapsed_age_years':rs[-1]['star_age']-rs[0]['star_age'],'Teff_range_K':[float(min(t)),float(max(t))],'largest_adjacent_Teff_change_K':float(max(abs(np.diff(t)))),'retries':rs[-1]['num_retries']-rs[0]['num_retries'],'active_diffusion_saved_states':sum(r['diffusion_solver_steps']>0 for r in rs),'wall_seconds':done['wall_seconds'],'process_cpu_seconds':sum(cpu),'endpoint':rs[-1],'profiles':profiles,'input_sha256':{str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in [d/'plan.json',d/'completion.json',d/'LOGS/history.data',d/'run.log']}}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);p.add_argument('cases',nargs='+',type=Path);args=p.parse_args()
    results=[review(d) for d in args.cases]
    report={'created_utc':datetime.now(timezone.utc).isoformat(),'cases':results,'selected_for_evolution':False,'limitations':['Equal starting model and accepted-model cap; different achieved ages are reported, not compared as matched endpoints.','Short controls locate numerical sensitivities, not global temporal convergence.','Removing or automatically suppressing diffusion is a diagnostic, not proof that diffusion is physically negligible at every later stage.'],'review_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    write_result(args.output,report)
    for r in results:print(r['case'],r['Teff_range_K'],r['elapsed_age_years']/1e6,'Myr',r['retries'],'retries',r['active_diffusion_saved_states'],'active diffusion states')
