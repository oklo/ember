#!/usr/bin/env python3
"""Summarize a completed forward track and optional numerical refinements.

The requested age and immutable-run receipts are checked explicitly. These
are segment comparisons, not hydrogen-exhaustion lifetime uncertainties.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def record(path,age):
    path=Path(path);raw=path.read_bytes();data=json.loads(raw)
    if not data['converged'] or data['history'][-1][0]!=age:raise ValueError('requested age was not reached')
    if any(not math.isfinite(v) for key in ['history','profile'] for row in data[key] for v in row):
        raise ValueError('nonfinite stellar history or final profile')
    receipt=json.loads(path.with_suffix('.receipt.json').read_text())
    checksum=hashlib.sha256(raw).hexdigest()
    if (receipt['returncode'] or receipt['data_changed_during_run'] or receipt.get('restart_changed_during_run')
            or receipt['output_sha256']!=checksum):
        raise ValueError('immutable run receipt does not validate the output')
    rows=[dict(zip(data['columns'],row,strict=True)) for row in data['history']]
    accepted=rows[1:]
    if not accepted:raise ValueError('no accepted evolutionary steps')
    peak=max(rows,key=lambda r:r['central_Y3'])
    result={'input':str(path),'sha256':checksum,'points':data['points'],
        'step_error_tolerances':data['step_error_tolerances'],'accepted_macrosteps':len(accepted),
        'rejected_macrosteps':data['rejected_steps'],'final':rows[-1],
        'peak_He3':{'age_yr':peak['age_yr'],'mass_fraction':peak['central_Y3']},
        'max_last_halfstep_luminosity_imbalance':max(abs(r['discrete_luminosity_balance']) for r in accepted),
        'max_last_halfstep_nuclear_mass_imbalance':max(abs(r['nuclear_rest_mass_balance']) for r in accepted),
        'minimum_accepted_convective_mass_fraction':min(r['convective_mass_fraction'] for r in accepted),
        'provenance':receipt}
    if 'restart' in data:
        restart=data['restart']
        if rows[0]['age_yr']!=restart['history_start_age_yr'] or not receipt.get('restart_input_sha256'):
            raise ValueError('restart history or input provenance is incomplete')
        result['history_interval_yr']=[rows[0]['age_yr'],rows[-1]['age_yr']]
        result['total_accepted_macrosteps']=len(accepted)+restart['accepted_steps_before_restart']
        result['total_rejected_macrosteps']=data['rejected_steps']
        result['rejected_macrosteps']=data['rejected_steps']-restart['rejected_steps_before_restart']
        result['peak_He3']['scope']='maximum in this recorded continuation segment only'
    return data,result


def difference(first,second):
    a,b=first['final'],second['final']
    if a['age_yr']!=b['age_yr']:raise ValueError('different comparison ages')
    return {'reference':first['input'],'comparison':second['input'],
            'relative':{k:b[k]/a[k]-1 for k in ['R_Rsun','L_Lsun','Teff_K','central_T_K','central_rho']},
            'absolute':{k:b[k]-a[k] for k in ['central_X','central_Y3']}}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('reference',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--age',type=float,required=True);p.add_argument('--mesh',type=Path);p.add_argument('--tight',type=Path)
    p.add_argument('--plot',action='store_true')
    a=p.parse_args();data,base=record(a.reference,a.age)
    compact={k:v for k,v in data.items() if k not in ['history','profile','profile_columns','columns']}
    compact.update(base)
    names=['age_yr','R_Rsun','L_Lsun','Teff_K','central_X','central_Y3','central_T_K','central_rho']
    indices=[data['columns'].index(k) for k in names]
    compact['history_columns']=names;compact['history']=[[r[i] for i in indices] for r in data['history']]
    compact['scope']='Completed main-sequence segment; no hydrogen-exhaustion or lifetime convergence claim.'
    a.output.parent.mkdir(parents=True,exist_ok=True)
    checks={'description':'Same-physics numerical segment comparisons, not physical error bars','runs':[base]}
    for label,path in [('mesh_comparison',a.mesh),('timestep_comparison',a.tight)]:
        if path is None:continue
        other,result=record(path,a.age)
        keys=['mass_Msun','mass_basis','nuclear_model','transport_model','atmosphere_model','eos_model',
              'convection_criterion','secular_mixing','metal_inventory','age_origin']
        if any(data.get(k)!=other.get(k) for k in keys):raise ValueError('comparison changes physical prescriptions')
        if data['history'][0][0]!=other['history'][0][0]:raise ValueError('comparison changes the recorded starting age')
        if base['provenance']['executable_sha256']!=result['provenance']['executable_sha256']:
            raise ValueError('comparison changes the executable')
        before,after=base['provenance'],result['provenance']
        if any(after['data_sha256'].get(k)!=v for k,v in before['data_sha256'].items()):
            raise ValueError('comparison changes reference input data')
        allowed={0,1 if label=='mesh_comparison' else 4}
        if (before['working_directory']!=after['working_directory'] or
                len(before['command'])!=len(after['command']) or
                any(x!=y for i,(x,y) in enumerate(zip(before['command'],after['command'])) if i not in allowed)):
            raise ValueError('comparison changes other run arguments')
        if label=='mesh_comparison' and (other['points']<=data['points'] or other['step_error_tolerances']!=data['step_error_tolerances']):
            raise ValueError('mesh comparison must refine only the mesh')
        if label=='timestep_comparison' and (other['points']!=data['points'] or
                any(v>=data['step_error_tolerances'][k] for k,v in other['step_error_tolerances'].items())):
            raise ValueError('timestep comparison must tighten tolerances at fixed mesh')
        checks['runs'].append(result);checks[label]=difference(base,result)
    a.output.write_text(json.dumps(compact,indent=2,allow_nan=False)+'\n')
    if len(checks['runs'])>1:
        a.output.with_name(a.output.stem+'_convergence.json').write_text(json.dumps(checks,indent=2,allow_nan=False)+'\n')
    if a.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
        fig,axes=plt.subplots(2,2,figsize=(9,6),sharex=True,constrained_layout=True)
        values={k:[r[data['columns'].index(k)] for r in data['history']] for k in names}
        ages=[t/1e12 for t in values['age_yr']]
        for k,label,color in [('central_X','Hydrogen-1','#b34b24'),('central_Y3','Helium-3','#246e91')]:
            axes[0,0].plot(ages,values[k],label=label,color=color)
        profile=dict(zip(data['profile_columns'],data['profile'][0],strict=True))
        nonmetal=profile['X']+profile['Y3']+profile['Y4']
        axes[0,0].plot(ages,[nonmetal-x-y for x,y in zip(values['central_X'],values['central_Y3'])],label='Helium-4',color='#657344')
        axes[0,0].set_ylabel('Central baryonic mass fraction');axes[0,0].legend(frameon=False)
        for ax,k,label,scale in [(axes[0,1],'L_Lsun',r'Luminosity ($10^{-3}L_\odot$)',1000),
                                (axes[1,0],'R_Rsun',r'Radius ($R_\odot$)',1),
                                (axes[1,1],'Teff_K','Effective temperature (K)',1)]:
            ax.plot(ages,[scale*v for v in values[k]],color='#246e91');ax.set_ylabel(label)
        for ax in axes.flat:ax.grid(alpha=.18);ax.set_xlim(ages[0],a.age/1e12)
        for ax in axes[1]:ax.set_xlabel('Elapsed time (trillion years)')
        fig.suptitle(f"{data['mass_Msun']:g} solar mass · {data['points']} mass points · specified main-sequence initial model\n"
                     'Composition-dependent EOS and non-grey boundary · SFII/SVH · Ledoux transport')
        for suffix in ['.png','.pdf']:fig.savefig(a.output.with_suffix(suffix),dpi=180)
    print(json.dumps({'age_yr':a.age,'final':base['final'],'accepted_steps':base['accepted_macrosteps']}))


if __name__=='__main__':main()
