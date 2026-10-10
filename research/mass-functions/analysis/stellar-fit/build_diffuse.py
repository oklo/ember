"""Diffuse association templates, using exact Galactic orbits and existing IMF draws."""
import json,time
import numpy as np
from common import W,V2,KMS
from galactic_diffusion import evolve

def population(cluster,birth):
    masses=[];maps=[];pairs=[];relative=[];system_mass=[];offset=0;soffset=0
    for repeat in range(2):
        for seed in [41,42,43]:
            name=f'pleiades_r4_fieldlow_s{seed}' if birth<.6 else f'pleiades_r3_assort_s{seed}'
            p=V2/'runs'/name;meta=json.loads((p/'initial.json').read_text());d=np.load(p/'initial.npz')
            m=d['mass'];nb=meta['nbinary'];n=len(m)
            sm=np.r_[m[:2*nb].reshape(-1,2).sum(1),m[2*nb:]]
            idx=np.r_[np.repeat(np.arange(nb),2),np.arange(nb,len(sm))]
            delta=np.zeros((n,3))
            dr=d['pos'][0:2*nb:2]-d['pos'][1:2*nb:2]
            # New isotropic orbital phases preserve the separation distribution.
            rng=np.random.default_rng(562+repeat*10+seed)
            unit=rng.normal(size=dr.shape);unit/=np.linalg.norm(unit,axis=1)[:,None]
            dr=unit*np.linalg.norm(dr,axis=1)[:,None]
            delta[0:2*nb:2]=dr*(m[1:2*nb:2]/sm[:nb])[:,None]
            delta[1:2*nb:2]=-dr*(m[0:2*nb:2]/sm[:nb])[:,None]
            pairs.extend(np.column_stack([np.arange(0,2*nb,2)+offset,np.arange(1,2*nb,2)+offset]))
            masses.extend(m);maps.extend(idx+soffset);relative.extend(delta);system_mass.extend(sm)
            offset+=n;soffset+=len(sm)
    return np.array(masses),np.array(maps),np.array(pairs),np.array(relative),np.array(system_mass)

if __name__=='__main__':
    (W/'diffuse').mkdir(exist_ok=True)
    plan=json.loads((W/'plan.json').read_text());receipts=[]
    for cl in ['hyades','pleiades']:
        name='hyades_r2.5_assort_s41' if cl=='hyades' else 'pleiades_r3_assort_s41'
        meta=json.loads((V2/'runs'/name/'initial.json').read_text())
        for birth in ([.63] if cl=='hyades' else [.63,.55]):
            m,mapper,pairs,delta,sm=population(cl,birth)
            nref=np.sum((m>=.1)&(m<.2));nbd=np.sum((m>=.03)&(m<.08))
            # Common standardized draws reduce sampling differences across sizes.
            rng=np.random.default_rng(476+int(birth*100));z=rng.normal(size=(len(sm),6))
            for size in plan['diffuse_initial_one_axis_size_pc']:
                for speed in plan['diffuse_initial_one_axis_velocity_kms']:
                    tag=f'{cl}_R{birth:g}_r{size}_v{speed:g}';p=W/'diffuse'/tag
                    if (W/'diffuse'/f'{tag}.npz').exists():continue
                    initial=np.array(meta['orbit']['start'])+z*np.array([size]*3+[speed*KMS]*3)
                    final,cost=evolve(initial,meta['age_Myr'])
                    rel=final[mapper,:3]+delta-np.array(meta['orbit']['present'])[:3]
                    vel=final[mapper,3:]-np.array(meta['orbit']['present'])[3:]
                    np.savez_compressed(str(p)+'.npz',mass=m,rel=rel,vel=vel,pairs=pairs)
                    info=dict(cluster=cl,expected_birth_ratio=birth,realized_birth_ratio=nbd/nref,size_pc=size,speed_kms=speed,
                              N_reference_components=int(nref),initial_mass_Msun=float(m.sum()),mass_per_reference_star=float(m.sum()/nref),
                              n_components=len(m),n_systems=len(sm),meta=meta,**cost)
                    (W/'diffuse'/f'{tag}.json').write_text(json.dumps(info,indent=2)+'\n')
                    receipts.append(dict(tag=tag,**cost));print(tag,cost,flush=True)
    # Preserve original measured costs when completed templates are reused.
    receipts=[]
    for path in sorted((W/'diffuse').glob('*.json')):
        info=json.loads(path.read_text())
        receipts.append(dict(tag=path.stem,**{key:info[key] for key in ['wall_seconds','n_orbits','rhs_calls','rtol']}))
    (W/'diffuse_execution.json').write_text(json.dumps(receipts,indent=2)+'\n')
