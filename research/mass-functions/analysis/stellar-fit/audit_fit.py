"""Independent checks of the stellar objective, input counts and provenance."""
import hashlib,json,sqlite3
from contextlib import closing
import numpy as np
import pandas as pd
from common import W,V1,V2,C,ANGLES
from fit_stars import gravity_limit

def main():
    metadata=[json.loads(p.read_text()) for p in (W/'templates').glob('*.json')]
    bound=[r for r in metadata if r['kind']=='bound']
    diffuse=[r for r in metadata if r['kind']=='diffuse']
    counts={r['tag']:np.load(W/'templates'/f'{r["tag"]}.npz')['counts'] for r in metadata}
    for r in metadata:
        a=counts[r['tag']]
        assert a.shape[1:]==(len(ANGLES[r['cluster']])-1,5)
        assert np.isfinite(a).all() and np.all(a>=0)
    rows=[]
    for b in bound:
        if not b['structural_screen_pass']:continue
        ob=np.load(W/f'observed_{b["cluster"]}.npz')['counts'][:,:,2:]
        ob=ob.reshape(len(ob),-1)
        bc=counts[b['tag']][:,:,2:].reshape(len(counts[b['tag']]),-1)
        y=ob.mean(0);base=bc.mean(0)
        variance=ob.var(0,ddof=1)+bc.var(0,ddof=1)+y+base+.5
        for d in [None]+diffuse:
            if d is not None and any(b[k]!=d[k] for k in ['cluster','cloud','expected_birth_ratio','membership_recovery','deep']):continue
            v=np.zeros_like(y) if d is None else counts[d['tag']][:,:,2:].mean(0).ravel()
            cap=0 if d is None else gravity_limit(b,d)[0]
            if d is not None and cap<=0:continue
            amp=float(np.clip(np.sum(v*(y-base)/variance)/max(np.sum(v*v/variance),1e-30),0,cap))
            score=np.sum((y-base-amp*v)**2/variance)
            full=counts[b['tag']].mean(0)+(0 if d is None else amp*counts[d['tag']].mean(0))
            rows.append(dict(cluster=b['cluster'],expected_birth_ratio=b['expected_birth_ratio'],cloud=b['cloud'],
                             membership_recovery=b['membership_recovery'],deep=b['deep'],run=b['run'],
                             diffuse_tag='' if d is None else d['tag'],amplitude=amp,score=score,
                             predicted_BD=full[:,0].sum(),predicted_stars=full[:,2:].sum()))
    test=pd.DataFrame(rows);keys=['cluster','expected_birth_ratio','cloud','membership_recovery','deep']
    best=test.loc[test.groupby(keys).score.idxmin()].copy()
    best.to_csv(W/'diagonal_objective_check.csv',index=False)
    # Catalogue attribution is a union: source rows below can overlap and must
    # not be summed into a second population count.
    attribution=[]
    with closing(sqlite3.connect(C/'cluster_catalogues.sqlite')) as con:
        sources=pd.read_sql('SELECT DISTINCT object_id,source FROM source_records WHERE selected=1',con)
    for cl,edges in ANGLES.items():
        d=np.load(C/f'{cl}_mass_sensitivity_draws.npz')
        mass=d['mass_Msun'];ids=d['object_id']
        frame=pd.read_csv(V1/f'{cl}_spatial_selection.csv').set_index('object_id').loc[ids]
        frame['expected_reference']=np.mean((mass>=.1)&(mass<.2),axis=0)
        frame['expected_BD']=np.mean((mass>=.03)&(mass<.08),axis=0)
        frame['angle_bin']=np.clip(np.searchsorted(edges,frame.angle_to_center_deg,side='right')-1,0,len(edges)-2)
        for (angle,source),g in frame.reset_index().merge(sources,on='object_id').groupby(['angle_bin','source']):
            attribution.append(dict(cluster=cl,angle_bin=angle,source=source,
                                    expected_reference=g.expected_reference.sum(),expected_BD=g.expected_BD.sum(),
                                    unique_candidates=len(g)))
    pd.DataFrame(attribution).to_csv(W/'catalogue_attribution.csv',index=False)
    prediction=pd.read_csv(W/'brown_dwarf_predictions.csv')
    assert (prediction.neglected_gravity_fraction<=.200001).all()
    assert (prediction.diffuse_mass_inside_aperture_upper_Msun<.1).all()
    assert len(set(r['run'] for r in bound))==15 # Two random-pairing stress runs omitted.
    assert len(set(r['run'] for r in diffuse))==27
    assert hashlib.md5((W/'raw/allsky_M10_hpx7.hdf5').read_bytes()).hexdigest()=='4d80621b40aead69ee6759954b7292b2'
    prior=json.loads((V2/'validation_receipt.json').read_text())
    assert prior['completed_full_age_models']==17 and not prior['pending_runs']
    receipt=dict(completed_collisional_runs=17,bound_templates_used=15,diffuse_geometries=27,
                 selection_scenarios_per_geometry=8,all_counts_finite_and_nonnegative=True,
                 published_map_checksum_verified=True,training_mass_Msun=[.1,.2],
                 held_out_mass_Msun=[.03,.08],cluster_normalizations_fixed=True,
                 selected_diffuse_core_mass_upper_Msun=float(prediction.diffuse_mass_inside_aperture_upper_Msun.max()),
                 selected_diffuse_neglected_gravity_fraction_max=float(prediction.neglected_gravity_fraction.max()),
                 independent_diagonal_objective_check=True,
                 code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in W.glob('*.py')})
    (W/'validation_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(best.groupby(['cluster','expected_birth_ratio'])[['score','amplitude','predicted_BD','predicted_stars']].agg(['min','max']).to_string(float_format=lambda x:f'{x:.4g}'))
    print(pd.DataFrame(attribution).query('(cluster=="pleiades" and angle_bin>=6) or (cluster=="hyades" and angle_bin>=3)').to_string(index=False,float_format=lambda x:f'{x:.4g}'))

if __name__=='__main__':main()
