"""Cluster-composition atmosphere pilots from published SPHINX I profiles.

Atmospheric [M/H] is interpolated between the published 0 and +0.25 nodes;
the interior uses scaled GS98 metals at the same [Fe/H]. This is not an
element-by-element mixture match. Y=0.27 is a modelling assumption. The
published atmospheric helium mixture remains fixed. These are pilot inputs,
not a validated cluster isochrone grid.
"""
from pathlib import Path
import json, re, math, hashlib, subprocess
import numpy as np

W=Path(__file__).resolve().parents[1]
ROOT=W.parents[1]
profiles={}
for path in (W/'raw/sphinx_atmospheres').glob('*_atms.txt'):
    m=re.search(r'Teff_([\d.]+)_logg_([\d.]+)_logZ_([+\-\d.]+)',path.name)
    T,g,z=map(float,m.groups());data=np.loadtxt(path)
    profiles[T,g,z]=(data[:,1]*1e6,data[:,0],path)
Ts=sorted({k[0] for k in profiles});Gs=sorted({k[1] for k in profiles})
cfg=ROOT/'out/luhman16-independent-oct8-v1/nc_m0.034_fresh/configuration.txt'
receipt=dict(stage='pilot boundary; cluster grid not yet validated',source='https://zenodo.org/records/7416042',
             source_metallicities=[0,.25],CO=.5,Teff_K=Ts,logg=Gs,
             pressure_cgs=1e8,assumptions=__doc__,cases={})
for name,feh in [('metalpoor',-.3),('metalrich',.3)]:
    ratio=.0231*10**feh;Y=.27;X=(1-Y)/(1+ratio);Z=ratio*X
    records=[];queries=[]
    low,high=(-.5,0.) if feh<0 else (0.,.5)
    fraction=(feh-low)/(high-low)
    for T in Ts:
        for g in Gs:
            pa,ta,fa=profiles[T,g,low];pb,tb,fb=profiles[T,g,high]
            assert pa[0]<1e8<=pa[-1] and pb[0]<1e8<=pb[-1]
            # Interpolate temperatures on common pressure coordinates. Both
            # bracketing source profiles cover every used pressure.
            lp=np.linspace(max(np.log(pa[0]),np.log(pb[0])),np.log(1e8),300)
            lt=(1-fraction)*np.interp(lp,np.log(pa),np.log(ta))+fraction*np.interp(lp,np.log(pb),np.log(tb))
            crossings=np.where(np.diff((lt>=np.log(T)).astype(int))==1)[0]
            assert len(crossings)>=1,(T,g)
            i=crossings[0];lpph=np.interp(np.log(T),lt[i:i+2],lp[i:i+2])
            lgp=np.linspace(lpph,np.log(1e8),80);tgrid=np.exp(np.interp(lgp,lp,lt));pgrid=np.exp(lgp)
            queries.extend(f'{t:.17g} {p:.17g}\n' for t,p in zip(tgrid,pgrid))
            records.append(dict(Teff=T,logg=g,Tmatch=float(np.exp(lt[-1])),logP=lgp,
                                files=[str(fa),str(fb)],Pphotosphere=float(np.exp(lpph))))
    cmd=[str(W/'physics/envelope_thermo'),str(cfg),str(X),str(Z)]
    result=subprocess.run(cmd,input=''.join(queries),capture_output=True,text=True,timeout=120)
    assert result.returncode==0,result.stderr
    lines=result.stdout.splitlines();assert len(lines)==len(queries)
    bad=[s for s in lines if s.startswith('ERROR')];assert not bad,bad[:4]
    state=np.array([list(map(float,s.split())) for s in lines]).reshape(len(records),80,6)
    assert np.all(state[:,:,2]>0) and np.all(state[:,:,3]>0)
    body=[]
    for record,a in zip(records,state):
        depth=float(np.trapezoid(a[:,1]/a[:,2]/10**record['logg'],record.pop('logP')))
        assert depth>0
        record['depth_cm']=depth
        body.append(f"{math.log10(record['Tmatch']):.17g} {depth:.17g}")
    dest=W/f'physics/sphinx_{name}_pm100.dat'
    dest.write_text('\n'.join(['EMBER_PRESSURE_MATCHED_ATMOSPHERE 1',
        f'scenario sphinxI_CO0.5_mh{feh:g}_GS98interior_Y0.27_pilot',
        'pressure_cgs 1e8',f'composition {X:.17g} {Z:.17g} 0.002 0.00001',
        f'log_teff {len(Ts)} '+' '.join(f'{math.log10(T):.17g}' for T in Ts),
        f'log_g {len(Gs)} '+' '.join(map(str,Gs)),'data']+body)+'\n')
    receipt['cases'][name]=dict(FeH=feh,X=X,Y=Y,Z=Z,file=str(dest),
        sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),nodes=records)
    print(name,'X,Y,Z',X,Y,Z,'nodes',len(records),flush=True)
(W/'physics/sphinx_boundary_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
