"""Pilot cluster boundaries: Sonora Diamondback below 2000 K, SPHINX above 2400 K.

Interpolate ln T at fixed gas pressure and [M/H], retaining separate clear and
f_sed=2 cloud scenarios. A smooth 2000--2400 K join is an explicit approximation
between different atmosphere calculations, not a new radiative-transfer solve.
Native Gibbs EOS supplies geometric depths and a deep adiabat where required.
Each extrapolated source profile has a withheld-pressure comparison recorded.
"""
from pathlib import Path
import io,json,zipfile,subprocess,hashlib,math
import numpy as np
from scipy.interpolate import RegularGridInterpolator

W=Path(__file__).resolve().parents[1];ROOT=W.parents[1]
cfg=ROOT/'out/luhman16-independent-oct8-v1/nc_m0.034_fresh/configuration.txt'
zip_path=ROOT/'out/substellar-cloud-benchmark-oct3-v1/pressure-temperature_profiles.zip'
base=json.loads((W/'physics/sphinx_boundary_receipt.json').read_text())
gravities=[31,100,316,1000,3160];Ts=np.arange(900,2401,100)
receipt=dict(source='https://doi.org/10.5281/zenodo.12735103',
    source_archive=str(zip_path),source_archive_sha256=hashlib.sha256(zip_path.read_bytes()).hexdigest(),
    assumptions=__doc__,cases={})
for cl in ['pleiades','hyades']:
    X,Z,feh=[base['cases'][cl][k] for k in ['X','Z','FeH']]
    for cloud in ['nc','f2']:
        records=[];requests=[]
        with zipfile.ZipFile(zip_path) as z:
            for t in Ts:
                for g in gravities:
                    for mh,label in [(0.,'0.0'),(.5,'+0.5')]:
                        name=f'pressure-temperature_profiles/t{t}g{g}{cloud}_m{label}_co1.0.pt'
                        raw=z.read(name);a=np.loadtxt(io.StringIO('\n'.join(raw.decode().splitlines()[2:])))
                        p,T=a[:,1]*1e6,a[:,2]
                        pstart=min(p[-1],1e8);tstart=np.exp(np.interp(np.log(pstart),np.log(p),np.log(T)))
                        trialP=p[-1]/1.2;trialT=np.exp(np.interp(np.log(trialP),np.log(p),np.log(T)))
                        requests.extend([(trialT,trialP,p[-1]),(tstart,pstart,1e8)])
                        records.append(dict(Teff=int(t),g=g,mh=mh,source=name,sha256=hashlib.sha256(raw).hexdigest(),
                                            source_Pmax_bar=float(p[-1]/1e6),p=p,T=T,deep_source_T=float(T[-1])))
        proc=subprocess.run([str(W/'physics/adiabat'),str(cfg),str(X),str(Z)],
            input=''.join(' '.join(f'{v:.17g}' for v in a)+'\n' for a in requests),text=True,capture_output=True,timeout=120)
        assert proc.returncode==0,proc.stderr
        lines=proc.stdout.splitlines();assert len(lines)==len(requests)
        assert not any(l.startswith('ERROR') for l in lines),[l for l in lines if l.startswith('ERROR')][:3]
        nodes={};diagnostics=[]
        for i,r in enumerate(records):
            check=float(lines[2*i].split()[0]);tm=float(lines[2*i+1].split()[0])
            p,T=r.pop('p'),r.pop('T');r['withheld_relative_T']=check/r.pop('deep_source_T')-1
            if p[-1]<1e8:
                p=np.append(p,1e8);T=np.append(T,tm)
            r['T100']=tm;nodes[r['Teff'],r['g'],r['mh']]=(p,T);diagnostics.append(r)
        # Interpolate source temperatures and evaluate native densities for
        # atmosphere depth from the T=Teff photosphere to the 100-bar match.
        requests=[];combined=[]
        for t in Ts:
            for g in gravities:
                pa,ta=nodes[t,g,0.];pb,tb=nodes[t,g,.5]
                lp=np.linspace(max(np.log(pa[0]),np.log(pb[0])),np.log(1e8),250)
                lt=(1-feh/.5)*np.interp(lp,np.log(pa),np.log(ta))+feh/.5*np.interp(lp,np.log(pb),np.log(tb))
                crossing=np.where(np.diff((lt>=np.log(t)).astype(int))==1)[0]
                assert len(crossing)>0,(t,g,cl,cloud)
                j=crossing[0];ph=np.interp(np.log(t),lt[j:j+2],lp[j:j+2]);pgrid=np.linspace(ph,np.log(1e8),80)
                tgrid=np.exp(np.interp(pgrid,lp,lt))
                requests.extend(zip(tgrid,np.exp(pgrid)))
                combined.append(dict(Teff=int(t),g=g,T100=float(np.exp(lt[-1])),pgrid=pgrid))
        proc=subprocess.run([str(W/'physics/envelope_thermo'),str(cfg),str(X),str(Z)],
            input=''.join(' '.join(f'{v:.17g}' for v in a)+'\n' for a in requests),text=True,capture_output=True,timeout=120)
        assert proc.returncode==0,proc.stderr
        lines=proc.stdout.splitlines();assert len(lines)==len(requests)
        assert not any(l.startswith('ERROR') for l in lines),[l for l in lines if l.startswith('ERROR')][:3]
        a=np.array([list(map(float,l.split())) for l in lines]).reshape(len(combined),80,6)
        assert (a[:,:,2]>0).all() and (a[:,:,3]>0).all()
        for r,q in zip(combined,a):r['depth']=float(np.trapezoid(q[:,1]/q[:,2]/(r['g']*100),r.pop('pgrid')))
        coolpath=W/f'physics/sonora_{cl}_{cloud}_pm100_v2.dat'
        coolpath.write_text('\n'.join(['EMBER_PRESSURE_MATCHED_ATMOSPHERE 1',f'scenario cluster_{cl}_{cloud}_sonora_only_pilot',
            'pressure_cgs 1e8',f'composition {X:.17g} {Z:.17g} 0.002 0.00001',
            f'log_teff {len(Ts)} '+' '.join(f'{np.log10(t):.17g}' for t in Ts),
            f'log_g {len(gravities)} '+' '.join(f'{np.log10(g)+2:.17g}' for g in gravities),'data']+
            [f"{np.log10(r['T100']):.17g} {r['depth']:.17g}" for r in combined])+'\n')
        # Evaluate both source families within their native gravity ranges.
        warm=base['cases'][cl]['nodes'];wt=np.array(base['Teff_K']);wg=np.array(base['logg'])
        wi=RegularGridInterpolator((np.log10(wt),wg),np.array([[np.log10(r['Tmatch']),r['depth_cm']] for r in warm]).reshape(len(wt),len(wg),2))
        output=[];join=[];allT=np.arange(900,4001,100);allg=np.log10(gravities[1:])+2
        for t in allT:
            for k,g in enumerate(allg):
                v=None
                if t<=2400:
                    r=combined[int((t-900)/100)*len(gravities)+k+1];v=np.array([np.log10(r['T100']),r['depth']])
                if t>=2000:
                    warmv=wi([[np.log10(t),g]])[0]
                    if t<2400:
                        z=(t-2000)/400;weight=z*z*(3-2*z)
                        join.append(dict(Teff=int(t),logg=float(g),T100_ratio=float(10**(warmv[0]-v[0])),depth_ratio=float(warmv[1]/v[1])))
                        v=(1-weight)*v+weight*warmv
                    else:v=warmv
                output.append(v)
        path=W/f'physics/cluster_{cl}_{cloud}_pm100_v2.dat'
        path.write_text('\n'.join(['EMBER_PRESSURE_MATCHED_ATMOSPHERE 1',f'scenario cluster_{cl}_{cloud}_sphinx_pilot',
            'pressure_cgs 1e8',f'composition {X:.17g} {Z:.17g} 0.002 0.00001',
            f'log_teff {len(allT)} '+' '.join(f'{np.log10(t):.17g}' for t in allT),
            f'log_g {len(allg)} '+' '.join(f'{g:.17g}' for g in allg),'data']+
            [f'{v[0]:.17g} {v[1]:.17g}' for v in output])+'\n')
        receipt['cases'][cl+'_'+cloud]=dict(file=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            cool_only_file=str(coolpath),cool_only_sha256=hashlib.sha256(coolpath.read_bytes()).hexdigest(),
            FeH=feh,X=X,Z=Z,source_profiles=diagnostics,join_comparisons=join,
            withheld_max_relative_T=max(abs(r['withheld_relative_T']) for r in diagnostics))
        print(cl,cloud,len(output),'nodes', 'max withheld T error',receipt['cases'][cl+'_'+cloud]['withheld_max_relative_T'],flush=True)
(W/'physics/cool_boundary_receipt_v2.json').write_text(json.dumps(receipt,indent=2)+'\n')
