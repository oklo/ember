"""Extend cluster test boundaries using existing Sonora Bobcat profiles.

Bobcat supplies the clear boundary below 900 K. Between 900 and 1200 K,
blend log(T100) and geometric thickness into the original Diamondback grid.
For the f2 branch this is an explicit cloud-clearing scenario, not a cloudy
calculation below Diamondback's published range. Run fresh tracks through the
whole transition. This is a sensitivity model, not a fitted cloud prescription.
Use the published density for geometric thickness; where a profile stops short
of 100 bar, extend along the native adiabat and use the native endpoint density.
Atmospheric helium remains the published value, so this geometry is not an
exact match to the chosen interior mixture.
"""
from pathlib import Path
import io, re, json, tarfile, subprocess, hashlib
import numpy as np
from scipy.interpolate import RegularGridInterpolator

W=Path(__file__).resolve().parents[1]; ROOT=W.parents[1]
source=ROOT/'out/sonora-bobcat-cool-oct4-v1'
cfg=ROOT/'out/luhman16-independent-oct8-v1/nc_m0.034_fresh/configuration.txt'
base=json.loads((W/'physics/sphinx_boundary_receipt.json').read_text())
gs=[31,100,316,1000,3160]; profiles={}; archives=[]
for mh in (0.,.5):
    p=source/f'structures_m+{mh:.1f}.tar.gz'
    archives.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    with tarfile.open(p) as tar:
        for member in tar:
            m=re.search(r't(\d+)g(\d+)nc_',member.name)
            if not m: continue
            t,g=map(int,m.groups())
            if not 500<=t<1200 or g not in gs: continue
            raw=tar.extractfile(member).read()
            a=np.loadtxt(io.StringIO(raw.decode()),skiprows=1)
            profiles[t,g,mh]=(a[:,1]*1e6,a[:,2],a[:,6],member.name,hashlib.sha256(raw).hexdigest())
ts=sorted({t for t,g,m in profiles if all((t,g,z) in profiles for g in gs for z in (0.,.5))})
assert ts[0]==500 and ts[-1]==1100

def run_helper(helper,queries,X,Z):
    p=subprocess.run([str(W/'physics'/helper),str(cfg),str(X),str(Z)],
        input=''.join(' '.join(f'{x:.17g}' for x in r)+'\n' for r in queries),
        text=True,capture_output=True,timeout=120)
    assert p.returncode==0,p.stderr
    rows=p.stdout.splitlines(); assert len(rows)==len(queries)
    assert not any(r.startswith('ERROR') for r in rows),[r for r in rows if r.startswith('ERROR')][:5]
    return np.array([list(map(float,r.split())) for r in rows])

def read_table(path):
    a=path.read_text().splitlines()
    tx=np.array([float(x) for x in a[4].split()[2:]])
    gx=np.array([float(x) for x in a[5].split()[2:]])
    return tx,gx,np.loadtxt(io.StringIO('\n'.join(a[7:]))).reshape(len(tx),len(gx),2)

receipt=dict(source='https://doi.org/10.5281/zenodo.5063476',archives=archives,
             assumptions=__doc__,cases={})
for cl in ('pleiades','hyades'):
    X,Z,feh=[base['cases'][cl][k] for k in ('X','Z','FeH')]
    requests=[]; nodes=[]
    for t in ts:
        for g in gs:
            for mh in (0.,.5):
                p,T,rho,name,h=profiles[t,g,mh]
                match=min(1e8,p[-1]);start=np.exp(np.interp(np.log(match),np.log(p),np.log(T)))
                requests.append((start,match,1e8));nodes.append((t,g,mh,p,T,rho,name,h))
    deep=run_helper('adiabat',requests,X,Z)
    endpoint=run_helper('envelope_thermo',[(a[0],1e8) for a in deep],X,Z)
    expanded={}; source_rows=[]
    for row,ans,end in zip(nodes,deep,endpoint):
        t,g,mh,p,T,rho,name,h=row
        if p[-1]<1e8:p=np.append(p,1e8);T=np.append(T,ans[0]);rho=np.append(rho,end[2])
        expanded[t,g,mh]=(p,T,rho)
        source_rows.append(dict(Teff=t,g_m_s2=g,metallicity=mh,name=name,sha256=h,
                                adiabatic_extension=bool(profiles[t,g,mh][0][-1]<1e8)))
    bc=[]
    for t in ts:
        for g in gs:
            pa,ta,ra=expanded[t,g,0.];pb,tb,rb=expanded[t,g,.5]
            lp=np.linspace(max(np.log(pa[0]),np.log(pb[0])),np.log(1e8),300)
            lt=(1-feh/.5)*np.interp(lp,np.log(pa),np.log(ta))+feh/.5*np.interp(lp,np.log(pb),np.log(tb))
            lr=(1-feh/.5)*np.interp(lp,np.log(pa),np.log(ra))+feh/.5*np.interp(lp,np.log(pb),np.log(rb))
            crossings=np.where(np.diff((lt>=np.log(t)).astype(int))==1)[0]
            assert len(crossings),(cl,t,g)
            j=crossings[0];photosphere=np.interp(np.log(t),lt[j:j+2],lp[j:j+2])
            lgp=np.linspace(photosphere,np.log(1e8),400)
            depth=float(np.trapezoid(np.exp(lgp-np.interp(lgp,lp,lr))/(g*100),lgp))
            assert depth>0
            bc.append([lt[-1]/np.log(10),depth])
    bi=RegularGridInterpolator((np.log10(ts),np.log10(gs)+2),np.array(bc).reshape(len(ts),len(gs),2))
    for cloud in ('nc','f2'):
        old=W/f'physics/sonora_{cl}_{cloud}_pm100_v2.dat';tx,gx,v=read_table(old)
        warm=RegularGridInterpolator((tx,gx),v)
        allT=sorted(set(ts)|set(np.round(10**tx).astype(int)));out=[];join=[]
        for t in allT:
            for g in gx:
                if t<1200:
                    b=bi([[np.log10(t),g]])[0]
                if t<=900: value=b
                elif t>=1200: value=warm([[np.log10(t),g]])[0]
                else:
                    a=warm([[np.log10(t),g]])[0];u=(t-900)/300;s=u*u*(3-2*u)
                    value=(1-s)*b+s*a
                    join.append(dict(Teff=t,logg=float(g),warm_to_bobcat_T100=float(10**(a[0]-b[0]))))
                out.append(value)
        p=W/f'physics/sonora_{cl}_{cloud}_cool500.dat'
        p.write_text('\n'.join(['EMBER_PRESSURE_MATCHED_ATMOSPHERE 1',
            f'scenario cluster_{cl}_{cloud}_bobcat_below900_cloudclearing_control',
            'pressure_cgs 1e8',f'composition {X:.17g} {Z:.17g} 0.002 0.00001',
            f'log_teff {len(allT)} '+' '.join(f'{np.log10(t):.17g}' for t in allT),
            f'log_g {len(gx)} '+' '.join(f'{g:.17g}' for g in gx),'data']+
            [f'{a:.17g} {b:.17g}' for a,b in out])+'\n')
        receipt['cases'][cl+'_'+cloud]=dict(file=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
            X=X,Z=Z,FeH=feh,join=join,source_profiles=source_rows)
        wx,wg,wv=read_table(W/f'physics/cluster_{cl}_{cloud}_pm100_v2.dat')
        cv=np.array(out).reshape(len(allT),len(gx),2)
        assert np.allclose(gx[1:],wg,rtol=0,atol=1e-12)
        high=wx>np.log10(2400)+1e-12
        finalx=np.concatenate([np.log10(allT),wx[high]])
        finalv=np.concatenate([cv[:,1:,:],wv[high]],axis=0)
        wp=W/f'physics/cluster_{cl}_{cloud}_cool500.dat'
        wp.write_text('\n'.join(['EMBER_PRESSURE_MATCHED_ATMOSPHERE 1',
            f'scenario cluster_{cl}_{cloud}_bobcat_diamondback_sphinx_control',
            'pressure_cgs 1e8',f'composition {X:.17g} {Z:.17g} 0.002 0.00001',
            f'log_teff {len(finalx)} '+' '.join(f'{x:.17g}' for x in finalx),
            f'log_g {len(wg)} '+' '.join(f'{g:.17g}' for g in wg),'data']+
            [f'{a:.17g} {b:.17g}' for a,b in finalv.reshape(-1,2)])+'\n')
        receipt['cases'][cl+'_'+cloud].update(warm_file=str(wp),warm_sha256=hashlib.sha256(wp.read_bytes()).hexdigest())
        print(cl,cloud,len(out),'nodes',flush=True)
(W/'physics/colder_boundary_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
