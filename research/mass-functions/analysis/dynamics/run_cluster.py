"""Private, reproducible PeTar cluster pilot. Units: Msun, pc, Myr.

Galactic orbit uses galpy MWPotential2014, ro=8 kpc, vo=220 km/s.
Stellar mass loss uses PeTar SSE/BSE; low-mass mock light uses Ember later.
No gravitational softening and no escaper removal. Initially no mass segregation.
"""
from pathlib import Path
import argparse, json, os, subprocess, time, resource
import numpy as np
from scipy.spatial.distance import pdist, squareform
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014
from astropy import units as u
from astropy.coordinates import SkyCoord, CartesianRepresentation, CartesianDifferential

W=Path(__file__).resolve().parent
G=.00449830997959438
KMS=1.022712165045695
CLUSTERS={'hyades':dict(ra=66.75,dec=15.87,distance_pc=47.03,age_Myr=650.,Z=.02307),
          'pleiades':dict(ra=56.75,dec=24.12,distance_pc=136.2,age_Myr=125.,Z=.01648)}

def integral(a,b,alpha):
    return (b**(1-alpha)-a**(1-alpha))/(1-alpha)

def masses(n,rng):
    # R(.03-.08)/R(.1-.2)=.63; above .08 use Kroupa slopes.
    bdamp=.63*integral(.1,.2,1.3)/integral(.03,.08,.4)
    seg=[(.01,.08,.4,bdamp),(.08,.5,1.3,1.),(.5,50.,2.3,.5)]
    prob=np.array([c*integral(a,b,s) for a,b,s,c in seg]);prob/=prob.sum()
    k=rng.choice(3,n,p=prob);m=np.empty(n)
    for i,(a,b,s,c) in enumerate(seg):
        z=rng.uniform(size=(k==i).sum());m[k==i]=(a**(1-s)+z*(b**(1-s)-a**(1-s)))**(1/(1-s))
    return m

def direction(n,rng):
    a=rng.normal(size=(n,3));return a/np.linalg.norm(a,axis=1)[:,None]

def plummer(m,rh,rng):
    n=len(m);a=rh/1.304766
    # Mild outer truncation, r<20 rh, then rescale actual half-mass radius.
    umax=(1+(a/(20*rh))**2)**(-1.5)
    r=a/np.sqrt(rng.uniform(0,umax,n)**(-2/3)-1)
    p=r[:,None]*direction(n,rng);p-=np.average(p,weights=m,axis=0)
    r=np.linalg.norm(p,axis=1);order=np.argsort(r)
    rhalf=np.interp(m.sum()/2,np.cumsum(m[order]),r[order]);p*=rh/rhalf
    q=np.empty(n);todo=np.arange(n)
    while len(todo):
        x=rng.uniform(size=len(todo));y=rng.uniform(0,.1,len(todo));ok=y<x*x*(1-x*x)**3.5
        q[todo[ok]]=x[ok];todo=todo[~ok]
    r=np.linalg.norm(p,axis=1)
    v=(q*np.sqrt(2*G*m.sum()/np.sqrt(r*r+a*a)))[:,None]*direction(n,rng)
    v-=np.average(v,weights=m,axis=0)
    # Exact virial normalization of cluster CM motion, not binary internal energy.
    d=squareform(pdist(p));np.fill_diagonal(d,np.inf)
    pot=-.5*G*np.sum(m[:,None]*m[None,:]/d)
    kin=.5*np.sum(m[:,None]*v*v);v*=np.sqrt(-.5*pot/kin)
    return p,v,dict(potential=pot,kinetic=-.5*pot,rhalf_pc=rh)

def cluster_orbit(cl):
    c=CLUSTERS[cl]
    if cl=='hyades':
        # Gaia Collaboration 2018/Gaia Hyades study: equatorial velocity.
        xyz=SkyCoord(c['ra']*u.deg,c['dec']*u.deg,distance=c['distance_pc']*u.pc).cartesian
        sc=SkyCoord(xyz.with_differentials(CartesianDifferential(np.array([-5.96,45.60,5.57])*u.km/u.s)),frame='icrs')
        obs=[c['ra'],c['dec'],c['distance_pc']/1000,sc.pm_ra_cosdec.to_value(u.mas/u.yr),sc.pm_dec.to_value(u.mas/u.yr),sc.radial_velocity.to_value(u.km/u.s)]
    else:
        obs=[c['ra'],c['dec'],c['distance_pc']/1000,19.997,-45.548,5.79]
    o=Orbit(obs,radec=True,ro=8.,vo=220.,zo=.0208,solarmotion='schoenrich')
    times=np.linspace(0,-c['age_Myr'],2001)*u.Myr
    o.integrate(times,MWPotential2014,method='dop853')
    tt=times[-1]
    start=[o.x(tt)*1000,o.y(tt)*1000,o.z(tt)*1000,o.vx(tt)*KMS,o.vy(tt)*KMS,o.vz(tt)*KMS]
    final=[o.x(0)*1000,o.y(0)*1000,o.z(0)*1000,o.vx(0)*KMS,o.vy(0)*KMS,o.vz(0)*KMS]
    return dict(observables=obs,start=start,present=final,ro_kpc=8,vo_kms=220,zo_kpc=.0208,solarmotion='schoenrich')

def prepare(args):
    path=W/'runs'/args.name;path.mkdir(parents=True,exist_ok=False)
    rng=np.random.default_rng(args.seed);m=masses(args.n,rng)
    # First paired components are random: binary control, not a fitted population.
    nb=int(args.binary_fraction*args.n/(1+args.binary_fraction))
    if args.binary_mode=='assortative':
        # Preserve every sampled component mass, but suppress star-BD pairs.
        # BD f_system=.15 with q~.7-1; stellar f_system=.30, q~.5-1 below
        # .5 Msun and .1-1 above. This is a sensitivity model, not a fit.
        pairs=[];remaining=[]
        for select,fbin,isbd in [(m<.08,.15,True),(m>=.08,.30,False)]:
            available=list(np.flatnonzero(select));npair=int(fbin*len(available)/(1+fbin))
            for _ in range(npair):
                first=int(rng.choice(available));available.remove(first)
                q=rng.uniform(.7 if isbd else (.5 if m[first]<.5 else .1),1.)
                partner=min(available,key=lambda j:abs(np.log(m[j]/(q*m[first]))))
                available.remove(partner);pairs.extend([first,partner])
            remaining+=available
        m=m[np.r_[pairs,remaining].astype(int)];nb=len(pairs)//2
    sm=np.r_[m[:2*nb].reshape(-1,2).sum(axis=1),m[2*nb:]]
    p,v,vir=plummer(sm,args.rh,rng)
    if nb:
        upper=np.where(np.max(m[:2*nb].reshape(-1,2),axis=1)<.08,20.,200.) if args.binary_mode=='assortative' else np.full(nb,200.)
        semis=10**rng.uniform(np.log10(.2),np.log10(upper),nb)/206264.806
        dr=direction(nb,rng)*semis[:,None]
        dv=np.cross(dr,direction(nb,rng));dv/=np.linalg.norm(dv,axis=1)[:,None]
        dv*=np.sqrt(G*sm[:nb]/semis)[:,None]
        f=m[:2*nb].reshape(-1,2)/sm[:nb,None]
        bp=np.stack([p[:nb]+f[:,1,None]*dr,p[:nb]-f[:,0,None]*dr],axis=1).reshape(-1,3)
        bv=np.stack([v[:nb]+f[:,1,None]*dv,v[:nb]-f[:,0,None]*dv],axis=1).reshape(-1,3)
        p=np.concatenate([bp,p[nb:]]);v=np.concatenate([bv,v[nb:]])
    orbit=cluster_orbit(args.cluster)
    np.savetxt(path/'particles.dat',np.column_stack([m,p,v]),fmt='%.16e')
    np.savez_compressed(path/'initial.npz',mass=m,pos=p,vel=v)
    meta=dict(vars(args));meta.update(CLUSTERS[args.cluster]);meta.update(orbit=orbit,virial=vir,actual_binary_fraction_by_system=nb/(len(m)-nb),nbinary=nb,total_mass=m.sum(),birth_ratio=float(((m>=.03)&(m<.08)).sum()/((m>=.1)&(m<.2)).sum()))
    (path/'initial.json').write_text(json.dumps(meta,indent=2)+'\n')
    cmd=[str(W/'install/bin/petar.init'),'-c',','.join(str(x) for x in orbit['start']),'-t','-s','bse','-v','1','-f','input','particles.dat']
    subprocess.run(cmd,cwd=path,stdout=(path/'init.log').open('w'),stderr=subprocess.STDOUT,check=True)
    assert sum(1 for _ in (path/"input").open()) == len(m)+1, "Incomplete PeTar input"
    return path,meta

def run(args):
    path,meta=prepare(args)
    env=os.environ.copy();env.update(OMP_STACKSIZE='128M',OMP_NUM_THREADS=str(args.threads))
    cmd=[str(W/'install/bin/petar'),'-u','1','-b',str(meta['nbinary']), '--galpy-set','MWPotential2014','--bse-metallicity',str(meta['Z']),'-t',str(args.end or meta['age_Myr']),'-o',str(args.output),'-f','data','--stellar-evolution',str(args.stellar_evolution),'-T',str(args.theta),'input']
    if args.step: cmd[-1:-1]=['-s',str(args.step)]
    cmd[-1:-1]=['--number-group-limit',str(args.group)]
    (path/'command.json').write_text(json.dumps(dict(command=cmd,threads=args.threads),indent=2))
    start=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
    with (path/'run.log').open('w') as f:
        p=subprocess.Popen(cmd,cwd=path,env=env,stdout=f,stderr=subprocess.STDOUT)
        (path/'pid').write_text(str(p.pid));code=p.wait()
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    receipt=dict(exit_code=code,wall_seconds=time.monotonic()-start,cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime)
    (path/'receipt.json').write_text(json.dumps(receipt,indent=2));print(args.name,receipt,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('name');p.add_argument('--cluster',choices=CLUSTERS,default='hyades');p.add_argument('--n',type=int,default=3000);p.add_argument('--rh',type=float,default=2);p.add_argument('--seed',type=int,default=41);p.add_argument('--end',type=float,default=0);p.add_argument('--threads',type=int,default=4);p.add_argument('--output',type=float,default=25);p.add_argument('--binary-fraction',type=float,default=0)
    p.add_argument('--step',type=float,default=0);p.add_argument('--theta',type=float,default=.3);p.add_argument('--stellar-evolution',type=int,default=1)
    p.add_argument('--group',type=int,default=512)
    p.add_argument('--binary-mode',choices=['random','assortative'],default='random')
    run(p.parse_args())
