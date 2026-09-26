#!/usr/bin/env python3
"""Independent source, export and native-domain review of the ion extension."""
from pathlib import Path
import hashlib,json,re,subprocess,time
from datetime import datetime,timezone
import numpy as np
from scipy.interpolate import PPoly
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expit
from yukawa_scattering_control_v2 import collision_integrals


def main():
    start=time.monotonic();sealed=Path('docs/research/fable/results/fable-ion-table-extension-v1')
    output=Path('docs/results/fable_ion_extension_review_v1.json')
    if output.exists():raise ValueError('preserve completed review')
    hashes={}
    def pin(p):
        p=Path(p);h=hashlib.sha256(p.read_bytes()).hexdigest();hashes[str(p)]=h;return h
    ready=json.loads((sealed/'READY.json').read_text())
    for name,h in ready['files_sha256'].items():assert pin(sealed/name)==h
    for name,h in ready['input_sha256'].items():assert pin(name)==h
    old=Path('/tmp/ember-collision-transport-table-v1.dat')
    new=sealed/'collision_transport_table_candidate_v2.dat'
    before=old.read_text().splitlines();after=new.read_text().splitlines()
    assert before[:52]==after[:52] and before[55:]==after[55:]
    def polynomial(lines,at):
        n,order=map(int,lines[at].split());x=np.fromstring(lines[at+1],sep=' ')
        c=np.fromstring(lines[at+2],sep=' ').reshape(order,n-1,4)
        assert len(x)==n;return PPoly(c,x,extrapolate=False)
    a,b=polynomial(before,52),polynomial(after,52)
    assert np.array_equal(a.x,b.x[:len(a.x)]) and np.array_equal(a.c,b.c[:,:a.c.shape[1]])
    join=a.x[-1];at=a.c.shape[1]
    following=PPoly(b.c[:,at:at+1],b.x[at:at+2],extrapolate=False)
    continuity=[float(abs((a(join,d)-following(join,d))/(1+abs(a(join,d)))).max()) for d in range(3)]
    assert max(continuity)<1e-10
    controls=[]
    for logg in [2.4832,2.9011]:
        low=np.array(list(collision_integrals(10**logg).values()))
        fine=np.array(list(collision_integrals(10**logg,thermal_order=48,angle_order=160,radius_order=320,thermal_max=10).values()))
        actual=np.exp(b(logg));error=float(abs(actual/fine-1).max());quadrature=float(abs(low/fine-1).max())
        assert error<1e-3 and quadrature<2e-5
        controls.append(dict(log10_strength=logg,direct=fine.tolist(),exported=actual.tolist(),interpolation_error=error,quadrature_change=quadrature))
    header=Path('include/ember/gs98_mixture.hpp');pin(header)
    metals=np.array([list(map(float,s.split(','))) for s in re.findall(r'^\s*\{([^{}]+)\}, //',header.read_text(),re.M)])
    moment=np.sum(metals[:,0]*metals[:,3]/metals[:,1])
    kb,amu,hp,me=1.380649e-16,1.66053906660e-24,6.62607015e-27,9.1093837015e-28
    queries=[];labels=[]
    for file,indices in [('out/evolution-metal-512-2500gyr-gas-checkpoint.json',[0,128,263,300,340,379]),('out/evolution-cold-remnant-forward-512-3400gyr-v1.json',[0,240,316,340,380,394])]:
        pin(file);m=np.array(json.loads(Path(file).read_text())['profile'])
        for i in indices:
            T,rho,X,Y3=m[i,3],m[i,2],m[i,5],m[i,6];Y4=.98-X-Y3
            ne=rho/amu*(X+2*Y3/3+Y4/2+.02*moment)
            prefactor=4*np.pi*(2*me*kb*T)**1.5/hp**3
            integral=lambda eta:quad(lambda x:np.sqrt(x)*expit(eta-x),0,max(60,eta+40),epsabs=1e-11,epsrel=1e-11)[0]
            eta=brentq(lambda eta:prefactor*integral(eta)/ne-1,-40,100)
            derivative=quad(lambda x:np.sqrt(x)*expit(eta-x)*expit(x-eta),0,max(60,eta+40),epsabs=1e-11,epsrel=1e-11)[0]
            stiffness=kb*T*integral(eta)/derivative
            queries.append(' '.join(format(x,'.17g') for x in [T,rho,X,Y3,.02,0,stiffness,1]));labels.append([file,i])
    native=Path('/tmp/ember-collision-derivative-build-v1/tests/collision_transport_probe');pin(native)
    replies=[]
    for path in [old,new]:
        p=subprocess.run([str(native),str(path),'derivatives'],input='\n'.join(queries)+'\n',text=True,capture_output=True,check=True)
        replies.append([json.loads(line) for line in p.stdout.splitlines()]);assert len(replies[-1])==len(queries)
    states=[]
    for label,q,previous,current in zip(labels,queries,*replies):
        assert 'error' not in current
        if 'error' not in previous:
            previous.pop('seconds');current.pop('seconds');assert previous==current
        states.append(dict(profile=label[0],node=label[1],query=q,previously_supported='error' not in previous,response=current))
    # An explicit beyond-domain screening length still rejects.
    T,rho,X,Y3,*_=map(float,queries[0].split())
    # Keep electron degeneracy fixed and its dimensionless screening length
    # inside the table, so this rejects specifically on the ion coordinate.
    # This is a numerical-domain test of the prescribed ionized mixture.
    rho *= (3e5/T)**1.5;T=3e5
    length=784*2.307077552341736e-19/(kb*T*10**3.)
    q=' '.join(format(v,'.17g') for v in [T,rho,X,Y3,.02,length])+'\n'
    p=subprocess.run([str(native),str(new)],input=q,text=True,capture_output=True,check=True)
    assert 'ion collision query outside table' in json.loads(p.stdout)['error']
    for name in [__file__,'scripts/yukawa_scattering_control_v2.py','scripts/yukawa_scattering_control.py']:pin(name)
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome='passed_independent_source_export_and_native_controls',seconds=time.monotonic()-start,
                old_coefficient_blocks_preserved=True,join_errors_through_second_derivative=continuity,new_source_controls=controls,native_states=states,
                unsupported_still_rejected=True,sha256=hashes,scope='Same classical repulsive Yukawa and fully stripped mixture prescription; extends numerical coverage only. Stellar heat response remains a separate comparison.')
    output.write_text(json.dumps(report,indent=2)+'\n');print(report['outcome'],len(states),'native states',controls,flush=True)


if __name__=='__main__':main()
