"""Collisionless association expansion in the same Galactic potential as PeTar.

Positions are pc, velocities pc/Myr; galpy is evaluated in ro=8 kpc, vo=220 km/s.
No local linearization of Galactic forces is used.
"""
from pathlib import Path
import json,time
import numpy as np
from scipy.integrate import solve_ivp
from galpy.potential import MWPotential2014,evaluateRforces,evaluatezforces

W=Path(__file__).resolve().parent
KMS=1.022712165045695;RO=8000.;VO=220*KMS;TIME=RO/VO

def evolve(start,age,rtol=2e-9):
    start=np.asarray(start,float);n=len(start)
    scale=np.array([RO]*3+[VO]*3)
    def rhs(t,flat):
        y=flat.reshape(n,6);p=y[:,:3];r=np.hypot(p[:,0],p[:,1])
        fr=evaluateRforces(MWPotential2014,r,p[:,2],use_physical=False)
        fz=evaluatezforces(MWPotential2014,r,p[:,2],use_physical=False)
        return np.column_stack([y[:,3:],fr*p[:,0]/r,fr*p[:,1]/r,fz]).ravel()
    s=time.monotonic()
    sol=solve_ivp(rhs,[0,age/TIME],(start/scale).ravel(),method='DOP853',
                  rtol=rtol,atol=rtol*.02,t_eval=[age/TIME])
    if not sol.success:raise RuntimeError(sol.message)
    return sol.y[:,-1].reshape(n,6)*scale,dict(wall_seconds=time.monotonic()-s,n_orbits=n,rhs_calls=sol.nfev,rtol=rtol)

if __name__=='__main__':
    m=json.loads((W.parent/'cluster-dissolution-oct9-v2/runs/hyades_r2.5_assort_s41/initial.json').read_text())
    rng=np.random.default_rng(805)
    initial=np.array(m['orbit']['start'])+rng.normal(size=(1024,6))*[30,30,30,KMS*.7,KMS*.7,KMS*.7]
    initial[0]=m['orbit']['start']
    final,receipt=evolve(initial,m['age_Myr'])
    check,precise=evolve(initial[:32],m['age_Myr'],rtol=2e-11)
    error=np.linalg.norm(check[:,:3]-final[:32,:3],axis=1)
    center_error=np.linalg.norm(final[0,:3]-np.array(m['orbit']['present'])[:3])
    assert center_error<.001 and error.max()<.01,(center_error,error.max())
    receipt.update(center_position_error_pc=center_error,max_position_error_tighter_control_pc=float(error.max()),control=precise)
    (W/'galactic_integration_check.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
