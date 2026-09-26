#!/usr/bin/env python3
"""Independent spline, analytic thermodynamics and mask checks of the native EOS.

Synthetic potentials exercise the implementation; they do not accept the
physical FreeEOS family. SciPy constructs independent composition splines.
The thermodynamic oracle differentiates an explicit polynomial potential.
All generated source fixtures and probe responses are retained locally.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path
import subprocess
import numpy as np
from scipy.interpolate import CubicSpline
from composition_potential_v2 import (RGAS,mixing_value,mixing_derivatives,
                                      source_mixing_value)

ARAD=4*5.670374419e-5/2.99792458e10
X=np.array([0.,.0003,.001,.009,.07,.2,.45,.75])
Y=np.array([0.,.06,.12])
LT=np.array([5.6,5.8,6.,6.2,6.4])
LQ=np.array([-2.,-1.,0.,1.,2.])
T0=6*np.log(10.)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def coefficients(x,y):
    return np.array([.7+.5*x+.13*y+.025*np.exp(3*x)*(1+y*y),
                     1.2+.4*x+.02*np.sin(4*x)*(1+y),
                     .01*np.cos(8*x)*y*y,
                     .0002*(1+.5*np.cos(4*x))*(1+y*y)])


def tq_derivative(c,t,q,dt,dq):
    a,b,constant,d=c
    result=0.
    if dt==dq==0:result=a*(q+1.5*(t-T0))-b*t+constant
    elif (dt,dq)==(1,0):result=1.5*a-b
    elif (dt,dq)==(0,1):result=a
    if dt<=3 and dq<=2:
        result+=d*math.factorial(3)/math.factorial(3-dt)*math.factorial(2)/math.factorial(2-dq)*(t-T0)**(3-dt)*q**(2-dq)
    return RGAS*result


def physical_jet(c,T,rho):
    t=np.log(T);q=np.log(rho)-1.5*(t-T0)
    result=np.zeros((4,4))
    for i in range(4):
        for j in range(4-i):
            result[i,j]=sum(math.comb(i,k)*(-1.5)**k*tq_derivative(c,t,q,i-k,j+k) for k in range(i+1))
    return result


def fixture(work,masked=False):
    directory=work/('masked' if masked else 'complete');directory.mkdir()
    # Match the table's carried inert metal slots exactly; these slots are
    # merely bookkeeping for the fixed GS98 source approximation.
    source=Path('data/eos/numerical_electron_base_v1/freeeos300_gs98_x000_he3000_potential.dat')
    with source.open() as stream:
        for line in stream:
            if line.startswith('composition '):
                metals=list(map(float,line.split()[1:]))[3:];break
        else:raise ValueError('missing source composition')
    files=[]
    for ix,x in enumerate(X):
        for iy,y in enumerate(Y):
            name=f'x{ix}-y{iy}.dat';files.append(name)
            text=['EMBER_HELMHOLTZ 2','source "Synthetic potential for independent runtime checks"',
                  'composition_proxy "synthetic GS98 bookkeeping"','basis baryon_mass',
                  'metal_inventory gs98','composition '+' '.join(map(repr,[float(x),float(y),float(.98-x-y),*metals])),
                  'log_t '+str(len(LT))+' '+' '.join(map(str,LT)),
                  'log_q '+str(len(LQ))+' '+' '.join(map(str,LQ)),'data']
            c=coefficients(x,y)
            for lt in LT:
                for lq in LQ:
                    # An absent middle H plane leaves a three-node fragment
                    # at high density. The native slope support must reject it.
                    valid=not(masked and lq>=0 and ix==3)
                    node=[tq_derivative(c,lt*np.log(10),lq*np.log(10),i,j)
                          +(source_mixing_value(x,y) if i==j==0 else 0.)
                          for i in range(3) for j in range(3)]
                    text.append(str(int(valid))+' '+' '.join(format(v,'.17g') for v in node))
            (directory/name).write_text('\n'.join(text)+'\n')
    family=directory/'family.dat'
    family.write_text('EMBER_METAL_HELMHOLTZ 1\nhydrogen '+str(len(X))+' '+' '.join(map(str,X))+
                     '\nhelium3 '+str(len(Y))+' '+' '.join(map(str,Y))+'\n'+
                     '\n'.join(json.dumps(f) for f in files)+'\n')
    return family


class Oracle:
    def __init__(self,xaxis=X):
        self.xaxis=xaxis
        self.spline=CubicSpline(xaxis,np.array([[coefficients(x,y) for y in Y] for x in xaxis]),axis=0)
        self.by=CubicSpline(Y,np.eye(3),axis=0)

    def c(self,x,y,dx=0,dy=0):return self.by(y,dy)@self.spline(x,dx)

    def state(self,x,y,T,rho):
        f=physical_jet(self.c(x,y),T,rho);f[0,0]+=mixing_value(x,y)
        rad=-ARAD*T**3/(3*rho)
        for i in range(4):
            for j in range(4-i):f[i,j]+=rad*3**i*(-1)**j
        P=rho*T*f[0,1];E=-T*f[1,0];S=-f[0,0]-f[1,0]
        cv=-f[1,0]-f[2,0];ct=1+f[1,1]/f[0,1];cr=1+f[0,2]/f[0,1]
        delta=ct/cr;cp=cv+f[0,1]*ct*delta;ad=f[0,1]*delta/cp
        return np.array([P,E,S,cv,cp,ct,cr,delta,ad,cr*cp/cv])

    def expected(self,x,y,T,rho):
        s=self.state(x,y,T,rho);f=physical_jet(self.c(x,y),T,rho)
        # Independently difference the analytic oracle for third-derivative
        # responses. Two spacings verify that this measurement is resolved.
        derivatives=[]
        for h in [1e-5,5e-6]:
            t=(self.state(x,y,T*np.exp(h),rho)-self.state(x,y,T*np.exp(-h),rho))/(2*h)
            r=(self.state(x,y,T,rho*np.exp(h))-self.state(x,y,T,rho*np.exp(-h)))/(2*h)
            derivatives.append(np.array([r[1],t[4],r[4],t[7],r[7],t[8],r[8]]))
        source_spacing=np.max(np.abs(derivatives[0]-derivatives[1])/np.array([max(abs(s[1]),1),s[4],s[4],s[7],s[7],s[8],s[8]]))
        g=[physical_jet(self.c(x,y,int(k==0),int(k==1)),T,rho) for k in range(2)]
        pe=[rho*T*g[0][0,1],rho*T*g[1][0,1],-T*g[0][1,0],-T*g[1][1,0]]
        grad,hess=mixing_derivatives(x,y)
        for k in range(2):
            grad[k]+=g[k][0,0]
            for l in range(2):
                hess[k,l]+=physical_jet(self.c(x,y,int(k==0)+int(l==0),int(k==1)+int(l==1)),T,rho)[0,0]
        chemical=[f[0,0]+mixing_value(x,y),*grad,*hess.flatten(),g[0][1,0],g[1][1,0],g[0][0,1],g[1][0,1]]
        return np.array([*s,*derivatives[1],*pe,*chemical]),source_spacing


def query(probe,family,rows,work,label,mode='smooth'):
    inp=work/(label+'.input');out=work/(label+'.output');err=work/(label+'.stderr')
    inp.write_text(''.join(' '.join(map(str,row))+'\n' for row in rows))
    with inp.open('rb') as i,out.open('xb') as o,err.open('xb') as e:
        p=subprocess.run([str(probe),str(family),mode],stdin=i,stdout=o,stderr=e,timeout=180)
    if p.returncode:raise RuntimeError(f'{label} probe failed: '+err.read_text())
    result=[json.loads(line) for line in out.read_text().splitlines()]
    if len(result)!=len(rows):raise ValueError('missing native replies')
    return result


def main():
    global ARAD
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('probe',type=Path);p.add_argument('work',type=Path);p.add_argument('report',type=Path);args=p.parse_args()
    import re
    constants=Path('include/ember/constants.hpp').read_text()
    def constant(name):return float(re.search(r'\b'+name+r'\s*=\s*([\d.eE+-]+)',constants).group(1))
    ARAD=4*constant('sigma_SB')/constant('c')
    args.work.mkdir()
    family=fixture(args.work);masked=fixture(args.work,True);oracle=Oracle()
    rows=[]
    for x,y in [(1e-5,.002),(.0006,.03),(.005,.09),(.12,.06),(.36,.11),(.7,.04)]:
        for lt,lq in [(5.73,-1.4),(6.13,-.4),(6.27,1.3)]:
            T=10**lt;rho=10**lq*(T/1e6)**1.5;rows.append([x,y,T,rho,1])
    for x in X[1:-1]:
        for dx in [-1e-10,0,1e-10]:rows.append([float(x+dx),.045,10**6.07,3.7,1])
    actual=query(args.probe,family,rows,args.work,'analytic')
    comparisons=[]
    for row,response in zip(rows,actual):
        if not response['ok']:raise ValueError(response)
        expected,spacing=oracle.expected(*row[:4]);v=np.array(response['values'][:32])
        scale=np.maximum(np.abs(expected),1.)
        errors=np.abs(v-expected)/scale
        comparisons.append(dict(query=row,maximum_relative_error=float(max(errors)),oracle_spacing=float(spacing),errors=errors.tolist(),passed=bool(max(errors)<2e-6 and spacing<1e-7)))
    # Endpoint values remain finite. Chemical derivatives of an absent species
    # are rejected explicitly, as are composition and material extrapolation.
    guards=[[0,.06,1e6,1,0],[.1,0,1e6,1,0],[.75,.12,1e6,1,1],
            [0,.06,1e6,1,1],[.1,0,1e6,1,1],[-1e-5,.06,1e6,1,0],
            [.751,.06,1e6,1,0],[.1,.121,1e6,1,0],[.1,.06,1e5,1,0],
            [.1,.06,1e6,101,0],[.1,.06,-1,1,0]]
    g=query(args.probe,family,guards,args.work,'guards')
    guard_pass=[a['ok']==(i<3) for i,a in enumerate(g)]
    masks=[[.0006,.03,1e6,.03,1],[.0006,.03,1e6,1,1],[.36,.03,1e6,1,1],
           [.005,.03,1e6,1,1]]
    mr=query(args.probe,masked,masks,args.work,'masks')
    mask_pass=[a['ok']==(i in [0,2]) for i,a in enumerate(mr)]
    # The surviving four-node segment uses its own end conditions, without
    # bringing any value from the masked source plane across the gap.
    local_expected,_=Oracle(X[4:]).expected(*masks[2][:4])
    local_error=float(np.max(np.abs(np.array(mr[2]['values'][:32])-local_expected)/np.maximum(np.abs(local_expected),1))) if mr[2]['ok'] else None
    # Density inversion's advertised support must exclude the same fragment.
    range_pass=mr[0]['ok'] and mr[0]['values'][-1]<.101 and mr[2]['ok'] and mr[2]['values'][-1]>99
    accepted=all(x['passed'] for x in comparisons) and all(guard_pass) and all(mask_pass) and range_pass and local_error is not None and local_error<2e-6
    report=dict(created_utc=datetime.now(timezone.utc).isoformat(),scope=__doc__,outcome='passed' if accepted else 'failed',accepted_for_evolution=False,
                criteria=dict(relative_oracle=2e-6,oracle_spacing=1e-7),comparisons=comparisons,
                guard_pass=guard_pass,mask_pass=mask_pass,masked_segment_relative_error=local_error,density_range_pass=bool(range_pass),
                probe_sha256=sha(args.probe),inputs_sha256={str(p):sha(p) for p in [Path(__file__),Path('scripts/composition_potential_v2.py'),Path('src/eos_smooth_mixture.cpp'),Path('src/eos_helmholtz.cpp'),Path('include/ember/eos_smooth_mixture.hpp')]},
                artifacts={str(p):sha(p) for p in args.work.rglob('*') if p.is_file()},remaining_work=['Physical source-family assembly and source comparisons','Full stellar profile and atmosphere checks','Runtime cost and memory measurement','Transport integration and evolution'])
    args.report.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ['outcome','guard_pass','mask_pass','masked_segment_relative_error','density_range_pass']}))
    if not accepted:raise SystemExit(1)


if __name__=='__main__':main()
