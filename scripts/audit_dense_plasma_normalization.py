#!/usr/bin/env python3
"""Bound the continuous plasma normalization effect at fixed dense table states.

Complete ionization supplies the electron-density upper bound. Source frequency
bin rounding, frequency-dependent refraction and evolution are excluded.
"""
import argparse
from pathlib import Path
import json,gzip,hashlib,math
import numpy as np
from scipy.integrate import quad
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('transport_report',type=Path);parser.add_argument('output',type=Path)
args=parser.parse_args()
source=args.transport_report;meta=json.loads(source.read_text());raw=Path(meta['raw_queries']);assert hashlib.sha256(raw.read_bytes()).hexdigest()==meta['raw_queries_sha256']
records=json.loads(gzip.decompress(raw.read_bytes()))['records']
points=np.array([[r[k] for k in ['X','Z','Y3','temperature_K','density','radiative_interpolated','conductive_opacity']] for r in records])
x,z,y3,T,rho,kr,kc=points.T;ye=x+2*y3/3+(1-x-z-y3)/2+z/2
u=1.054571817e-27*np.sqrt(4*np.pi*4.80320471257e-10**2*rho*ye/1.66053906660e-24/9.1093837015e-28)/(1.380649e-16*T)
def removed(cutoff,n):
 nodes,weights=np.polynomial.legendre.leggauss(n);answer=np.empty(len(cutoff))
 for i in range(0,len(cutoff),1024):
  b=np.minimum(cutoff[i:i+1024],40);v=b[:,None]*(nodes+1)/2
  integrand=v**4*np.exp(-v)/np.expm1(-v)**2
  answer[i:i+1024]=b/2*np.dot(integrand,weights)/(4*np.pi**4/15)
 return answer
r64=removed(u,64);r128=removed(u,128);error=float(np.max(abs(r128-r64)));assert error<1e-12
fraction=kc/(kc+kr);increase=1/(1-r128*fraction)-1
result=dict(scope='Continuous classical plasma-normalization diagnostic at dense opacity source/interpolation states. Fully ionized H/He and Z/A <= 1/2 give an electron-density upper bound. This excludes source-bin rounding and frequency-dependent refractive index and is not an evolutionary error bound.',input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,raw,Path(__file__)]},states=len(records),quadrature_64_128_maximum_absolute_difference=error,groups={})
for name,mask in [('all',np.ones(len(records),dtype=bool)),('one_to_twenty_MK',(T>=1e6)&(T<=2e7)),('hydrogen_poor_z020',(T>=1e6)&(T<=2e7)&(z==.02)&(x<=.2))]:
 ids=np.flatnonzero(mask);i=int(ids[np.argmax(increase[ids])]);result['groups'][name]=dict(states=int(mask.sum()),maximum_combined_opacity_relative_increase=float(increase[i]),classical_cutoff_u=float(u[i]),removed_rosseland_weight=float(r128[i]),state=records[i])
args.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v['maximum_combined_opacity_relative_increase'] for k,v in result['groups'].items()}))
