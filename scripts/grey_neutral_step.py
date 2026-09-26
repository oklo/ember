"""Conservative mixing at a marginal convective face in the grey CN model.

The usual finite step is retained when no partial exchange is needed.
A nonnegative common species flux may hold one interface at Schwarzschild
neutrality, the stiff limit of the same mixing-length prescription. The
caller must independently check timestep accuracy and conservation. This
closure is selected only for the grey model without microscopic diffusion.
"""
from scipy.optimize import brentq
from evolve_cn_transport import audit
from audit_cn_thermal import YEAR

def neutral_step(m,years,previous,request,tolerance=1e-12):
 cache={};trace=[]
 def call(lam,face):
  key=(face,lam)
  if key in cache:return cache[key]
  r=request(face,lam,m,years,tolerance)
  a,c,b=audit(m,r,years*YEAR);cache[key]=(a,r,c,b);trace.append(dict(parameter=lam,face=face,stability=r.get('boundary_stability'),cpu_seconds=r['cpu_seconds']))
  return cache[key]
 # A prior neutral interface is the best nonlinear-solve starting guess.
 # Test its unmixed and strongly mixed limits before any global retry.
 boundaries=[]
 for lo,hi in previous['mixing_regions']:
  if hi-lo<=1:continue
  if lo>0:boundaries.append(lo-1)
  if hi<len(m):boundaries.append(hi-1)
 boundaries=sorted(set(boundaries),reverse=True)
 previous_face=previous.get('neutral_face')
 if previous.get('mixing_parameter',0)>0 and previous_face in boundaries:
  boundaries.remove(previous_face);boundaries.insert(0,previous_face)
 faces=[]
 for offset in [0,-1,1,-2,2,-3,3,-4,4]:
  for face in boundaries:
   candidate=face+offset
   if 0<=candidate<len(m)-1 and candidate not in faces:faces.append(candidate)
 failures=[]
 for face in faces:
  try:
   first=call(0.,face)
   if first[1]['boundary_stability']<=0:
    return first,trace
   def f(lam):return call(lam,face)[1]['boundary_stability']
   hi=1.;upper=f(hi)
   if upper>0:hi=100.;upper=f(hi)
   if upper>0:hi=1e4;upper=f(hi)
   if upper>=0:continue
   root=brentq(f,0.,hi,xtol=1e-10,rtol=1e-10)
   result=call(root,face);rec=result[1]
   assert abs(rec['boundary_stability'])<1e-9
   assert rec['mixing_rate']>=0
   assert rec['mlt_required_radiative_excess']/rec['boundary_grad_ad']<1e-8
   assert any((a==face+1 or b==face+1) and b-a>1 for a,b in rec['mixing_regions'])
   return result,trace
  except Exception as e:failures.append(dict(face=face,error=repr(e)))
 try:return call(0.,-1),trace
 except Exception as e:raise RuntimeError('No conservative neutral interface found: '+repr(failures)+'; ordinary: '+repr(e))
