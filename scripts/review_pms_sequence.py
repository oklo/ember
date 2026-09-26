"""Audit actual saved PMS structures and isotope changes without running evolution.

Accepts either a birth model or an explicitly pinned parent checkpoint. This
review covers homogeneous, deuterium-dominated intervals only. It evaluates
independent source powers, implicit pp contributions and storage precision;
it neither changes acceptance limits nor selects a physical atmosphere.
"""
from pathlib import Path
import sys,json,math,hashlib,subprocess,datetime,resource
import numpy as np
import argparse
parser=argparse.ArgumentParser(description='Independently audit saved homogeneous, deuterium-dominated PMS intervals, including pp isotope production.')
parser.add_argument('run',type=Path)
parser.add_argument('work',type=Path)
parser.add_argument('output',type=Path)
parser.add_argument('--probe',type=Path,required=True,help='Compiled source probe emitting total heat, neutrinos, D heat and nine pp isotope rates.')
parser.add_argument('--probe-source',type=Path,required=True)
parser.add_argument('--parent-checkpoint',type=Path)
args=parser.parse_args()
R=Path(__file__).resolve().parents[1];P=args.run.resolve();W=args.work.resolve()
W.mkdir(exist_ok=False)
sys.path.insert(0,str(R/'scripts'))
from review_pms_initial_contraction import checkpoint,BINDING,D_HEAT,YEAR,LSUN
probe=args.probe.resolve();probe_source=args.probe_source.resolve()
plan=json.loads((P/'plan.json').read_text());receipt=json.loads((P/'result.json').read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(P/'plan.json')==receipt['plan_sha256']
for f,h in plan['input_sha256'].items():assert sha(f)==h,f
attempts=[json.loads(s)for s in (P/'attempts.jsonl').read_text().splitlines()]
accepted=[a for a in attempts if a['accepted']]
assert len(accepted)==receipt['native_result']['accepted_intervals']
models={};files=[]
for i in range(len(accepted)):
 for name in ['seed','full','half1','final']:
  f=P/f'interval-{i}'/(name+'.checkpoint');models[i,name]=checkpoint(f);files.append(f)
assert np.array_equal(models[0,'seed']['data'],checkpoint(P/'seed.checkpoint')['data'])
assert np.array_equal(models[len(accepted)-1,'final']['data'],checkpoint(P/'final.checkpoint')['data'])
if args.parent_checkpoint:
 parent=args.parent_checkpoint.resolve()
 restart_flags=[f for f in ['--resume','--refine-abundances-from'] if f in plan['command']]
 assert len(restart_flags)==1
 flag=restart_flags[0]
 assert Path(plan['command'][plan['command'].index(flag)+1]).resolve()==parent
 if flag=='--refine-abundances-from':
  assert receipt['native_result']['explicit_abundance_refinement']==1
  assert receipt['native_result']['abundance_tolerance']==1e-15
  assert '--source-executable' in plan['command']
  original=plan['command'][plan['command'].index('--source-executable')+1]
  assert original==receipt['native_result']['source_executable']
  assert sha(original)==plan['input_sha256'][original]
 previous=checkpoint(parent)
 assert np.array_equal(models[0,'seed']['data'],previous['data']) and models[0,'seed']['age']==previous['age']
 files.append(parent)
else:
 assert '--resume' not in plan['command'] and '--refine-abundances-from' not in plan['command']
 birth=np.loadtxt(plan['command'][10]);assert np.all(models[0,'seed']['data'][:,7:]==birth)
 assert models[0,'seed']['age']==0

for i in range(1,len(accepted)):
 assert np.array_equal(models[i,'seed']['data'],models[i-1,'final']['data'])
 assert models[i,'seed']['age']==models[i-1,'final']['age']
keys=[(i,n)for i in range(len(accepted))for n in ['full','half1','final']]
def limit():resource.setrlimit(resource.RLIMIT_CPU,(120,125))
with (W/'nuclear.txt').open('x')as output,(W/'nuclear.stderr').open('x')as err:
 p=subprocess.Popen([str(probe)],stdin=subprocess.PIPE,stdout=output,stderr=err,text=True,preexec_fn=limit)
 for k in keys:
  m=models[k];a=m['data'];p.stdin.write(f"{len(a)} {m['mass']:.17g}\n")
  table=np.c_[a[:,0],np.exp(a[:,3]),np.exp(a[:,2]),a[:,7:]]
  np.savetxt(p.stdin,table,fmt='%.17g')
 p.stdin.close();assert p.wait(timeout=180)==0
source_values=np.loadtxt(W/'nuclear.txt')
assert source_values.shape==(len(keys),12) and np.isfinite(source_values).all()
powers=dict(zip(keys,source_values))
rows=[]
for i,a in enumerate(accepted):
 assert a['converged'] and a['audit_pass'] and max(abs(v)for v in a['first_law_errors'])<2e-7
 seed=models[i,'seed'];mass=seed['mass'];m=seed['data'][:,0]
 w=np.diff(np.r_[0,(m[:-1]+m[1:])/2,mass])
 dt=a['dt_years']*YEAR;step=[]
 for start,end,ds in [('seed','full',dt),('seed','half1',dt/2),('half1','final',dt/2)]:
  before=models[i,start];after=models[i,end];x0=before['data'][:,7:];x1=after['data'][:,7:];dx=x1-x0
  assert np.array_equal(before['data'][:,0],after['data'][:,0])
  assert np.all(x1==x1[0]) and np.all(x1[:,3:8]==x0[:,3:8])
  assert abs(after['age']-before['age']-ds)<=8*math.ulp(after['age'])
  nuc,nu,D=powers[i,end][:3];release=float(nuc+nu);fuel=math.fsum(float(wi)*float(-d)for wi,d in zip(w,dx[:,8]));heat_error=float(D*ds/(fuel*D_HEAT)-1)
  u=.5*(np.spacing(x0)+np.spacing(x1))
  rp=math.fsum(float(w[q])*math.fsum(float(u[q,j])*abs(float(BINDING[j]))for j in range(9))for q in range(len(w)))/ds
  mp=-math.fsum(float(w[q])*math.fsum(float(dx[q,j])*float(BINDING[j])for j in range(9))for q in range(len(w)))/ds
  diff=mp-release;L=after['data'][-1,4]
  checks=dict(D_heat=abs(heat_error)<2e-6,D_dominated=abs(release-D)<2e-6*release,
   mass_energy=abs(diff)<=2e-6*release+rp,total_energy=abs(diff)<2e-7*L,
   H=bool(np.all(abs(dx[:,0]-.5*dx[:,8]-ds*powers[i,end][3])<=u[:,0]+.5*u[:,8])),
   He3=bool(np.all(abs(dx[:,1]+1.5*dx[:,8]-ds*powers[i,end][4])<=u[:,1]+1.5*u[:,8])))
  checks={k:bool(v)for k,v in checks.items()}
  step.append(dict(start=start,end=end,checks=checks,fuel_heat_relative=heat_error,mass_energy_relative=diff/release,representation_relative=rp/release,mass_error_surface_fraction=float(diff/L)))
 full=models[i,'full']['data'];final=models[i,'final']['data']
 structure=float(np.max(abs(full[:,1:4]-final[:,1:4])));composition=float(np.max(abs(full[:,7:]-final[:,7:])))
 energy=float(abs(2*powers[i,'full'][0]/(powers[i,'half1'][0]+powers[i,'final'][0])-1))
 norm=max(structure/1e-4,composition/1e-8,energy/.005)
 assert abs(norm-a['error_norm'])<1e-8
 rows.append(dict(interval=i,passed=all(all(s['checks'].values())for s in step)and norm<=1,error_norm=norm,steps=step))
history=[json.loads(s)for s in (P/'history.jsonl').read_text().splitlines()]
assert len(history)==len(accepted)+1
assert all(history[i+1]['radius_Rsun']<history[i]['radius_Rsun'] and history[i+1]['mean_entropy']<history[i]['mean_entropy'] and history[i+1]['D_mass_g']<history[i]['D_mass_g']for i in range(len(accepted)))
result=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),all_checks_pass=all(r['passed']for r in rows),accepted_intervals=len(accepted),rejected_intervals=len(attempts)-len(accepted),years=history[-1]['years'],seed=history[0],final=history[-1],maximum_error_norm=max(r['error_norm']for r in rows),maximum_independent_D_heat_error=max(abs(s['fuel_heat_relative'])for r in rows for s in r['steps']),maximum_mass_error_surface_fraction=max(abs(s['mass_error_surface_fraction'])for r in rows for s in r['steps']),native_cpu_seconds=receipt['cpu_seconds'],native_wall_seconds=receipt['wall_seconds'],intervals=rows,selected_for_physical_PMS_track=False,limitations=['Homogeneous, D-dominated PMS audit with actual pp isotope contributions. Atmosphere source/interpolation acceptance and later physical regimes require separate assessment.','Independent replay evaluates source powers from saved structures; the unchanged first-law diagnostic is reported by the coupled structure solve.','No microscopic transport, nonuniform composition or later-burning regime is covered.'],input_sha256={str(p):sha(p)for p in files+[P/'plan.json',P/'result.json',P/'history.jsonl',P/'attempts.jsonl',probe_source,probe,W/'nuclear.txt',Path(__file__)]})
result['parent_checkpoint_verified']=str(args.parent_checkpoint.resolve()) if args.parent_checkpoint else None
result['stop_reason']=receipt['native_result']['stop_reason']
result['requested_age_reached']=receipt['returncode']==0 and result['stop_reason']=='requested duration'
path=args.output.resolve()
serialized=json.dumps(result,indent=2)+'\n'
with path.open('x')as f:f.write(serialized)
print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','intervals']},indent=2))
assert result['all_checks_pass']
