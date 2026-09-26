"""Bounded grey/no-diffusion comparison using Ember's audited CN step control."""
from pathlib import Path
from datetime import datetime,timezone
import gzip,hashlib,json,math,resource,subprocess,sys,time
import numpy as np
R=Path('/Users/greglaughlin/Projects/ember');sys.path.insert(0,str(R/'scripts'))
from evolve_cn_transport import audit,compare,abundance_step_change
from audit_cn_thermal import weights,physical,YEAR
from grey_neutral_step import neutral_step
W=Path(sys.argv[1]);plan=json.loads((W/'plan.json').read_text());command=plan['native_command']
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):
 tmp=p.with_suffix(p.suffix+'.pending');tmp.write_text(json.dumps(d,indent=2,allow_nan=False,default=lambda v:v.item())+'\n');tmp.replace(p)
for p,h in plan['inputs_sha256'].items():assert sha(p)==h,p
assert not (W/'checkpoint.json').exists()
reference_path=Path(plan['accuracy_reference'])
reference=json.loads(reference_path.read_text())
control_age=reference['age_seconds'];control_checked=True
fuel_tolerance=plan['relative_hydrogen_time_tolerance']
assert fuel_tolerance==1e-4
assert sha(reference_path)==plan['accuracy_reference_sha256']
stderr=(W/'evolution.stderr').open('w');cpu0=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.monotonic()
p=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=stderr,text=True)
write(W/'process.json',dict(pid=p.pid,created_utc=datetime.now(timezone.utc).isoformat(),command=command))
def request(query):
 p.stdin.write(query+'\n');p.stdin.flush();line=p.stdout.readline()
 if not line:raise RuntimeError('native exited')
 return json.loads(line)
resume=json.loads(Path(plan['resume_checkpoint']).read_text())
record=resume['model_record']
write(W/'seed.json',record)
assert record.get('converged'),record
old=np.array(record['model']);w=weights(old[:,0]);mass=old[-1,0]
assert old.shape==(512,10) and abs(mass/1.98847e33-.1)<1e-5
seed_seconds=record['seconds'];seed_cpu=record['cpu_seconds']
history=[];attempts=[];solver_cpu=0.;solver_seconds=0.;calls=0;elapsed=0.;step=1e6*YEAR
def row(record,age,dt,error):
 m=np.array(record['model'])
 return dict(age_years=age/YEAR,interval_years=dt/YEAR,error_norm=error,Teff=record['Teff'],
   radius=m[-1,1],luminosity=m[-1,4],central_H=m[0,5],surface_H=m[-1,5],
   central_T=m[0,3],central_density=m[0,2],global_species=(w@physical(m)/mass).tolist(),
   mixed_mass_fraction=record['convective_mass_fraction'],solver_cpu_seconds=solver_cpu,
   cumulative_wall_seconds=time.monotonic()-start)
history=json.loads(Path(plan['previous_completion']).read_text())['history']
elapsed=resume['age_seconds'];step=5e5*YEAR
solver_cpu=json.loads(Path(plan['previous_completion']).read_text())['solver_cpu_seconds']
prior_solver_cpu_seconds=solver_cpu
with (W/'history.jsonl').open('w') as f:
 for prior_row in history:f.write(json.dumps(prior_row)+'\n')
def finite(model,dt,tol,previous):
 def exchange_request(face,parameter,m,years,tolerance):
  global solver_cpu,solver_seconds,calls
  q=' '.join(format(x,'.17g') for x in [1,face,parameter,years*YEAR,tolerance,len(m),*m.flat])
  r=request(q);calls+=1;solver_cpu+=r.get('cpu_seconds',0.);solver_seconds+=r.get('seconds',0.)
  if r.get('converged'):
   assert max(r['secular_diffusivity'],default=0)==0
  return r
 try:
  result,trace=neutral_step(model,dt/YEAR,previous,exchange_request,tol)
 except Exception as e:
  raise ValueError(str(e)) from e
 with (W/'boundary_trials.jsonl').open('a') as f:
  f.write(json.dumps(dict(age_years=elapsed/YEAR,interval_years=dt/YEAR,trials=trace))+'\n')
 return result

failure=None;outcome='running';accumulated=np.array(resume['accumulated_species_balance']);limit=10e12*YEAR
profiles=gzip.open(W/'accepted_models.jsonl.gz','wt')
try:
 while elapsed<limit:
  if len(history)>6000 or time.monotonic()-start>3600:
   outcome='bounded_comparison_limit';break
  if (W/'stop.request').exists():outcome='requested_checkpoint_stop';break
  h=min(step,limit-elapsed)
  if not control_checked:h=min(h,control_age-elapsed)
  attempt=dict(start_years=elapsed/YEAR,interval_years=h/YEAR)
  success=False;norm=4.
  for tol in [1e-12,1e-13]:
   try:
    _,full,c1,_=finite(old,h,tol,record)
    mid,half,c2,b1=finite(old,h/2,tol,record)
    new,fine,c3,b2=finite(mid,h/2,tol,half)
    diagnostics={}
    norm,criteria=compare(full,fine,'moving-boundary',diagnostics,record,half,.001)
    criteria['hydrogen']*=1e-6/fuel_tolerance
    norm=max(criteria.values())
    attempt.update(error_norm=norm,criteria=criteria,diagnostics=diagnostics,conservation=[c1,c2,c3],abundance_tolerance=tol)
    success=norm<=1
    break
   except ValueError as e:attempt.setdefault('solve_failures',[]).append(str(e))
  attempt['accepted']=success;attempts.append(attempt)
  with (W/'attempts.jsonl').open('a') as f:f.write(json.dumps(attempt)+'\n')
  if not success:
   recent=attempts[-4:]
   if len(recent)==4 and all(a.get('solve_failures') and not a['accepted'] for a in recent):
    msgs=[a['solve_failures'][-1] for a in recent]
    if len(set(msgs))==1:raise RuntimeError('unchanged rejection after four reductions: '+msgs[0])
   step=h*max(.1,min(.5,.8/math.sqrt(norm)))
   if step<YEAR:raise RuntimeError('minimum timestep; '+str(attempt))
   print('retry',f'{step/YEAR:.4g}',str(attempt.get('solve_failures',criteria if 'criteria' in locals() else ''))[-300:],flush=True)
   continue
  change=abundance_step_change(old,new);old=new;record=fine;elapsed+=h;accumulated+=b1+b2
  if not control_checked and abs(elapsed-control_age)<1:
   ref=np.array(reference['model_record']['model'])
   result=dict(relative_total_hydrogen=float((w@old[:,5])/(w@ref[:,5])-1),
     temperature_K=record['Teff']-reference['model_record']['Teff'],
     relative_luminosity=old[-1,4]/ref[-1,4]-1,
     relative_radius=old[-1,1]/ref[-1,1]-1,
     maximum_absolute_species=float(abs(physical(old)-physical(ref)).max()),
     comparison_age_years=elapsed/YEAR,reference_steps=len(json.loads(reference_path.with_name('completion.json').read_text())['history'])-1,
     reference_cpu_seconds=json.loads(reference_path.with_name('completion.json').read_text())['solver_cpu_seconds'],
     candidate_cpu_seconds=solver_cpu)
   passed=abs(result['relative_total_hydrogen'])<.001 and abs(result['temperature_K'])<2 and abs(result['relative_luminosity'])<.002
   result.update(passed=passed,scope='Matched evolution from identical initial model and microphysics; only total-hydrogen timestep comparison changes from 1e-6 to 1e-4. Species-mass and energy conservation limits unchanged. This is a finite early-phase check, not a whole-track bound.')
   write(W/'accuracy_control.json',result)
   write(R/'docs/results/ember_grey_time_accuracy_v2.json',result)
   print('ACCURACY CONTROL',json.dumps(result,default=lambda v:v.item()),flush=True)
   if not passed:raise RuntimeError('coarser time control failed the declared matched-interval limits')
   control_checked=True
  factor=min(2.,max(.5,.85/math.sqrt(max(norm,1e-6))),.85*.005/max(change,1e-30))
  step=h*factor;history.append(row(record,elapsed,h,norm))
  checkpoint=dict(age_seconds=elapsed,model_record=record,numerical_history_row=history[-1],next_step_seconds=step,
    accumulated_species_balance=accumulated.tolist(),plan_sha256=sha(W/'plan.json'))
  write(W/'checkpoint.json',checkpoint)
  profiles.write(json.dumps(checkpoint)+'\n');profiles.flush()
  with (W/'history.jsonl').open('a') as f:f.write(json.dumps(history[-1])+'\n')
  print(f"accepted {len(history)-1}: {elapsed/YEAR/1e12:.4g} Tyr, Teff {record['Teff']:.4g} K, Xc {old[0,5]:.4g}, CPU {solver_cpu/3600:.4g} h",flush=True)
 else:outcome='requested_age_reached'
except Exception as e:outcome='stopped';failure=str(e);print('STOP',failure,flush=True)
finally:
 profiles.close()
 try:p.stdin.close()
 except BrokenPipeError:pass
 code=p.wait(timeout=120);stderr.close()
cpu1=resource.getrusage(resource.RUSAGE_CHILDREN)
report=dict(created_utc=datetime.now(timezone.utc).isoformat(),outcome=outcome,failure=failure,native_exit_code=code,
 accepted_intervals=len(history)-1,rejected_intervals=sum(not a['accepted'] for a in attempts),
 history=history,seed_cpu_seconds=seed_cpu,seed_wall_seconds=seed_seconds,
 native_solves=calls,solver_cpu_seconds=solver_cpu,solver_wall_seconds=solver_seconds,
 total_wall_seconds=time.monotonic()-start,total_native_cpu_seconds=cpu1.ru_utime+cpu1.ru_stime-cpu0.ru_utime-cpu0.ru_stime,
 maximum_accumulated_species_balance=float(abs(accumulated.sum(axis=0)).max()),
 accepted_as_primary_track=False,comparison_validation='Independent step and conservation checks; EOS source controls and retained-profile review required before scientific selection.',
 prior_solver_cpu_seconds=prior_solver_cpu_seconds,plan=plan,driver_sha256=sha(__file__),endpoint=history[-1])
for path,h in plan['inputs_sha256'].items():assert sha(path)==h,path
write(W/'completion.json',report)
write(R/'docs/results/ember_grey_no_diffusion_v8.json',report)
print(json.dumps({k:report[k] for k in ['outcome','failure','accepted_intervals','solver_cpu_seconds','total_native_cpu_seconds','total_wall_seconds']},indent=2))
