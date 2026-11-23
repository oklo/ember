"""One bounded birth control for a mass whose initial atmosphere was out of range."""
import fcntl,json,math,os,re,subprocess,sys,tempfile
from pathlib import Path
from coordinate import root,out,reference,read,write,solver,campaign,dirs,locked,next_step
def run(m):
    cases=dirs(m)
    if not cases or next_step(m)[0] not in ('inspect interrupted run','atmosphere birth control'):raise ValueError('case is not a failed birth')
    if any((c/'evolution/history.jsonl').is_file() for c in cases):raise ValueError('existing history')
    if any(read(c/'launch.json')['action']=='atmosphere birth control' for c in cases):return
    if not any('pressure-matched atmosphere' in (c/'run.log').read_text() for c in cases):raise ValueError('unrelated failure')
    cfg=(reference/'solar-birth.cfg').read_text()
    def setting(s,key,val):return re.sub('(?m)^'+key+' ".*"$',f'{key} "{val}"',s)
    cfg=setting(setting(cfg,'mass_Msun',m/100),'zone_threads',1)
    # Keep the initial gravity guess inside the existing solar atmosphere for the light BDs.
    cfg=setting(cfg,'initial_radius_Rsun',math.sqrt(27420*(m/100)/10**3.82))
    if m>=9:
        # This optical-depth family explicitly interpolates in Z. Only the atmosphere
        # stack is transferred. Stellar composition and interior come from the solar run.
        legacy=out/'solar-1k-grid-nov17/m0300/cases/fresh-first-step/config.txt'
        keys=('atmosphere_main_sequence','atmosphere_join_','atmosphere_metal_chain','atmosphere_metal_join_','atmosphere_hydrogen_')
        cfg=re.sub('(?m)^atmosphere_pressure_matched .*\n','',cfg)
        cfg+='\n'+'\n'.join(l for l in legacy.read_text().splitlines() if l.startswith(keys))+'\n'
        mode='metal-interpolated optical-depth atmosphere; provisional model comparison'
    else:mode='original solar atmosphere; initial radius guess at log g=3.82'
    staging=Path(tempfile.mkdtemp(prefix='.birth-',dir=cases[-1].parent))
    dest=cases[-1].parent/f's{len(cases):04}'
    with (staging/'worker.lock').open('a') as lock,(staging/'run.log').open('w') as log:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        (staging/'config.txt').write_text(cfg)
        cmd=[str(solver),'--lifetime',str(dest/'config.txt'),'1e11',str(dest/'evolution'),'--cpu-seconds','600']
        info=dict(mass_Msun=m/100,initial_Z=campaign['initial_Z'],short_step=False,start_age_years=0,
                  action='atmosphere birth control',control=mode,command=cmd)
        write(staging/'launch.json',info);staging.rename(dest)
        child=subprocess.Popen(cmd,stdout=log,stderr=log,pass_fds=[lock.fileno()],start_new_session=True,env=dict(os.environ,OMP_NUM_THREADS='1'))
        write(dest/'launch.json',dict(info,pid=child.pid))
if __name__=='__main__':run(int(sys.argv[1]))
