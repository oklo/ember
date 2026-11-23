"""Bounded diagnostic of a 50 K extension past the 2800 K solar birth grid."""
import hashlib,json,math,re,shlex,subprocess,sys
from pathlib import Path
from coordinate import root,out,reference,solver,read,write
from atmosphere_repair import score_column

d=root/'birth-atmosphere-2850'
s=read(d/'spec.json');settings=(root/'m09/s0000/config.txt').read_text()
family=Path(dict(shlex.split(l) for l in settings.splitlines() if l.strip())['atmosphere_pressure_matched'])
def extend():
    body=[]
    for line in family.read_text().splitlines():
        if not line.startswith('file '):
            body.append(line);continue
        p=Path(shlex.split(line)[1]);v=p.read_text().splitlines()
        t=list(map(float,v[4].split()[2:]));g=list(map(float,v[5].split()[2:]));ng=len(g)
        assert abs(10**t[-1]-2800)<.01 and len(v[7:])==len(t)*ng
        end=math.log10(2850);w=(end-t[-1])/(t[-1]-t[-2])
        v[4]='log_teff '+str(len(t)+1)+' '+' '.join(map(str,t+[end]))
        new=[]
        for i in range(ng):
            a,b=[list(map(float,v[7+(len(t)-j)*ng+i].split())) for j in [2,1]]
            if sys.argv[1]=='predict-c1':
                prev=list(map(float,v[7+(len(t)-3)*ng+i].split()))
                h0,h1,dx=t[-1]-t[-2],t[-2]-t[-3],end-t[-1]
                vals=[]
                for j in range(3):
                    y0,y1,y2=[math.log(z[j]) if j==2 else z[j] for z in [prev,a,b]]
                    s0,s1=(y2-y1)/h0,(y1-y0)/h1
                    slope=((2*h0+h1)*s0-h0*s1)/(h0+h1)
                    if slope*s0<=0: slope=0.
                    elif s0*s1<=0 and abs(slope)>3*abs(s0):slope=3*s0
                    w1,w2=2*dx+h0,dx+2*h0
                    # Invert the reader's harmonic derivative at the new interior knot.
                    sec=w2/((w1+w2)/slope-w1/s0) if slope else 0.
                    assert not slope or sec*slope>0
                    y=y2+sec*dx
                    vals.append(math.exp(y) if j==2 else y)
                new.append(' '.join(format(y,'.17g') for y in vals))
            else:
                new.append(f'{b[0]+w*(b[0]-a[0]):.17g} {b[1]+w*(b[1]-a[1]):.17g} {b[2]*(b[2]/a[2])**w:.17g}')
        dest=d/('extended-'+p.name);dest.write_text('\n'.join(v+new)+'\n')
        body.append(f'file "{dest}"')
    f=d/'trial-family.dat';f.write_text('\n'.join(body)+'\n')
    args=read(out/'highz-native-opacity-oct12-v1/matching_depth_controls.json')['command'][1:]
    def query(table,te):
        a=s['arms']['anchor']
        q=subprocess.run([out/'boundary-corrected-item1-oct13/probe',*args,table],text=True,capture_output=True,
            input=f"thin {te} {s['logg']} {a['X']} {a['He3']} {a['Z']}\n",check=True)
        val=list(map(float,q.stdout.split()[6:9]));assert len(val)==3 and all(x>0 for x in val);return val
    for te in (2750,2799):
        assert all(abs(a/b-1)<1e-9 for a,b in zip(query(family,te),query(f,te)))
    write(d/'prediction.json',query(f,s['teff']))
    (d/'trial.cfg').write_text(settings.replace(str(family),str(f)))
if sys.argv[1] in ('predict','predict-c1'):extend()
elif sys.argv[1]=='assess':
    result=score_column(d);write(d/'assessment.json',result);sys.exit(0 if result['passes'] else 1)
elif sys.argv[1]=='pilot':
    # A diagnostic remains separate from the adopted track. Run even if the
    # independent atmosphere is still computing; promotion requires both checks.
    q=subprocess.run([solver,'--lifetime',d/'trial.cfg','1e11',d/'pilot','--max-steps','12','--cpu-seconds','180'])
    report=d/'pilot/report.json'
    if report.exists():r=read(report)
    else:r={'reason':'solver exited before a model was accepted'}
    write(d/'pilot-result.json',dict(returncode=q.returncode,report=r))
    sys.exit(0 if r.get('accepted_this_invocation',0)>0 else 1)
