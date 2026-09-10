#!/usr/bin/env python3
"""Compare local heat transport in saved 0.1-solar-mass models.

Native logarithms are converted only for an offline diagnostic. This script
never constructs or changes a restart file. The printed nodal gradients do
not replace the evolution driver's face-based Ledoux mixing criterion.
"""
import json,math,subprocess,hashlib
from pathlib import Path
work=Path('/tmp/ember-convection-explain-v1');work.mkdir(exist_ok=True)
checkpoint=Path('out/evolution-cold-remnant-x015-transition-512-3560gyr-v1.restart')
lines=checkpoint.read_text().splitlines();offset=8+int(lines[7]);head=lines[offset].split();n=int(head[0]);rows=[]
for line in lines[offset+1:offset+1+n]:
 v=list(map(float,line.split()));rows.append([v[0],math.exp(v[1]),math.exp(v[2]),math.exp(v[3]),v[4],*v[7:10]])
profiles=[]
for age in [3300,3400]:
 p=Path(f'out/evolution-cold-remnant-forward-512-{age}gyr-v1.json');d=json.loads(p.read_text());profiles.append((str(p),d['history'][-1][0],d['profile']))
profiles.append((str(checkpoint),float(head[2])/31557600,rows))
records=[]
for file,age,rows in profiles:
 text=''.join(' '.join(format(v,'.17g') for v in row)+'\n' for row in rows)
 command=['/tmp/ember-evolution-physics-exhaustion-probe-v1','--gs98','--eos-family','data/eos/exhaustion_refined_v2/freeeos300_gs98_z020.dat','--opacity-directory','data/opacity/hydrogen_poor_refined_v4']
 r=subprocess.run(command,input=text,text=True,capture_output=True,check=True);values=[json.loads(l) for l in r.stdout.splitlines()];idx=min(range(len(rows)),key=lambda i:abs(rows[i][0]/rows[-1][0]-.21));v=values[idx];row=rows[idx]
 record={'input':file,'input_sha256':hashlib.sha256(Path(file).read_bytes()).hexdigest(),'age_yr':age,'central_X':rows[0][5],'row':idx,'mass_fraction':row[0]/rows[-1][0],'radius_fraction':row[1]/rows[-1][1],'rho':row[2],'T_K':row[3],'L':row[4],'X':row[5],'P':v[0],'grad_ad':v[3],'radiative_opacity':v[4],'conductive_opacity':v[5],'combined_opacity':v[6],'required_gradient':v[7],'required_over_adiabatic':v[7]/v[3],'conductive_share_of_diffusive_transport':v[6]/v[5],'central_required_over_adiabatic':values[0][7]/values[0][3]}
 (work/(Path(file).stem+'.profile')).write_text(text);(work/(Path(file).stem+'.physics')).write_text(r.stdout)
 records.append(record)
 if file==str(checkpoint):
  cmd=['/tmp/ember-evolution-mixing-probe-v1','data/eos/exhaustion_refined_v2/freeeos300_gs98_z020.dat','data/opacity/hydrogen_poor_refined_v4']
  mixing=subprocess.run(cmd,input=text,text=True,capture_output=True,check=True);regions=[json.loads(l) for l in mixing.stdout.splitlines()];combined=[]
  for x in regions:
   x['inner_radius_fraction']=rows[x['begin']][1]/rows[-1][1];x['outer_radius_fraction']=rows[x['end_exclusive']-1][1]/rows[-1][1]
   if combined and not combined[-1]['convective'] and not x['convective']:
    c=combined[-1];c['end_exclusive']=x['end_exclusive'];c['mass_fraction']+=x['mass_fraction'];c['outer_enclosed_fraction']=x['outer_enclosed_fraction'];c['outer_radius_fraction']=x['outer_radius_fraction']
   else:combined.append(x)
  record['regions']=combined
print(json.dumps(records,indent=2));(work/'diagnostics.json').write_text(json.dumps(records,indent=2)+'\n')
