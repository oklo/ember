"""Plot the preserved contraction sequence; do not join it to later tracks."""
from pathlib import Path
import csv,hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
columns=['accepted','years','Teff_K','luminosity_Lsun','radius_Rsun','nuclear_fraction']
csvfile=HERE/'pms_main_sequence.csv'
inputs=[ROOT/'out/convective-transport-sept26-v1/v4/contraction/history.jsonl']
if inputs[0].exists():
 data=[{c:json.loads(line)[c] for c in columns} for line in inputs[0].read_text().splitlines()]
 assert [r['accepted'] for r in data]==list(range(1020))
 assert data[0]['years']==0 and data[-1]['years']==1e9
 assert all(b['years']>a['years'] for a,b in zip(data,data[1:]))
 with csvfile.open('w') as f:
  writer=csv.DictWriter(f,fieldnames=columns,lineterminator='\n');writer.writeheader();writer.writerows(data)
 (HERE/'pms_figure_inputs.json').write_text(json.dumps({'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},'count':len(data),'scope':'Common Hayashi-started D/pp/CN calculation with plasma losses, checked whole-star convection and composition heat after initial D; no radiative settling regime yet.'},indent=2)+'\n')
else:
 with csvfile.open() as f:data=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
a={k:np.array([r[k]for r in data])for k in columns}
plt.rcParams.update({'font.size':11,'axes.labelsize':12,'xtick.labelsize':10,'ytick.labelsize':10,'pdf.fonttype':42})
fig,(hr,power)=plt.subplots(1,2,figsize=(9.6,3.8),constrained_layout=True)
orange='#cb651f'
x=a['Teff_K'];y=np.log10(a['luminosity_Lsun']);age=a['years']/1e6
hr.plot(x,y,color=orange,lw=1.8);hr.set(xlabel=r'$T_{\rm eff}$ (K)',ylabel=r'$\log_{10}(L/L_\odot)$',xlim=(3060,2700),ylim=(-3.3,-.65))
hr.scatter([x[0],x[-1]],[y[0],y[-1]],c=orange,s=19,zorder=3)
for t,offset in [(0,(12,5)),(10,(14,7)),(100,(-70,4)),(1000,(6,6))]:
 i=int(np.argmin(abs(age-t)))
 hr.annotate('Start'if t==0 else ('1 Gyr'if t==1000 else f'{t} Myr'),(x[i],y[i]),xytext=offset,textcoords='offset points',fontsize=9,color='black')
for t0,t1 in [(3,5),(40,60),(130,190)]:
 i=int(np.argmin(abs(age-t0)));j=int(np.argmin(abs(age-t1)))
 hr.annotate('',(x[j],y[j]),(x[i],y[i]),arrowprops=dict(arrowstyle='-|>',color=orange,lw=1.5,mutation_scale=11))
positive=age>0
power.plot(age[positive],a['nuclear_fraction'][positive],color=orange,lw=1.8)
power.set(xscale='log',xlabel='Age from Hayashi start (Myr)',ylabel=r'$L_{\rm nuc}/L$',xlim=(.001,1500),ylim=(-.03,1.08))
power.axhline(1,color='.5',ls=':',lw=.8)
power.annotate('Initial D\nburning',(1.0,.91),xytext=(.012,.72),fontsize=9,arrowprops=dict(arrowstyle='-|>',lw=.7,color='black'))
power.annotate('H burning',(650,.93),xytext=(14,.79),fontsize=9,arrowprops=dict(arrowstyle='-|>',lw=.7,color='black'))
for ax in (hr,power):
 ax.spines[['top','right']].set_visible(False);ax.tick_params(direction='in');ax.grid(alpha=.13,lw=.4)
fig.savefig(HERE/'pms_main_sequence.pdf');fig.savefig(HERE/'pms_main_sequence.png',dpi=180)
print(json.dumps({'models':len(data),'nuclear_peak':float(a['nuclear_fraction'].max()),'age_at_peak_Myr':float(age[a['nuclear_fraction'].argmax()]),'endpoint':data[-1]}))
