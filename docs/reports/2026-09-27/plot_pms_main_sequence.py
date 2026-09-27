"""Plot the preserved continuous Hayashi-origin trajectory from its compact CSV."""
from pathlib import Path
import csv,hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
columns=['accepted','years','Teff_K','luminosity_Lsun','radius_Rsun','nuclear_fraction']
csvfile=HERE/'pms_main_sequence.csv'
with csvfile.open() as f:
 data=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
assert all(b['years']>a['years'] for a,b in zip(data,data[1:]))
a={k:np.array([r[k]for r in data])for k in columns}
plt.rcParams.update({'font.size':11,'axes.labelsize':12,'xtick.labelsize':10,'ytick.labelsize':10,'pdf.fonttype':42})
fig,(hr,power)=plt.subplots(1,2,figsize=(9.6,3.8),constrained_layout=True)
orange='#cb651f'
x=a['Teff_K'];y=np.log10(a['luminosity_Lsun']);age=a['years']/1e6
hr.plot(x,y,color=orange,lw=1.8);hr.set(xlabel=r'$T_{\rm eff}$ (K)',ylabel=r'$\log_{10}(L/L_\odot)$',xlim=(3340,2680),ylim=(-3.3,-.65))
hr.scatter([x[0],x[-1]],[y[0],y[-1]],c=orange,s=19,zorder=3)
for t,offset in [(0,(12,5)),(10,(14,7)),(1000,(-48,-17))]:
 i=int(np.argmin(abs(age-t)))
 hr.annotate('Start'if t==0 else ('1 Gyr'if t==1000 else f'{t} Myr'),(x[i],y[i]),xytext=offset,textcoords='offset points',fontsize=9,color='black',arrowprops=(dict(arrowstyle='-|>',lw=.65,color='black') if t==1000 else None))
for t0,t1 in [(3,5),(40,60),(130,190),(1.5e6,1.8e6),(2.8e6,3.0e6)]:
 i=int(np.argmin(abs(age-t0)));j=int(np.argmin(abs(age-t1)))
 hr.annotate('',(x[j],y[j]),(x[i],y[i]),arrowprops=dict(arrowstyle='-|>',color=orange,lw=1.5,mutation_scale=11))
hr.annotate(f"{age[-1]/1e6:.4g} Tyr",(x[-1],y[-1]),xytext=(12,15),textcoords='offset points',fontsize=9,ha='left',arrowprops=dict(arrowstyle='-|>',lw=.7,color='black'))
positive=(age>0)&(age<=1000)
power.plot(age[positive],a['nuclear_fraction'][positive],color=orange,lw=1.8)
power.set(xscale='log',xlabel='Age from Hayashi start (Myr)',ylabel=r'$L_{\rm nuc}/L$',xlim=(.001,1500),ylim=(-.03,1.08))
power.axhline(1,color='.5',ls=':',lw=.8)
power.annotate('Initial D\nburning',(1.0,.91),xytext=(.012,.72),fontsize=9,arrowprops=dict(arrowstyle='-|>',lw=.7,color='black'))
power.annotate('H burning',(650,.93),xytext=(14,.79),fontsize=9,arrowprops=dict(arrowstyle='-|>',lw=.7,color='black'))
for ax in (hr,power):
 ax.spines[['top','right']].set_visible(False);ax.tick_params(direction='in');ax.grid(alpha=.13,lw=.4)
fig.savefig(HERE/'pms_main_sequence.pdf');fig.savefig(HERE/'pms_main_sequence.png',dpi=180)
print(json.dumps({'models':len(data),'nuclear_peak':float(a['nuclear_fraction'].max()),'age_at_peak_Myr':float(age[a['nuclear_fraction'].argmax()]),'endpoint':data[-1]}))
