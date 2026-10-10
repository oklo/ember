"""Inspect observed populations before assigning evolutionary masses."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u

W=Path(__file__).resolve().parent
d=pd.read_csv(W/'catalogues/objects_with_luminosity.csv',low_memory=False)
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                     'savefig.dpi':180,'font.family':'DejaVu Sans'})
fig,axes=plt.subplots(2,2,figsize=(12,8),gridspec_kw={'width_ratios':[1.15,1]})
summary={}; colours=['#1c6b80','#d79138','#777777']
for row,cl in enumerate(['hyades','pleiades']):
    s=d[(d.cluster==cl)&d.selected_by_any].copy()
    inside=s.associations.str.contains('inside_tidal_radius',na=False)
    outer=s.associations.str.contains('tail|outside_tidal|greater_pleiades',na=False)&~inside
    s['population']=np.select([inside,outer],['Inside a published tidal radius','Tails / other associated candidates'],default='Other cluster candidates')
    valid=s.distance_for_luminosity_pc.between(0,1000)
    a=axes[row,0]
    # Heliocentric Galactic Cartesian coordinates, not a projected sky-area
    # correction. Missing individual distances remain absent from this panel.
    c=SkyCoord(ra=s.loc[valid,'ra_2016_deg'].to_numpy()*u.deg,
        dec=s.loc[valid,'dec_2016_deg'].to_numpy()*u.deg,
        distance=s.loc[valid,'distance_for_luminosity_pc'].to_numpy()*u.pc).galactic.cartesian
    for label,col in zip(['Other cluster candidates','Tails / other associated candidates','Inside a published tidal radius'],colours[::-1]):
        mask=s.loc[valid,'population'].eq(label).to_numpy()
        a.scatter(c.x.value[mask],c.y.value[mask],s=8,c=col,alpha=.45,rasterized=True,label=label)
    a.set_xlabel('Galactic X from the Sun (pc)');a.set_ylabel('Galactic Y from the Sun (pc)')
    a.text(.03,.96,cl.capitalize(),transform=a.transAxes,ha='left',va='top',fontsize=15,weight='bold')
    if row==0:a.legend(loc='lower left',fontsize=8.5,frameon=False)
    a.set_aspect('equal',adjustable='datalim')
    bins=np.arange(-5.2,-.999,.2);b=axes[row,1]
    hist={}
    for label,col in zip(['Inside a published tidal radius','Tails / other associated candidates','Other cluster candidates'],colours):
        lum=s.loc[s.population.eq(label),'logL_observed'].dropna()
        counts,_=np.histogram(lum,bins);hist[label]=counts.tolist()
        b.stairs(counts,bins,color=col,lw=1.6)
    b.set_xlabel(r'$\log_{10}(L/L_\odot)$');b.set_ylabel('Candidates per 0.2 dex')
    b.set_yscale('symlog',linthresh=2)
    b.set_ylim(0,2*max(max(v) for v in hist.values()))
    b.set_xlim(bins[0],bins[-1]);b.grid(axis='y',alpha=.15)
    b.text(.02,.97,'Observed counts; completeness correction pending',transform=b.transAxes,
           ha='left',va='top',fontsize=9)
    summary[cl]=dict(selected_candidates=len(s),with_individual_distance=int(valid.sum()),
        with_luminosity=int(s.logL_observed.notna().sum()),luminosity_bin_edges=bins.tolist(),counts=hist)
fig.tight_layout(pad=1.5)
fig.savefig(W/'catalogue_coverage.png');fig.savefig(W/'catalogue_coverage.pdf')
(W/'catalogue_coverage.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:{a:b for a,b in v.items() if a not in ['counts','luminosity_bin_edges']} for k,v in summary.items()},indent=2))
