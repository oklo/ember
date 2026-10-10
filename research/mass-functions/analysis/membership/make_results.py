"""Plot the membership audit and preserve queryable numerical results."""
from pathlib import Path
import json,sqlite3,hashlib,shutil,platform
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter,ScalarFormatter
from scipy.interpolate import PchipInterpolator
from scipy.special import ndtr
W=Path(__file__).resolve().parent;ROOT=W.parents[1]
DOC=ROOT/'docs/research/independent';STEM='CLUSTER_MEMBERSHIP_TEST_20261009'
BLUE='#315c80';ORANGE='#c57530';TEAL='#408a84';GRAY='#89939b'

def main():
    iso=pd.read_csv(W/'isochrones_with_temperature.csv')
    rec=pd.read_csv(W/'candidate_recovery_counts.csv')
    sim=pd.read_csv(W/'cp_injection_summary.csv')
    group=pd.read_csv(W/'pleiades_group_partition.csv',dtype={'GaiaID':str})
    totals=pd.read_csv(W/'catalogue_count_audit.csv')
    plt.rcParams.update({'font.size':12,'axes.labelsize':12,'axes.spines.top':False,'axes.spines.right':False,
                         'legend.frameon':False,'savefig.facecolor':'white','font.family':'DejaVu Sans'})
    fig,axes=plt.subplots(2,2,figsize=(13.5,9.5),gridspec_kw={'hspace':.37,'wspace':.28})
    ax=axes[0,0]
    ax.axvspan(.03,.08,color=ORANGE,alpha=.09)
    ax.axvspan(.1,.2,color=BLUE,alpha=.08)
    for cloud,color,style in [('nc',BLUE,'-'),('f2',TEAL,'--')]:
        q=iso[(iso.cluster=='pleiades')&(iso.cloud==cloud)].sort_values('mass_Msun')
        m=np.geomspace(.03,.25,500)
        ax.plot(m,PchipInterpolator(q.mass_Msun,q.Teff_K)(m),color=color,lw=2.2,ls=style,label='Clear' if cloud=='nc' else 'Cloud clearing')
    ax.axhline(3000,color='black',lw=1,ls=':')
    ax.text(.033,3095,r'Complex survey: estimated $T_{\mathrm{eff}}>3000$ K',fontsize=11)
    ax.annotate(r'$0.1071\,M_\odot$',xy=(.1071245,3000),xytext=(.144,2540),fontsize=11,
                arrowprops={'arrowstyle':'->','lw':.9,'color':'black'})
    ax.text(.047,1430,'BD interval',ha='center',fontsize=11)
    ax.text(.145,1430,'Reference stars',ha='center',fontsize=11)
    ax.set(xscale='log',xlim=(.028,.27),ylim=(1250,3575),xlabel=r'Mass ($M_\odot$)',ylabel='Effective temperature at 125 Myr (K)')
    ax.set_xticks([.03,.05,.08,.1,.15,.2,.25]);ax.xaxis.set_major_formatter(ScalarFormatter());ax.minorticks_off()
    ax.legend(loc='upper left',bbox_to_anchor=(0,1),fontsize=10,ncol=2)

    ax=axes[0,1]
    x=group.xp+8122;y=group.yp
    ax.scatter(x,y,s=8,color=GRAY,alpha=.30,linewidths=0,rasterized=True)
    p=group.expected_reference.fillna(0).to_numpy();only=group.boyle_only.fillna(False).to_numpy(dtype=bool)
    chosen=only&(p>.01)
    ax.scatter(x[chosen],y[chosen],s=16+45*p[chosen],color=ORANGE,alpha=.65,linewidths=.2,edgecolors='white',rasterized=True)
    centers={'Pleiades':(-120.37,29.33),'UPK 303':(-156.10,110.72),'AB Dor':(-25.37,16.53),
             'HSC 1964':(-33.85,-106.80),'UPK 545':(28.53,-320.95),'GPC-1':(10.15,118.14),'GPC-2':(148.94,-101.76)}
    offsets={'Pleiades':(-12,-27),'UPK 303':(-15,17),'AB Dor':(16,-5),'HSC 1964':(-18,-31),'UPK 545':(14,-12),'GPC-1':(10,18),'GPC-2':(-2,19)}
    for label,(xx,yy) in centers.items():
        ax.plot(xx,yy,'+',color='black',ms=7,mew=1)
        dx,dy=offsets[label]
        ax.annotate(label,(xx,yy),xytext=(dx,dy),textcoords='offset points',fontsize=9.7,ha='center',
                    arrowprops={'arrowstyle':'-','lw':.5,'color':'#555555'})
    ax.plot(0,0,'*',color='black',ms=8)
    ax.set(xlim=(-325,325),ylim=(-415,330),xlabel='Position toward Galactic center (pc)',ylabel='Position along Galactic rotation (pc)')
    ax.text(.97,.96,'Greater Pleiades Complex',transform=ax.transAxes,ha='right',va='top',fontsize=11)
    ax.text(.025,.12,'Orange: additional reference stars\nfrom complex catalogue: 86.23',transform=ax.transAxes,ha='left',va='top',fontsize=10,
            bbox={'facecolor':'white','edgecolor':'none','alpha':.85,'pad':3})

    ax=axes[1,0];stages=[1,4,6,7]
    for cl,style,marker in [('hyades','-','o'),('pleiades','--','s')]:
        for pop,color in [('reference',BLUE),('BD',ORANGE)]:
            q=rec[(rec.cluster==cl)&(rec.population==pop)].set_index('stage').expected_count
            name=cl.capitalize()+(' stars' if pop=='reference' else ' BDs')
            ax.plot(range(4),q.loc[stages].to_numpy()/q[1],ls=style,color=color,marker=marker,ms=6,lw=1.8,label=name)
    ax.set_xticks(range(4),['Gaia\ncounterparts','Quality\ncuts','Distance +\nsearch volume','Automatic\nmembership'])
    ax.set(ylim=(0,1.12),xlim=(-.1,3.13),ylabel='Fraction of catalogued Gaia counterparts')
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.legend(loc='lower left',fontsize=9.7,ncol=2,columnspacing=1)
    ax.text(.04,.95,'Observed candidates',transform=ax.transAxes,va='top',fontsize=11)
    ax.grid(axis='y',alpha=.14)

    ax=axes[1,1]
    runs=['hyades_r2.5_assort_s43','pleiades_r3_assort_s43','pleiades_r4_fieldlow_s43']
    for pop,col,offset in [('reference',BLUE,-.17),('BD',ORANGE,.17)]:
        y=[];lo=[];hi=[]
        for run in runs:
            q=sim[(sim.run==run)&sim.astrometric_noise&(sim.region=='all')&(sim.population==pop)].CP_recovery_given_Gaia
            y.append(q.mean());lo.append(q.mean()-q.min());hi.append(q.max()-q.mean())
        xx=np.arange(3)+offset
        ax.bar(xx,y,width=.3,color=col,alpha=.8,label='Reference stars' if pop=='reference' else 'Brown dwarfs')
        ax.errorbar(xx,y,yerr=[lo,hi],fmt='none',color='black',capsize=4,lw=1)
        for x0,y0 in zip(xx,y):ax.text(x0,y0+.035,f'{100*y0:.1f}%',ha='center',fontsize=10)
    ax.set_xticks(range(3),['Hyades\n'+r'$R_{\rm birth}=0.63$','Pleiades\n'+r'$R_{\rm birth}=0.63$','Pleiades\n'+r'$R_{\rm birth}=0.55$'])
    ax.set(ylim=(0,1.28),ylabel='Automatic recovery of simulated Gaia detections')
    ax.set_yticks(np.arange(0,1.01,.2));ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.legend(loc='upper left',fontsize=10,ncol=2)
    ax.text(.03,.855,'Faint-source quality losses not yet applied',transform=ax.transAxes,fontsize=9.7)
    ax.grid(axis='y',alpha=.14);ax.set_axisbelow(True)
    for letter,ax in zip('abcd',axes.flat):ax.text(-.12,1.03,letter,transform=ax.transAxes,fontweight='bold',fontsize=14)
    fig.subplots_adjust(left=.085,right=.975,bottom=.09,top=.955)
    for ext in ['png','pdf']:fig.savefig(DOC/f'{STEM}.{ext}',dpi=190)
    plt.close(fig)

    # Illustrative temperature-error sensitivity, not a measured TIC error model.
    scatter=[]
    for cloud,g in iso[iso.cluster=='pleiades'].groupby('cloud'):
        m=np.geomspace(.03,.08,2000);t=PchipInterpolator(g.mass_Msun,g.Teff_K)(m)
        for sigma in [100,200]:
            eligibility=ndtr((t-3000)/sigma)
            scatter.append(dict(cloud=cloud,sigma_T_assumed_K=sigma,
                                BD_probability_logmass_average=float(np.trapezoid(eligibility,np.log(m))/np.log(.08/.03)),
                                probability_at_008=float(eligibility[-1])))
    pd.DataFrame(scatter).to_csv(W/'temperature_scatter_sensitivity.csv',index=False)
    # Save the proper comparison to the CP catalogue, without pretending that
    # automatic mock counts already include its quality and manual selections.
    published=pd.read_csv(W/'published_CP_counts.csv')
    comparison=sim[(sim.region=='all')&sim.astrometric_noise].merge(published,on=['cluster','population'])
    comparison.to_csv(W/'conditional_CP_comparison.csv',index=False)
    tablepaths=['catalogue_count_audit','pleiades_groups','isochrones_with_temperature','cp_injection_summary',
                'candidate_recovery_counts','published_CP_counts','conditional_CP_comparison','temperature_scatter_sensitivity']
    for name in tablepaths:shutil.copyfile(W/f'{name}.csv',DOC/f'cluster_membership_{name}_20261009.csv')
    tables=tablepaths+['candidate_selection_audit','pleiades_group_partition','cp_injection_draws','candidate_recovery_details']
    with sqlite3.connect(W/'cluster_membership_test.sqlite') as con:
        for name in tables:
            pd.read_csv(W/f'{name}.csv',dtype={'GaiaID':str,'gaia_dr3_id':str,'source_id':str}).to_sql(name,con,if_exists='replace',index=False)
        for cl in ['hyades','pleiades']:
            pd.read_csv(W/f'{cl}_cp_reproduction.csv',dtype={'GaiaDR3':str}).to_sql(cl+'_CP_reproduction',con,if_exists='replace',index=False)
        con.execute('CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY,value TEXT)')
        entries={'scope':'Conditional selection tests, not a recovered birth mass function',
                 'target':'All associated members and tidal tails remain in the retained candidate union',
                 'CP_missing_stages':'Mock faint-source quality, manual isolation; Hyades SCCP omitted',
                 'group_partition':'Nearest published 3D center only, not published group membership',
                 'temperature_scatter':'100/200 K scenarios assumed, not measured; uniform log mass weighting',
                 'population_bins':'BD=.03-.08 Msun; reference=.1-.2 Msun; fractional candidate counts are mass-draw means',
                 'parent_selection':'Gaia quality-selected field; remove all known cluster candidates before injection',
                 'CP_reproduction':(W/'cp_reproduction.json').read_text(),
                 'execution':(W/'cp_injection_execution.json').read_text()}
        con.executemany('INSERT OR REPLACE INTO metadata VALUES (?,?)',entries.items())
        con.execute('CREATE INDEX IF NOT EXISTS candidates_by_object ON candidate_recovery_details(object_id)')
        assert con.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
    validate(sim,rec)

def validate(sim,rec):
    d=pd.read_csv(W/'cp_injection_draws.csv')
    assert len(d)==2880
    counts=d[['intrinsic_sources','gaia_detected','automatic_window','CP_recovered']].to_numpy()
    assert (counts>=0).all() and (np.diff(counts,axis=1)<=0).all()
    assert d.groupby(['run','cloud','astrometric_noise','population','trial']).region.nunique().eq(5).all()
    for key,g in rec.groupby(['cluster','population']):assert np.all(np.diff(g.sort_values('stage').expected_count)<=1e-10)
    audit=pd.read_csv(W/'catalogue_count_audit.csv')
    for pop in ['BD','reference']:
        a=audit[(audit.cluster=='pleiades')&(audit.kind==pop)].set_index('sample').expected_count
        assert abs(a['all_associated']-a['boyle_only']-a['other_catalogues'])<1e-9
    q=json.loads((W/'candidate_gaia_receipt.json').read_text());assert q['requested']==q['returned']==5415
    repro=json.loads((W/'cp_reproduction.json').read_text())
    assert repro[0]['automatic_recovered_CP_members']==922 and repro[1]['automatic_recovered_CP_members']==2028
    raw_verified=0
    for receipt in (W/'raw').glob('*.json'):
        metadata=json.loads(receipt.read_text());raw=receipt.with_suffix('.csv')
        if raw.exists() and 'sha256' in metadata:
            assert hashlib.sha256(raw.read_bytes()).hexdigest()==metadata['sha256']
            raw_verified+=1
    # Low-precision source IDs are catastrophic: ensure exact round trip in SQLite.
    with sqlite3.connect(W/'cluster_membership_test.sqlite') as con:
        ids=pd.read_csv(W/'pleiades_group_partition.csv',dtype={'GaiaID':str}).GaiaID
        got={r[0] for r in con.execute('SELECT GaiaID FROM pleiades_group_partition')}
        assert set(ids)==got
    paths=list(W.glob('*.py'))+list(W.glob('*.csv'))+list((W/'raw').glob('*.json'))
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    receipt=dict(completed_UTC=pd.Timestamp.now(tz='UTC').isoformat(),injection_trials=288,count_rows=len(d),
                 candidate_Gaia_requested=5415,candidate_Gaia_returned=5415,count_order_checks=True,
                 exact_ID_roundtrip=True,source_union_partition_check=True,database_integrity='ok',
                 raw_download_hashes_verified=raw_verified,
                 no_new_collisional_runs=True,CP_members_recovery=[922/975,2028/2029],
                 dependencies=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,matplotlib=matplotlib.__version__),
                 hashes=hashes)
    (W/'validation_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('Saved figure, tables and database. Validation passed.')

if __name__=='__main__':main()
