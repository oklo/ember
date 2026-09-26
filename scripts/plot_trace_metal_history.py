#!/usr/bin/env python3
"""Plot the passive O/Fe transport controls; no stellar trajectory is altered."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prepare_nongrey_sources import digest
from write_scientific_result import write_result

def main():
    p=argparse.ArgumentParser();p.add_argument('pdf',type=Path);p.add_argument('receipt',type=Path);a=p.parse_args()
    if a.pdf.exists() or a.receipt.exists():raise FileExistsError('preserve figures')
    report=Path('docs/results/trace_metal_history_v1.json');d=json.loads(report.read_text())
    endpoint=Path('/tmp/ember-trace-metal-history-v1/endpoint.npz')
    assert digest(endpoint)==d['artifacts_sha256'][str(endpoint)]
    data=np.load(endpoint);profiles=data['fine_ratios'];q=data['mass']/data['mass'][-1]
    ages=np.array([d['initial_age_years']]+[r['age_years'] for r in d['history']])/1e12
    surface=np.r_[np.ones((1,2,8)),np.array([r['surface_ratio'] for r in d['history']])]
    plt.rcParams.update({'font.size':9,'axes.labelsize':10,'legend.fontsize':8,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.1),layout='constrained')
    for h,style in [(0,'-'),(1,'--')]:
        for k,name,color in [(2,'O','#3274a1'),(7,'Fe','#c17d25')]:
            label=f"{name}, {'with' if h==0 else 'without'} heat flow"
            axes[0].plot(ages,surface[:,h,k],style,color=color,lw=1.5,label=label)
            axes[1].plot(q,profiles[h,k],style,color=color,lw=1.5)
    axes[0].set(xlabel='Age (Tyr)',ylabel='Surface abundance / initial abundance',yscale='log',ylim=(.05,1.1))
    axes[0].set_xticks([3.55,3.60,3.65,3.70]);axes[0].legend(loc='lower left',frameon=False)
    axes[1].set(xlabel='Enclosed mass / stellar mass',ylabel='Abundance / initial abundance',xlim=(0,1),ylim=(0,2.5))
    for ax in axes:
        ax.tick_params(direction='in',top=True,right=True);ax.grid(alpha=.15)
    fig.savefig(a.pdf);plt.close(fig)
    write_result(a.receipt,dict(scope=__doc__,output=str(a.pdf),output_sha256=digest(a.pdf),
        input_sha256={str(p.resolve()):digest(p) for p in [Path(__file__),report,endpoint]},
        caption_scope='Passive fully stripped O/Fe on prescribed stellar structures; heat treatments are separate approximations, not uncertainty bounds.'))

if __name__=='__main__':main()
