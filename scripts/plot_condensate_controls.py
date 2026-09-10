#!/usr/bin/env python3
"""Plot archived coupled condensate controls and their gas atmospheres."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from import_nongrey_grid import read_text
from validate_condensate_model import validate

ROOT=Path(__file__).resolve().parents[1]


def profile(directory):
    rows=[]
    for line in read_text(directory/'run.log.gz').rsplit('FINAL MODEL ATMOSPHERE',1)[1].splitlines():
        words=line.replace('D','E').split()
        if len(words)==11 and words[0].isdigit():rows.append(list(map(float,words)))
    return np.array(rows)


def main():
    cases=[('condensation_atmosphere_x700_2600_g515',r'$X_\mathrm{H}=0.70,\ X_3=0$ · 2600 K'),
           ('condensation_atmosphere_x700_2800_g515',r'$X_\mathrm{H}=0.70,\ X_3=0$ · 2800 K'),
           ('condensation_atmosphere_x300_he3120_2800_g515',r'$X_\mathrm{H}=0.30,\ X_3=0.12$ · 2800 K')]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(3,len(cases),figsize=(13,8),sharex=True,constrained_layout=True)
    for j,(name,title) in enumerate(cases):
        directory=ROOT/'data/atmosphere/sources'/name;result=validate(directory)
        cond=profile(directory);gas=profile(directory/'gas-control')
        axes[0,j].semilogx(gas[:,2],gas[:,3],label='Gas only',color='#7b7b7b',ls='--')
        axes[0,j].semilogx(cond[:,2],cond[:,3],label='Equilibrium depletion',color='#bb562e')
        axes[0,j].set_title(title);axes[0,j].set_ylabel('Temperature (K)')
        lower=max(gas[0,2],cond[0,2]);upper=min(gas[-1,2],cond[-1,2])
        tau=np.geomspace(lower,upper,600)
        for row,k,label,color in [(1,3,'Temperature','#bb562e'),(2,6,'Gas pressure','#246e91')]:
            c=np.exp(np.interp(np.log(tau),np.log(cond[:,2]),np.log(cond[:,k])))
            g=np.exp(np.interp(np.log(tau),np.log(gas[:,2]),np.log(gas[:,k])))
            axes[row,j].semilogx(tau,100*(c/g-1),label=label,color=color)
            boundary=np.exp(np.interp(np.log(100),np.log(cond[:,2]),np.log(cond[:,k])))/np.exp(
                np.interp(np.log(100),np.log(gas[:,2]),np.log(gas[:,k])))-1
            axes[row,j].text(.98,.92,rf'$\tau=100$: {100*boundary:+.3f}%',ha='right',va='top',
                             transform=axes[row,j].transAxes,color=color,
                             bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
            axes[row,j].set_ylabel(label+' change (%)')
            axes[row,j].axhline(0,color='black',lw=.5)
        depth=max(r['largest_tau'] for r in result['chemistry']['condensates'])
        for ax in axes[:,j]:
            ax.axvspan(lower,depth,color='#bb562e',alpha=.09)
            ax.axvline(100,color='black',lw=.8,ls=':');ax.grid(alpha=.15)
            ax.set_xlim(lower,upper)
        axes[2,j].set_xlabel(r'Rosseland optical depth $\tau$')
        axes[0,j].legend(frameon=False,fontsize=8)
    fig.suptitle('Coupled atmosphere controls · zero grain opacity · log g = 5.15\nShading extends to the deepest condensate; dotted line marks the interior boundary',fontsize=12)
    destination=ROOT/'docs/results/condensation_atmosphere_controls'
    fig.savefig(destination.with_suffix('.png'),dpi=180)
    fig.savefig(destination.with_suffix('.pdf'))


if __name__=='__main__':main()
