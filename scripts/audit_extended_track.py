#!/usr/bin/env python3
"""Audit an extended track's transport, source coverage and plasmon losses."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
from stellar_composition import interior_composition


def plasmon(T,rho,ye):
    # Haft, Raffelt & Weiss 1994, eqs. 23--27; source sin^2(theta_W)=.23.
    # Fully ionized approximation. This evaluates only the plasma process.
    lam=1.686e-10*T
    gamma=math.sqrt(1.1095e11*rho*ye/(T*T*math.sqrt(1+(1.019e-6*rho*ye)**(2/3))))
    ft=2.4+.6*math.sqrt(gamma)+.51*gamma+1.25*gamma**1.5
    fl=(8.6*gamma**2+1.35*gamma**3.5)/(225-17*gamma+gamma**2)
    x=(17.5+math.log10(2*rho*ye)-3*math.log10(T))/6
    y=(-24.5+math.log10(2*rho*ye)+3*math.log10(T))/6
    fxy=1 if abs(x)>.7 or y<0 else 1.05+(.39-1.25*x-.35*math.sin(4.5*x)-.3*math.exp(-(4.5*x+.9)**2))*math.exp(-(min(0,y-1.6+1.25*x)/(.57-.25*x))**2)
    return .9248*3e21*lam**9*gamma**6*math.exp(-gamma)*(ft+fl)*fxy/rho


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('probe',type=Path);ap.add_argument('track',type=Path);ap.add_argument('output',type=Path);a=ap.parse_args()
    track=json.loads(a.track.read_text())
    if not track['converged']:raise ValueError('audit requires a completed track')
    profile=track['profile'];request=''.join(' '.join(f'{v:.17g}' for v in row)+'\n' for row in profile)
    metal=track.get('metal_inventory','').lower().startswith('gs98')
    p=subprocess.run([str(a.probe.resolve())]+(['--gs98'] if metal else []),input=request,text=True,capture_output=True,check=True)
    local=[json.loads(line) for line in p.stdout.splitlines()]
    if len(local)!=len(profile):raise ValueError('incomplete transport audit')
    weights=[0.]*len(profile);weights[0]=profile[0][0]
    for i in range(1,len(profile)):
        dm=profile[i][0]-profile[i-1][0];weights[i-1]+=.5*dm;weights[i]+=.5*dm
    loss=0;maxloss=0;rest=0;minimum_margin=math.inf
    eos_source_h=[];opacity_source_h=[];opacity_source_z=[]
    initial=interior_composition();integers=[1,3,4,12,13,14,16,20];atomic=[1.00782503,3.01602932,4.00260325,12,13.00335484,14.003074,15.9949146,20]
    if metal:
        from generate_nongrey_grid import composition,SOURCES
        numbers,_=composition(.7,0,initial[3:])
        elements=json.loads((SOURCES/'synple-elements.json').read_text())
        metal_ye=sum((i+1)*.7*v for i,v in enumerate(numbers) if i>=2)
        source_metal_mass=sum(.7*v*m for v,m in zip(numbers[2:],elements['mass'][2:],strict=True))
    for i,(row,state,w) in enumerate(zip(profile,local,weights,strict=True)):
        m,r,rho,T,L,X,Y3,Y4=row;c=interior_composition(X);c[1]=Y3;c[2]=Y4
        ye=sum(z*x/aa for z,x,aa in zip([1,2,2,6,6,7,8,10],c,integers,strict=True))
        if metal:ye=X+2*Y3/3+Y4/2+metal_ye
        q=plasmon(T,rho,ye);loss+=w*q;maxloss=max(maxloss,q)
        rest+=w*sum((aa/ai-1)*(old-new) for aa,ai,old,new in zip(atomic,integers,initial,c,strict=True))
        eos_scale=atomic[0]*X+atomic[2]*(Y3/3+(Y4+sum(c[3:]))/4)
        metal_mass=sum(aa*x/ai for aa,x,ai in zip(atomic[3:],c[3:],integers[3:],strict=True))
        if metal:metal_mass=source_metal_mass
        opacity_scale=atomic[0]*X+atomic[2]*(Y3/3+Y4/4)+metal_mass
        eos_source_h.append(X if metal else atomic[0]*X/eos_scale)
        opacity_source_h.append(atomic[0]*X/opacity_scale)
        opacity_source_z.append(metal_mass/opacity_scale)
        if i:
            old=profile[i-1];s=local[i-1];Tb=.5*(T+old[3]);mb=.5*(m+old[0]);Pb=.5*(state[0]+s[0]);kb=.5*(state[6]+s[6]);Lb=.5*(L+old[4]);ad=.5*(state[3]+s[3]);gr=3*kb*Lb*Pb/(16*math.pi*7.565733250033928e-15*2.99792458e10*6.6743e-8*mb*Tb**4)
            minimum_margin=min(minimum_margin,gr/ad)
    summary={
        'description':'Final-profile transport and source-domain audit; plasmon fit only, not total thermal neutrinos',
        'input_sha256':hashlib.sha256(a.track.read_bytes()).hexdigest(),'points':track['points'],
        'age_yr':track['history'][-1][0],
        'metal_inventory':track.get('metal_inventory','carried_isotopes'),
        'central':dict(zip(['P','E','S','grad_ad','radiative_opacity','conductive_opacity','combined_opacity','gradient_if_all_flux_diffusive','mlt_gradient','estimated_conductive_flux_fraction','opacity_rho_min','opacity_rho_max','classic_conductive_opacity','undamped_conductive_opacity'],local[0],strict=True)),
        'minimum_midpoint_diffusive_gradient_over_adiabatic':minimum_margin,
        'maximum_estimated_local_conductive_flux_fraction':max(s[9] for s in local),
        'minimum_opacity_density_above_floor_factor':min(r[2]/s[10] for r,s in zip(profile,local,strict=True)),
        'minimum_opacity_density_below_ceiling_factor':min(s[11]/r[2] for r,s in zip(profile,local,strict=True)),
        'mapped_source_composition':{
            'EOS_H_min_max':[min(eos_source_h),max(eos_source_h)],
            'opacity_H_min_max':[min(opacity_source_h),max(opacity_source_h)],
            'opacity_Z_min_max':[min(opacity_source_z),max(opacity_source_z)],
        },
        'final_boundary':{'Teff_K':track['history'][-1][4],
                          'log10_g':math.log10(6.6743e-8*profile[-1][0]/profile[-1][1]**2)},
        'max_species_spread':max(max(r[j] for r in profile)-min(r[j] for r in profile) for j in [5,6,7]),
        'min_accepted_convective_mass_fraction':min(r[11] for r in track['history'][1:]),
        'max_recorded_last_halfstep_luminosity_imbalance':max(abs(r[8]) for r in track['history']),
        'max_recorded_last_halfstep_nuclear_mass_imbalance':max(abs(r[9]) for r in track['history']),
        'plasmon_luminosity_over_surface':loss/profile[-1][4],'maximum_plasmon_epsilon_erg_g_s':maxloss,
        'nuclear_rest_mass_loss_over_baryonic_mass':rest/profile[-1][0],
        'plasmon_source':'https://arxiv.org/abs/astro-ph/9309014; equations 23--27, fully ionized approximation; does not include bremsstrahlung/photo/pair/recombination'
    }
    a.output.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
