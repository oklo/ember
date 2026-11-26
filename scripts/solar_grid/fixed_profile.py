"""Recompute the evolution boundary and fluxes of the saved 850 K iterate at its own optical depths."""
import math
import os
from pathlib import Path
import sys
from collections import namedtuple

from cooling_interval import MATCH_LIMITS, digest, read, save

root=Path(sys.argv[1]).resolve()
ws=root.parent/'picaso-pilot-oct11-fable'
physics=root.parent/'highz-native-opacity-oct12-v1'
sys.path[:0]=list(map(str,[physics/'normalized_native',physics/'column_gravity',physics/'warm_atmosphere_controls',
                ws/'hhe_pos_p2671',ws/'multiregion_p2678',ws,root/'thermal-holdout/m06']))
job=root/'recovery/m06-850-birth-depth10-column'
receipt_path=Path(read(job/'selected-column.json')['receipt'])
native=read(receipt_path)
assert native['quality']['usable']
spec=native['spec']
mix=spec['arms']['anchor']
he3=mix['He3']/3.01603
helium=he3+mix['He4']/4.002602
he3_share=he3/helium
os.environ.update(PICASO_EFIX='3',GS98_PER_ELEMENT='1',HEH_OVERRIDE=repr(helium/(mix['X']/1.00794)),
                  ALLOW_NO100='1',picaso_refdata=str(ws/'picaso_src/reference'),PYSYN_CDBS=str(ws/'pysyn_empty'),
                  NUMBA_NUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
os.environ.pop('ROUTE',None)
import numpy as np
import pandas as pd
import astropy.units as u
import endpoint_gravity
import run_hhe
import picaso.climate as C
import picaso.justdoit as jdi
from picaso.photochem import EquilibriumChemistry
from picaso.justdoit import LODDERS2020_C_TO_O
from normalization import install as normalize
from tau_diag_p2535 import weights, means
from multi_region_p2678 import audit
from archive_matching_points import matching_points

endpoint_gravity.install()
normalize(run_hhe.jdi,add_lithium=True)
run_hhe.GASES=list(dict.fromkeys([*run_hhe.GASES,'Li']))
r=read(Path(native['column']))
path=job/'column-01'/('latest_guess_anchor'+spec['tag']+'.json')
guess=read(path)
assert guess['accepted'] is False
regions=guess['solved_regions']
P,T=np.asarray(guess['pressure_bar']),np.asarray(guess['temperature_K'])
assert len(P)==len(T)==spec['levels'] and np.all(np.diff(P)>0)
grav=r['grav_mps2']
assert math.isclose(grav,10**spec['logg']/100,rel_tol=1e-10)
opa = jdi.opannection(method='resortrebin', preload_gases=run_hhe.GASES)
run = jdi.inputs(calculation='browndwarf', climate=True); run.gravity(gravity=grav, gravity_unit=u.Unit('m/(s**2)')); run.effective_temp(float(r['teff']))
run.inputs_climate(temp_guess=T, pressure=P, rcb_guess=60, rfacv=0.0)
s = EquilibriumChemistry(thermofile=os.path.join(os.environ['picaso_refdata'], 'chemistry', 'thermo_data', 'thermo-sonora-component.yaml'), method='sonora-approx')
chem_calls=[]
original_solve=s.solve_metallicity
def trace_solve(*args,**kwargs):
    ok=original_solve(*args,**kwargs)
    chem_calls.append(bool(ok))
    return ok
s.solve_metallicity=trace_solve
s.set_atom_to_H_ratios(dict(r['solver_atom_to_H_ratios']))
gas, _ = s.equilibrate_atmosphere(P, T, float(np.log10(r['mh_linear'])), r['cto_absolute'] / LODDERS2020_C_TO_O)
# Each layer may retry up to nine times in this campaign's chemistry source.
cursor=0
failed_layers=[]
for i in range(len(P)):
    for attempt in range(9):
        ok=chem_calls[cursor]
        cursor+=1
        if ok:break
    if not ok:failed_layers.append(i)
assert cursor==len(chem_calls)
mu_chem = np.empty(len(P)); s.use_prev_guess = False
point_chemistry=[]
for i in range(len(P)):
    point_chemistry.append(s.solve_metallicity(P[i]*1e6,T[i],r['mh_linear'],r['cto_absolute']/LODDERS2020_C_TO_O))
    mu_chem[i]=s.mubar
df = pd.DataFrame({'pressure': P, 'temperature': T}); [df.__setitem__(k, np.maximum(v, 0.0)) for k, v in gas.items()]; run.atmosphere(df=df)
nsum = np.sum([np.maximum(v, 0.0) for v in gas.values()], axis=0); xHe = np.asarray(gas['He'])
mu_true = mu_chem - xHe * he3_share * (4.002602 - 3.01603)
OW, ON, Sc, Di, Atm, ho = C.calculate_atm(run, opa); D = np.asarray(ON.DTAU, float); D = D if D.ndim == 3 else D[:, :, None]
w, dw, gw = np.asarray(opa.wno, float), np.asarray(opa.delta_wno, float), np.asarray(opa.gauss_wts, float)
Bw, dBw = weights(w, dw, np.sqrt(T[1:] * T[:-1])); dR, _, bad, zero = means(D, gw, Bw, dBw); tau = np.concatenate([[0.0], np.cumsum(dR)])
tau_bottom_formal = float(tau[-1]); tau_bottom_trend = float(tau[-2] + dR[-2] * dR[-2] / dR[-3])   # P2679: endpoint layers carry 2x column mass (g/2 bug); trend-corrected last-layer dtau
Og = namedtuple('Opagrid', ['nwno', 'delta_wno', 'wno', 'ngauss', 'gauss_wts', 'tmin', 'tmax'])(opa.nwno, opa.delta_wno, opa.wno, opa.ngauss, opa.gauss_wts, min(opa.temps) * 0.7, 10000)
fl = C.get_fluxes(Atm, OW, ON, Sc, Di, Og, np.ones(opa.nwno), False, True); Fr = np.asarray(fl[5], float) / (5.670374419e-5 * r['teff'] ** 4)
n=len(P)
STRICT=np.zeros(n,bool)
STRICT[audit(regions,n)]=True
interface=np.zeros(n,bool)
for top,bot in regions:
    interface[top]=True
for top,bot in regions[:-1]:
    interface[bot+2]=True
lP, lT = np.log(P), np.log(T)
ci = run.inputs['climate']
AB = namedtuple('AB', ['t_table', 'p_table', 'grad', 'cp'])(ci['t_table'], ci['p_table'], ci['grad'], ci['cp'])
gad = np.array([C.did_grad_cp(.5*(T[j]+T[j+1]), math.sqrt(P[j]*P[j+1]), AB)[0]
                for j in range(len(T)-1)])
grad_ratio = np.diff(lT)/np.diff(lP)/gad
rad=~(STRICT|interface)
checks=dict(radiative_flux=float(max(abs(Fr[rad]-1))),interface_error=float(max(abs(Fr[interface]-1))),
            inward_flux=float(max(0,max(Fr[STRICT|interface])-1)),radiative_gradient=float(max(grad_ratio[~STRICT[1:]])),
            surface_flux=float(Fr[0]),bottom_tau=float(tau[-1]),chemistry_point_failures=sum(not v for v in point_chemistry),
            opacity_in_domain=bool(min(T)>=min(opa.temps) and max(T)<=max(opa.temps) and min(P)>=min(opa.pressures) and max(P)<=max(opa.pressures) and not np.any(bad)))
profile=dict(P=P.tolist(),T=T.tolist(),mu=mu_true.tolist(),tau=tau.tolist())
early=matching_points(dict(analysis=dict(profile=profile),fresh=dict(at_tau={}),logg=spec['logg']))['10']
final=matching_points(native)['10']
keys=['T','P_bar','depth_cm_mu_He3corr']
ratios=[early[k]/final[k] for k in keys]
save(root/'atmosphere-validation/early850/fixed-profile.json',dict(early_over_final=ratios,
     early=early,final=final,checks=checks,matching_within_limits=all(abs(v-1)<=lim for v,lim in zip(ratios,MATCH_LIMITS,strict=True)),
     physical_checks_passed=bool(checks['radiative_flux']<=.01 and checks['interface_error']<=.02 and
     checks['inward_flux']<=.01 and checks['radiative_gradient']<=1.02 and abs(checks['surface_flux']-1)<=.01 and
     min(tau_bottom_formal,tau_bottom_trend)>=30 and checks['chemistry_point_failures']==0 and not failed_layers and
     checks['opacity_in_domain'] and np.all(np.diff(tau)>0) and np.all(np.isfinite(mu_true)) and np.min(mu_true)>0),
     failed_chemistry_layers=failed_layers,iterate_sha256=digest(path),accepted_receipt_sha256=digest(receipt_path),production_eligible=False,
     scope='Retrospective fixed-profile calibration at recomputed optical depth. No solver convergence flag or production atmosphere awarded.'))
