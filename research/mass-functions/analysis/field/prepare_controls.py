"""Private field-age boundary controls, keeping the common code unchanged."""
from pathlib import Path
import re,json
W=Path(__file__).resolve().parent
(W/'physics').mkdir(exist_ok=True)
def prepare(parent,case,changes,scope_changes=None):
    src=W/'tracks'/parent;dst=W/'tracks'/case;dst.mkdir(exist_ok=True)
    if (dst/'launch.json').exists():return
    cfg=(src/'configuration.txt').read_text()
    atm=Path(re.search(r'^atmosphere_pressure_matched "([^"]+)"',cfg,re.M)[1])
    new=W/'physics'/atm.name
    if not new.exists():
        lines=atm.read_text().splitlines();tokens=lines[3].split()
        tokens[-2:]=['0.01','0.0003'];lines[3]=' '.join(tokens)
        new.write_text('\n'.join(lines)+'\n')
    changes={'atmosphere_pressure_matched':str(new),'atmosphere_maximum_delta_Z':'2e-5',
             'envelope_mass_fraction':'.0003','envelope_thermal_fraction_limit':'.01','fingering_oscillation':'reject',**changes}
    for k,v in changes.items():
        line=f'{k} "{v}"'
        cfg=re.sub(r'^'+re.escape(k)+r'\s+.*$',lambda _:line,cfg,flags=re.M) if re.search(r'^'+re.escape(k)+r'\s+',cfg,re.M) else cfg+line+'\n'
    (dst/'configuration.txt').write_text(cfg)
    scope=json.loads((src/'scope.json').read_text());scope.update(scope_changes or {})
    scope.update(parent_case=parent,private_changes=changes,
       approximation='Fixed initial-composition atmosphere permitted within absolute dX=.01 and dZ=.0003; native-envelope omitted thermal reservoir bounded at 1 percent of L. Explicit controls, not new composition-dependent atmosphere physics.')
    (dst/'scope.json').write_text(json.dumps(scope,indent=2)+'\n')
    return case
cases=[]
for cloud in ['f2','nc']:
    for mass in [.025,.035,.050,.060]:
        cases.append(prepare(f'field_{cloud}_m{mass:.3f}_v1',f'field_{cloud}_m{mass:.3f}_v5',{'fingering':'none'}))
    for mass in [.065,.070,.080,.085]:
        cases.append(prepare(f'field_{cloud}_m0.075_v1',f'field_{cloud}_m{mass:.3f}_v5',{'mass_Msun':str(mass),'fingering':'none'},{'mass_Msun':mass}))
    # Small-mass, young-age support below the main reported mass range.
    for mass,radius in [(.003,.10),(.005,.15)]:
        cases.append(prepare(f'field_{cloud}_m0.008_v1',f'field_{cloud}_m{mass:.3f}_v5',{'mass_Msun':str(mass),'initial_radius_Rsun':str(radius),'initial_Teff_K':'1500','fingering':'none'}, {'mass_Msun':mass,'initial_contraction_scale_years':0}))
    # Envelope depth controls: deliberately not another relaxation of its bound.
    cases.append(prepare(f'field_{cloud}_m0.025_v1',f'field_{cloud}_m0.025_thin_v5',{'envelope_mass_fraction':'.0001','fingering':'none'}))
    cases.append(prepare(f'field_{cloud}_m0.050_v1',f'field_{cloud}_m0.050_thin_v5',{'envelope_mass_fraction':'.0001','fingering':'none'}))
(W/'control_cases.json').write_text(json.dumps([c for c in cases if c],indent=2)+'\n')
print(len([c for c in cases if c]),'controls prepared')
