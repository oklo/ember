"""Fill the pilot mass grid; these clocks still require initial-age calibration."""
from pathlib import Path
import json, re

W = Path(__file__).resolve().parent
cases = []
for cluster in ('pleiades', 'hyades'):
    template = W / f'tracks/cluster_{cluster}_nc_m0.150_shallow/configuration.txt'
    for mass, radius, temperature, loss in (
        (.045, .26, 2650, 2e-7), (.075, .33, 2850, 2e-7),
        (.090, .37, 2900, 2e-7), (.120, .40, 3050, 1.5e-7),
        (.200, .42, 3250, 7e-8),
    ):
        name = f'cluster_{cluster}_nc_m{mass:.3f}_grid_v1'
        d = W / 'tracks' / name
        if d.exists():
            raise FileExistsError(d)
        changes = dict(mass_Msun=str(mass), initial_radius_Rsun=str(radius),
                       initial_Teff_K=str(temperature), initial_entropy_loss=str(loss),
                       envelope_mass_fraction='.001', envelope_thermal_fraction_limit='.005')
        text = template.read_text()
        for key, value in changes.items():
            text, count = re.subn(rf'(?m)^{key} .*$', f'{key} "{value}"', text)
            assert count == 1, key
        d.mkdir()
        (d/'configuration.txt').write_text(text)
        (d/'scope.json').write_text(json.dumps(dict(
            cluster=cluster, mass_Msun=mass, cloud='clear',
            target_years=200000000 if cluster == 'pleiades' else 850000000,
            maximum_cpu_seconds=600,
            scope='Post-deuterium pilot. Relative clock; initial age not calibrated. '
                  'Do not interpret as a complete Hayashi-started sequence.',
        ), indent=2)+'\n')
        cases.append(name)
(W/'mass_grid_cases.json').write_text(json.dumps(cases, indent=2)+'\n')
print(' '.join(cases))
