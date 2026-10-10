"""Inventory actual cluster runs without mistaking elapsed age for cluster age."""
from pathlib import Path
import json, shlex, math
import numpy as np

W = Path(__file__).resolve().parent
rows = []
for d in sorted((W/'tracks').iterdir()):
    if not (d/'configuration.txt').exists():
        continue
    cfg = dict(shlex.split(line) for line in (d/'configuration.txt').read_text().splitlines()
               if line.strip() and not line.startswith('#'))
    scope = json.loads((d/'scope.json').read_text()) if (d/'scope.json').exists() else {}
    done = json.loads((d/'completion.json').read_text()) if (d/'completion.json').exists() else None
    history = []
    if (d/'evolution/history.jsonl').exists():
        lines = (d/'evolution/history.jsonl').read_text().splitlines()
        for i, line in enumerate(lines):
            try:
                history.append(json.loads(line))
            except json.JSONDecodeError:
                if i != len(lines)-1:
                    raise
    cluster = scope.get('cluster') or next((c for c in ('hyades','pleiades') if c in d.name), None)
    row = dict(case=d.name, cluster=cluster, mass_Msun=float(cfg['mass_Msun']),
               cloud=scope.get('cloud', 'nc'),
               scope=scope.get('scope', 'Pilot; see original configuration and receipts'),
               status='complete' if done and done['returncode']==0 else 'failed' if done
                      else 'running_or_interrupted' if (d/'launch.json').exists() else 'prepared',
               target_elapsed_years=scope.get('target_years', 800000000),
               post_deuterium_start='post_deuterium' in cfg['composition'],
               elapsed_age_is_not_cluster_age=True,
               envelope_mass_fraction=float(cfg.get('envelope_mass_fraction',0)),
               accepted_records=len(history),
               wall_seconds=done['wall_seconds'] if done else None)
    if history:
        first, last = history[0], history[-1]
        row.update(start_radius_Rsun=first['radius_Rsun'], start_Teff_K=first['Teff_K'],
                   final_elapsed_years=last['years'], final_Teff_K=last['Teff_K'],
                   final_luminosity_Lsun=last['luminosity_Lsun'],
                   accepted_steps=last['accepted'], rejected_steps=last['rejected'],
                   completed_target=last['years'] >= row['target_elapsed_years']*(1-1e-12))
        # n=1.5 constant-Teff contraction age scale, not an assigned age correction.
        row['initial_contraction_scale_years']=(6.67430e-8*(row['mass_Msun']*1.98847e33)**2 /
            (7*first['radius_Rsun']*6.957e10*first['luminosity_Lsun']*3.828e33)/31557600)
        elapsed = np.array([r['years'] for r in history])
        row['at_elapsed_years'] = {}
        for age in (100000000,125000000,150000000,600000000,650000000,750000000):
            if age > elapsed[-1]:
                continue
            row['at_elapsed_years'][str(age)] = {key:float(np.interp(age,elapsed,[r[key] for r in history]))
                for key in ('luminosity_Lsun','radius_Rsun','Teff_K','central_X')}
    rows.append(row)
(W/'track_inventory.json').write_text(json.dumps(rows,indent=2)+'\n')
successful=[r for r in rows if r['status']=='complete' and r.get('completed_target') and r['cluster']]
print(json.dumps(dict(successful_cluster_runs=len(successful),
    running=[r['case'] for r in rows if r['status']=='running_or_interrupted'],
    prepared=[r['case'] for r in rows if r['status']=='prepared'],
    warning='Relative model clocks remain uncalibrated.'),indent=2))
