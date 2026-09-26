#!/usr/bin/env python3
"""Extract retained plasma-ratio comparisons relevant to the opacity connection."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
source = root/'docs/results/tops_cool_ratio_hydrogen_comparison_v1.json'
connection = root/'docs/results/opacity_connections_v1.json'
output = root/'docs/results/opacity_bridge_existing_ratios_v1.json'
if output.exists():
    raise FileExistsError(output)
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
inputs = {str(p):sha(p) for p in [Path(__file__),source,connection]}
original = json.loads(source.read_text())
rows = [r for r in original['comparisons']
        if .01 <= r['temperature_keV'] <= .03 and .01 <= r['density_atomic_g_cm3'] <= 1]
assert len(rows)==72
summaries = []
for x in sorted({r['X'] for r in rows}):
    group = [r for r in rows if r['X']==x]
    worst = max(group,key=lambda r:abs(r['relative_error']))
    summaries.append(dict(X=x,comparisons=len(group),maximum_relative_error=abs(worst['relative_error']),worst=worst))
report = dict(outcome='retained_connection_regime_controls_identified',accepted_for_stellar_opacity=False,
    candidate_X=original['candidate_X'],candidate_Z=original['candidate_Z'],
    temperature_interval_keV=[.01,.03],atomic_density_interval_g_cm3=[.01,1.],
    comparisons=rows,summaries=summaries,maximum_relative_error=max(s['maximum_relative_error'] for s in summaries),
    input_sha256=inputs,
    limitations=['These are selected comparisons from a previously completed audit, not new independent source calculations.',
        'The density/temperature interval brackets the saved gap; full interpolation support and a smooth connection have not been checked.',
        'The comparisons test plasma-ratio interpolation at X=0.1,0.3,0.7 and Z=0.02, not ordinary opacity or the actual stellar X/Z.',
        'Using MESA mean opacity still requires resolving its plasma convention; the retained TOPS uncut grey sources should also be inspected.',
        'No opacity, stellar input, atmosphere or trajectory is selected or modified.'])
for p,h in inputs.items():
    assert sha(Path(p))==h,p
output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps(dict(report=str(output),sha256=sha(output),comparisons=len(rows),maximum_relative_error=report['maximum_relative_error'])))
