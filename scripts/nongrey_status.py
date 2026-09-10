#!/usr/bin/env python3
"""Report saved source progress without inspecting or signalling processes."""
import argparse
import itertools
import json
import math
from pathlib import Path
import re


def iteration_progress(directory):
    """Read only complete depth sweeps; a partial sweep understates errors."""
    directory=Path(directory);parameters=directory/'tas';path=directory/'fort.9'
    if not parameters.exists() or not path.exists():return None
    match=re.search(r'\bND\s*=\s*(\d+)',parameters.read_text(),re.I)
    if not match:return None
    depths=int(match[1]);iterations={}
    for line in path.read_text().splitlines():
        words=line.replace('D','E').split()
        if len(words)!=9:continue
        try:iteration,depth=int(words[0]),int(words[1]);correction=float(words[2])
        except ValueError:continue
        if 1<=depth<=depths and math.isfinite(correction):
            iterations.setdefault(iteration,{})[depth]=correction
    complete=[]
    for iteration,values in sorted(iterations.items()):
        if len(values)!=depths:continue
        depth=max(values,key=lambda i:abs(values[i]))
        complete.append({'iteration':iteration,'maximum_relative_temperature_correction':abs(values[depth]),
                         'depth_of_maximum':depth})
    return {'completed_iterations':len(complete),'last_complete_sweeps':complete[-4:]}


def status(root):
    root=Path(root);spec=json.loads((root/"specification.json").read_text())
    rows=[]
    for i,(x,y) in enumerate(itertools.product(spec["hydrogen"],spec["helium3"])):
        p=root/f"plane-{i:03d}"
        active=[]
        for path in sorted(p.glob("**/running.json")):
            value=json.loads(path.read_text())
            progress=iteration_progress(path.parent)
            active.append({"directory":str(path.parent.relative_to(root)),**value,
                           **({'convergence_progress':progress} if progress is not None else {})})
        rows.append({"XH":x,"X3":y,
                     "saved_opacity_isotherms":len(list(p.glob("opacity/temperature-*/completed.json"))),
                     "required_opacity_isotherms":spec["log_temperature"][0],
                     "validated_atmospheres":len(list(p.glob("model-*/validated.json"))),
                     "rejected_atmospheres":sum(not (f.parent/'validated.json').exists() for f in p.glob('model-*/failure.json')),
                     "preserved_failure_receipts":len(list(p.glob('model-*/failure.json'))),
                     "required_atmospheres":len(spec["teff_K"])*len(spec["log_g"]),
                     "last_started_runs":active})
    return {"work":str(root.resolve()),"grid_imported":(root/"atmosphere.dat").exists(),
            "note":"Saved receipts are progress indicators; the importer revalidates physics and hashes. A running.json file can remain after interruption.",
            "planes":rows}


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("work",type=Path)
    print(json.dumps(status(p.parse_args().work),indent=2))
