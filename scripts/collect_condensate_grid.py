#!/usr/bin/env python3
"""Inventory accepted source cells for an explicitly specified rectangle.

Known cells in a plan and additional work roots are revalidated independently.
Missing cells and inconsistent duplicate solutions prevent completion. This
does not fill cells, change the requested axes, or install a runtime table.
"""
import argparse
import itertools
import json
from pathlib import Path

from validate_condensate_model import validate
from prepare_nongrey_sources import digest


def collect(prepared,specification,plan,roots):
    source=json.loads(Path(prepared).read_text());spec=json.loads(Path(specification).read_text())
    plan=json.loads(Path(plan).read_text())
    opacity_hashes={(r['XH'],r['X3']):digest(Path(r['directory'])/'opacity.bin') for r in plan['opacity_planes']}
    expected=set(itertools.product(spec['hydrogen'],spec['helium3'],spec['teff_K'],spec['log_g']))
    paths=[Path(r['directory']) for r in plan.get('models',[])]
    for root in roots:paths.extend(p.parent for p in sorted(Path(root).rglob('validated.json')))
    seen=set();models={};pending=[];rejected=[];duplicates=[];conflicts=[]
    for directory in paths:
        directory=directory.resolve()
        if directory in seen:continue
        seen.add(directory)
        if not (directory/'validated.json').exists():continue
        provenance=json.loads((directory/'provenance.json').read_text())
        key=tuple(provenance[k] for k in ['XH','X3','teff_K','log_g'])
        if key not in expected:continue
        if not (directory/'chemistry-audit.json').exists():
            pending.append({'coordinates':key,'directory':str(directory)});continue
        try:result=validate(directory,spec,source,key,opacity_hashes[key[:2]])
        except (ValueError,FileNotFoundError) as error:
            rejected.append({'coordinates':key,'directory':str(directory),'reason':str(error)});continue
        row={**dict(zip(['XH','X3','teff_K','log_g'],key)),'directory':str(directory),
             'diagnostics':result['diagnostics'],'chemistry':result['chemistry']}
        if key in models:
            previous=models[key]
            differences={k:row['diagnostics'][k]/previous['diagnostics'][k]-1
                         for k in ['T','Pgas','source_density']}
            check={'coordinates':key,'selected':previous['directory'],'other':str(directory),
                   'relative_differences':differences}
            duplicates.append(check)
            # A consistency gate for alternate numerical solutions, not a
            # physical uncertainty estimate or relaxed atmosphere tolerance.
            if max(abs(v) for v in differences.values())>1e-5:conflicts.append(check)
        else:models[key]=row
    missing=sorted(expected-set(models))
    return {'status':'complete' if not missing and not conflicts else 'incomplete',
            'expected_models':len(expected),'accepted_models':len(models),
            'specification':str(specification),'missing_coordinates':missing,
            'awaiting_chemistry':pending,'rejected_candidates':rejected,
            'duplicate_consistency_relative_tolerance':1e-5,
            'duplicate_comparisons':duplicates,'conflicting_solutions':conflicts,
            'models':[models[k] for k in sorted(models)]}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['prepared','specification','plan','output']:p.add_argument(name,type=Path)
    p.add_argument('roots',nargs='+',type=Path);p.add_argument('--require-complete',action='store_true')
    a=p.parse_args();result=collect(a.prepared,a.specification,a.plan,a.roots)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    temporary=a.output.with_suffix('.pending');temporary.write_text(json.dumps(result,indent=2)+'\n');temporary.replace(a.output)
    if result['status']=='complete':
        a.output.with_suffix('.models.json').write_text(json.dumps(result['models'],indent=2)+'\n')
    else:a.output.with_suffix('.models.json').unlink(missing_ok=True)
    print(json.dumps({k:result[k] for k in ['status','expected_models','accepted_models','missing_coordinates']},indent=2))
    if a.require_complete and result['status']!='complete':raise SystemExit(1)
