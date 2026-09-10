#!/usr/bin/env python3
"""Compare a complete depleted atmosphere grid with independent source models."""
import argparse
import json
from pathlib import Path
import tempfile

from audit_nongrey_family import interpolate
from import_nongrey_grid import import_grid,read_text,CONDENSATE_CALCULATION
from prepare_nongrey_sources import digest
from validate_condensate_model import validate,original,PHYSICAL_KEYS


def audit(manifest,opacity_archive,heldouts,output,maximum_boundary_error=.005):
    manifest=Path(manifest);spec=json.loads(manifest.read_text())
    if spec['calculation']!=CONDENSATE_CALCULATION:raise ValueError('a condensate grid is required')
    opacity_archive=Path(opacity_archive)
    receipt=json.loads((opacity_archive/'archive.json').read_text())
    for name,expected in receipt['files_sha256'].items():
        if digest(opacity_archive/name)!=expected:raise ValueError('independent opacity archive changed')
    op=json.loads(read_text(original(opacity_archive,'provenance.json')))
    if (op['mode']!='equilibrium' or op['prepared']!=spec['provenance'] or
            any(op['specification'].get(k)!=spec.get(k) for k in PHYSICAL_KEYS)):
        raise ValueError('independent opacity physical inputs differ')
    with tempfile.TemporaryDirectory() as temp:
        states=import_grid(manifest,Path(temp)/'grid.dat')
    comparisons=[]
    for directory in heldouts:
        direct=validate(directory,spec,spec['provenance'],opacity_sha256=receipt['opacity_sha256'])
        point=direct['coordinates']
        if point[:2]!=(op['XH'],op['X3']):raise ValueError('heldout opacity composition differs')
        if point in states:raise ValueError('heldout must be independent of the grid nodes')
        estimate,corners=interpolate(spec,states,point)
        relative={k:estimate[k]/direct['diagnostics'][k]-1 for k in estimate}
        comparisons.append({'coordinates':point,'source':str(directory),'corners':corners,
            'direct':{k:direct['diagnostics'][k] for k in estimate},'interpolated':estimate,
            'relative_difference':relative,'source_receipt_sha256':digest(original(directory,'completed.json'))})
    if not comparisons:raise ValueError('at least one independent source model is required')
    error=max(abs(r['relative_difference'][k]) for r in comparisons for k in ['T','Pgas'])
    report={'description':'Independent condensate interpolation checks; local sensitivities, not global stellar error bounds',
            'manifest_sha256':digest(manifest),'opacity_archive_receipt_sha256':digest(opacity_archive/'archive.json'),
            'comparisons':comparisons,'maximum_relative_boundary_error':error,
            'maximum_allowed_relative_boundary_error':maximum_boundary_error,'passed':error<=maximum_boundary_error,
            'density_note':'Source-density interpolation is diagnostic. Ember recomputes density with its own EOS.'}
    Path(output).write_text(json.dumps(report,indent=2)+'\n')
    if not report['passed']:raise ValueError('atmosphere interpolation requires refinement; report retained')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['manifest','opacity_archive','output']:p.add_argument(name,type=Path)
    p.add_argument('heldouts',nargs='+',type=Path)
    p.add_argument('--maximum-boundary-error',type=float,default=.005)
    a=p.parse_args()
    if not 0<a.maximum_boundary_error<1:raise ValueError('invalid comparison criterion')
    result=audit(a.manifest,a.opacity_archive,a.heldouts,a.output,a.maximum_boundary_error)
    print(json.dumps({'passed':result['passed'],'maximum_relative_boundary_error':result['maximum_relative_boundary_error']}))
