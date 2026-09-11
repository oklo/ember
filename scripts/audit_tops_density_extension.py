#!/usr/bin/env python3
"""Check dense-opacity source nodes, independent mixtures/densities and derivatives.

Only un-substituted source cells and complete runtime stencils enter source
comparisons. The independent checks start where the stellar opacity uses the
hot TOPS table alone. Results are local interpolation checks, not a physical
opacity uncertainty or a stellar-trajectory validation.
"""
import argparse
import gzip
import json
import math
from pathlib import Path
import subprocess

from import_tops_composition import read
from prepare_nongrey_sources import digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['family', 'nodes', 'heldouts', 'probe', 'output']:
        p.add_argument(name, type=Path)
    p.add_argument('--tolerance', type=float, default=.005,
                   help='relative opacity criterion for the independent sampled states')
    a = p.parse_args()
    if not 0 < a.tolerance < .1:
        raise ValueError('invalid comparison criterion')
    inputs = {str(f.resolve()): digest(f) for f in
              [a.nodes, a.heldouts, a.probe, Path(__file__), Path(__file__).with_name('import_tops_composition.py')]}
    for f in a.family.glob('*.dat'):
        inputs[str(f.resolve())] = digest(f)
    raw = []
    def query(z, points):
        table = a.family/f'tops_gs98_mixture_z{round(z*1000):03d}_high.dat'
        request = ''.join(' '.join(format(v,'.17g') for v in q)+'\n' for q in points)
        result = subprocess.run([str(a.probe.resolve()),str(table.resolve())],
                                input=request,text=True,capture_output=True,check=True)
        rows = [json.loads(line) for line in result.stdout.splitlines()]
        if len(rows) != len(points) or any(r['query'] != q for q,r in zip(points,rows,strict=True)):
            raise ValueError('probe changed or omitted queries')
        for row in rows:
            if type(row.get('covered')) is not bool:
                raise ValueError('missing coverage flag')
            if row['covered'] and (not row['kappa'] > 0 or any(not math.isfinite(row[k])
                    for k in ['kappa','dlnk_dlnT','dlnk_dlnrho','dlnk_dX'])):
                raise ValueError('invalid covered opacity response')
        raw.append({'Z':z,'input':request,'stdout':result.stdout})
        return rows
    summaries, derivative_points = [], []
    density_axes = {}
    for kind, path in [('nodes',a.nodes),('heldouts',a.heldouts)]:
        manifest = json.loads(path.read_text())
        if manifest.get('complete') is not True:
            raise ValueError('source calculations are unfinished')
        for entry in manifest['requests']:
            work = Path(entry['work'])
            receipt = json.loads((work/'receipt.json').read_text())
            coverage = json.loads((work/'coverage.json').read_text())
            dimensions = tuple(receipt['dimensions'])
            if (digest(work/'coverage.json') != entry['coverage_sha256']
                    or digest(work/'request.json') != receipt['request_sha256']
                    or (entry['X'],entry['Z']) != (receipt['X'],receipt['Z'])
                    or dimensions != tuple(coverage['dimensions'])
                    or len(dimensions) != 2 or dimensions[0] != 50
                    or not 2 <= dimensions[1] <= 100):
                raise ValueError('source request or coverage changed')
            tt,rr,cells,excluded = read(work/'source.txt',receipt,dimensions=dimensions)
            if kind == 'nodes':
                if entry['Z'] in density_axes and density_axes[entry['Z']] != rr:
                    raise ValueError('unaligned source-node density grids')
                density_axes[entry['Z']] = rr
            for f in [work/'receipt.json',work/'request.json',work/'coverage.json',work/'source.txt']:
                inputs[str(f.resolve())] = digest(f)
            material = [(t,r) for t in tt if math.log10(t*1e3*1.602176634e-12/1.380649e-16) >= 5.7 for r in rr]
            points = [[entry['X'],t*1e3*1.602176634e-12/1.380649e-16,r] for t,r in material]
            responses = query(entry['Z'],points)
            compared, unsupported = [], 0
            for (t,r),point,row in zip(material,points,responses,strict=True):
                if not row['covered']:
                    unsupported += 1
                    continue
                if (t,r) in excluded:
                    raise ValueError('runtime admitted a substituted source state')
                compared.append({'query':point,'source':cells[t,r],
                    'relative_difference':row['kappa']/cells[t,r]-1,
                    'intermediate_density': kind == 'heldouts' and r not in density_axes[entry['Z']]})
            if not compared:
                raise ValueError('source mixture has no supported comparisons')
            worst = max(compared,key=lambda r:abs(r['relative_difference']))
            summaries.append({'kind':kind,'X':entry['X'],'Z':entry['Z'],
                'compared_states':len(compared),'unsupported_runtime_states':unsupported,'worst':worst,
                'worst_at_existing_density':max((r for r in compared if not r['intermediate_density']),
                    key=lambda r:abs(r['relative_difference'])),
                'worst_at_intermediate_density':max((r for r in compared if r['intermediate_density']),
                    key=lambda r:abs(r['relative_difference']),default=None)})
            if kind == 'heldouts':
                derivative_points += [(entry['Z'],r['query']) for r in compared[::max(1,len(compared)//4)][:4]]
    maximum_derivative, derivative_count = 0., 0
    h = 1e-6
    for z, point in derivative_points:
        points = [point]
        for axis in range(3):
            for sign in [1,-1]:
                q = point.copy();q[axis] = q[axis]+sign*h if axis==0 else q[axis]*math.exp(sign*h)
                points.append(q)
        rows = query(z,points)
        if not all(r['covered'] for r in rows):
            continue
        for i,key in enumerate(['dlnk_dX','dlnk_dlnT','dlnk_dlnrho']):
            fd = math.log(rows[1+2*i]['kappa']/rows[2+2*i]['kappa'])/(2*h)
            maximum_derivative = max(maximum_derivative,abs(fd-rows[0][key])/max(1.,abs(fd)))
        derivative_count += 1
    source_max = max(abs(r['worst']['relative_difference']) for r in summaries if r['kind']=='nodes')
    independent_max = max(abs(r['worst']['relative_difference']) for r in summaries if r['kind']=='heldouts')
    passed = source_max < 1e-11 and independent_max <= a.tolerance and maximum_derivative < 3e-5 and derivative_count > 0
    if any(digest(f) != checksum for f,checksum in inputs.items()):
        raise ValueError('an input changed during the audit')
    raw_path = a.output.with_suffix('.probe.json.gz')
    raw_path.write_bytes(gzip.compress(json.dumps(raw,allow_nan=False).encode(),mtime=0))
    report = {'scope':__doc__,'passed':passed,'input_sha256':inputs,'comparisons':summaries,
              'maximum_source_node_relative_difference':source_max,
              'maximum_independent_relative_difference':independent_max,'independent_relative_tolerance':a.tolerance,
              'derivative_states':derivative_count,'derivative_boundary_states_skipped':len(derivative_points)-derivative_count,
              'maximum_derivative_difference':maximum_derivative,'raw_probe_sha256':digest(raw_path)}
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['scope','comparisons','input_sha256']}))
    if not passed:
        raise SystemExit('dense opacity failed its specified source comparison')


if __name__=='__main__':
    main()
