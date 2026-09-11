#!/usr/bin/env python3
"""Extend hot TOPS density coverage while retaining every original source cell.

Version-2 tables give the valid density prefix of each isotherm. All source
substitutions are excluded. Runtime additionally requires the full value and
derivative stencil. Assembly alone does not accept this candidate for a star.
"""
import argparse
import json
import math
from pathlib import Path
import shlex
import shutil

from assemble_hot_opacity_extension import planes
from import_tops_composition import read
from prepare_nongrey_sources import digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['previous', 'density_manifest', 'output']:
        p.add_argument(name, type=Path)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError('use a fresh candidate directory')
    manifest = json.loads(a.density_manifest.read_text())
    if manifest.get('complete') is not True or manifest.get('accepted_for_stellar_opacity') is not False:
        raise ValueError('a complete raw-source density calculation is required')
    sources, inputs = {}, {str(a.density_manifest.resolve()): digest(a.density_manifest),
                          str(Path(__file__).resolve()): digest(__file__)}
    for record in manifest['requests']:
        key = (record['X'], record['Z'])
        if key in sources:
            raise ValueError('duplicate density source mixture')
        work = Path(record['work'])
        coverage = work/'coverage.json'
        if digest(coverage) != record['coverage_sha256']:
            raise ValueError('source coverage report changed')
        report = json.loads(coverage.read_text())
        receipt_path = work/'receipt.json'
        receipt = json.loads(receipt_path.read_text())
        request = work/'request.json'
        raw = work/'source.txt'
        dimensions = tuple(report['dimensions'])
        if (len(dimensions) != 2 or dimensions[0] != 50 or not 2 <= dimensions[1] <= 100
                or tuple(receipt['dimensions']) != dimensions
                or (receipt['X'], receipt['Z']) != key
                or digest(request) != receipt['request_sha256']
                or digest(raw) != report['source_sha256']):
            raise ValueError('density source request or mixture differs')
        grid = read(raw, receipt, dimensions=dimensions)
        sources[key] = grid
        for path in [coverage, receipt_path, request, raw]:
            inputs[str(path.resolve())] = digest(path)
    family = a.previous/'tops_gs98_mixture_high.dat'
    replacements, summaries, used = {}, [], set()
    for line in family.read_text().splitlines()[1:]:
        z, name = shlex.split(line)
        z = float(z)
        if Path(name).name != name:
            raise ValueError('expected high table inside its source directory')
        path = a.previous/name
        (nt, nr, axis_lines), xs, rows = planes(path)
        old_density, old_temperature = [list(map(float, text.split())) for text in axis_lines]
        first = sources[(xs[0], z)]
        tt = [t for t in first[0] if t >= .025]
        new_density = first[1]
        # Preserve the original importer's arithmetic order as well as its
        # constants, including the last-bit kelvin conversion at 1.875 keV.
        if ([math.log10(t*1e3*1.602176634e-12/1.380649e-16) for t in tt] != old_temperature
                or math.log10(new_density[0]) != old_density[-1]):
            raise ValueError('density source does not share the old temperature/density overlap')
        log_density = old_density+[math.log10(r) for r in new_density[1:]]
        output = [f'EMBER_OPACITY_TABLE 2 {len(xs)} {nt} {len(log_density)} TOPS ATOMIC GS98; original density prefixes only',
                  ' '.join(format(v,'.17g') for v in log_density), axis_lines[1]]
        retained, added, counts = 0, 0, []
        for x, old_rows in zip(xs, rows, strict=True):
            key = (x,z); used.add(key)
            full_tt, rr, cells, excluded = sources[key]
            if full_tt != first[0] or rr != new_density:
                raise ValueError('unaligned new source grids')
            output.append(old_rows[0])
            support = []
            for t, row in zip(tt, old_rows[1:], strict=True):
                old = list(map(float,row.split()))
                if len(old) != nr or (t,rr[0]) in excluded or math.log10(cells[t,rr[0]]) != old[-1]:
                    raise ValueError('original hot-opacity cell changed or was substituted')
                extra = []
                missing = False
                for r in rr[1:]:
                    if (t,r) in excluded:
                        missing = True
                    elif missing:
                        raise ValueError('source support is not a contiguous density prefix')
                    else:
                        extra.append(math.log10(cells[t,r]))
                count = nr+len(extra)
                # Original tokens are retained; no resampling or replacement
                # of the original opacity values occurs during assembly.
                output.append(str(count)+' '+row+(' '+' '.join(format(v,'.17g') for v in extra) if extra else ''))
                retained += nr; added += len(extra); support.append(count)
            counts.append({'X':x,'density_prefix_sizes':support})
        replacements[name] = '\n'.join(output)+'\n'
        summaries.append({'Z':z,'file':name,'retained_source_cells':retained,
                          'added_source_cells':added,'support':counts})
    if used != set(sources):
        raise ValueError('density source mixtures differ from the selected hot family')
    # Make the copied family self-contained. The selected AESOPUS manifest
    # may refer to its parent's source files; copying that relative path to
    # an unrelated candidate directory would lose those physical tables.
    external = {}
    for name in ['aesopus21_gs98_mixture.dat', 'tops_gs98_mixture_low.dat',
                 'tops_gs98_mixture_high.dat']:
        path = a.previous/name
        lines = path.read_text().splitlines()
        header = lines[0].split()
        if (len(header) != 5 or header[:2] != ['EMBER_OPACITY_MIXTURE', '1']
                or len(lines)-1 != int(header[2])):
            raise ValueError('invalid source opacity manifest')
        rewritten = [lines[0]]
        changed = False
        for line in lines[1:]:
            z, relative = shlex.split(line)
            source = (a.previous/relative).resolve(strict=True)
            inputs[str(source)] = digest(source)
            if Path(relative).name != relative:
                target = source.name
                if ((a.previous/target).exists() or target in replacements
                        or target in external and external[target] != source):
                    raise ValueError('external opacity source filename collision')
                external[target] = source
                rewritten.append(z+' '+json.dumps(target))
                changed = True
            else:
                rewritten.append(line)
        if changed:
            replacements[name] = '\n'.join(rewritten)+'\n'
    for path in a.previous.rglob('*'):
        if path.is_file():
            inputs[str(path.resolve())] = digest(path)
    if any(digest(path) != checksum for path, checksum in inputs.items()):
        raise ValueError('source input changed during assembly')
    shutil.copytree(a.previous, a.output)
    for name, source in external.items():
        shutil.copy2(source, a.output/name)
    for name, text in replacements.items():
        (a.output/name).write_text(text)
    report = {'scope': __doc__, 'accepted_for_stellar_opacity':False,
              'input_sha256':inputs,'tables':summaries,
              'output_sha256':{p.name:digest(p) for p in a.output.iterdir() if p.is_file()}}
    (a.output/'density_extension_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'retained_cells':sum(r['retained_source_cells'] for r in summaries),
                      'added_cells':sum(r['added_source_cells'] for r in summaries)}))


if __name__=='__main__':
    main()
