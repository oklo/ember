#!/usr/bin/env python3
"""Measure opacity sensitivity to omitting He I 4471 A in a separate source.

Omission is a diagnostic only. It is not a proposed physical prescription.
Original source code, opacity and atmosphere files remain unchanged. The
prepared diagnostic carries initialization_only so atmosphere continuation
and assembly reject it. Completed source isotherms remain reusable.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess

from generate_nongrey_grid import opacity_inputs, execute, sequence
from nongrey_opacity import validate_table
from prepare_nongrey_sources import digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('prepared', type=Path)
    p.add_argument('specification', type=Path)
    p.add_argument('source_directory', type=Path)
    p.add_argument('work', type=Path)
    p.add_argument('--temperature', type=float, nargs='+', required=True)
    p.add_argument('--hydrogen', type=float, required=True)
    p.add_argument('--helium3', type=float, required=True)
    a = p.parse_args()
    source = a.source_directory/'synspec54.f'
    prepared = json.loads(a.prepared.read_text())
    spec = json.loads(a.specification.read_text())
    if digest(prepared['synspec']) != prepared['executables']['synspec']:
        raise ValueError('original source executable changed')
    if a.source_directory.resolve() != Path(prepared['synspec']).resolve().parent:
        raise ValueError('source directory must accompany the selected executable')
    a.work.mkdir(parents=True, exist_ok=True)
    control = a.work/'source'
    control.mkdir(exist_ok=True)
    text = source.read_text()
    begin = text.index('      FUNCTION PHE1(ID,FREQ,ILINE)')
    location = text.index('      T=TEMP(ID)+2.42E-8*VTURB(ID)', begin)
    insertion = ('C     Diagnostic omission only; never an atmosphere source.\n'
                 '      IF(ILINE.EQ.1) THEN\n'
                 '         PHE1=0.D0\n'
                 '         RETURN\n'
                 '      END IF\n')
    modified = text[:location]+insertion+text[location:]
    inputs = {str(path.resolve()): digest(path) for path in
              [a.prepared, a.specification, source, Path(__file__),
               *sorted(a.source_directory.glob('*.FOR'))]}
    recipe = {'scope': __doc__, 'input_sha256': inputs, 'insertion': insertion,
              'temperature_K': a.temperature, 'XH': a.hydrogen, 'X3': a.helium3}
    saved = a.work/'recipe.json'
    if saved.exists() and json.loads(saved.read_text()) != recipe:
        raise ValueError('changed diagnostic inputs require a new work directory')
    saved.write_text(json.dumps(recipe, indent=2)+'\n')
    for path in a.source_directory.glob('*.FOR'):
        shutil.copy2(path, control/path.name)
    (control/'synspec54.f').write_text(modified)
    executable = control/'synspec54'
    build = a.work/'build.json'
    if build.exists():
        record = json.loads(build.read_text())
        if (record['source_sha256'] != digest(control/'synspec54.f')
                or record['executable_sha256'] != digest(executable)):
            raise ValueError('diagnostic build changed')
    else:
        command = ['gfortran', '-O2', '-g', '-fno-automatic', '-std=legacy',
                   '-fallow-argument-mismatch', '-fcheck=bounds', '-fbacktrace',
                   '-o', str(executable.resolve()), 'synspec54.f']
        with (a.work/'build.log').open('w') as log:
            subprocess.run(command, cwd=control, stdout=log, stderr=subprocess.STDOUT, check=True)
        record = {'command': command, 'source_sha256': digest(control/'synspec54.f'),
                  'executable_sha256': digest(executable)}
        build.write_text(json.dumps(record, indent=2)+'\n')
    prepared['synspec'] = str(executable.resolve())
    prepared['executables']['synspec'] = digest(executable)
    prepared['initialization_only'] = True
    prepared['helium_line_omission_diagnostic'] = recipe
    (a.work/'prepared.json').write_text(json.dumps(prepared, indent=2)+'\n')
    completed = []
    for i, temperature in enumerate(a.temperature):
        directory = a.work/f'temperature-{i:03d}'
        abundance, _ = opacity_inputs(directory, prepared, spec, a.hydrogen,
                                       a.helium3, temperature)
        execute(prepared['synspec'], directory, ['fort.63', 'fort.29'])
        validate_table(directory/'fort.63', abundance, [temperature], sequence(spec['log_density']))
        completed.append({'temperature_K': temperature, 'table': str((directory/'fort.63').resolve()),
                          'table_sha256': digest(directory/'fort.63'),
                          'receipt_sha256': digest(directory/'completed.json')})
        print(json.dumps(completed[-1]), flush=True)
    if any(digest(path) != checksum for path, checksum in inputs.items()):
        raise ValueError('diagnostic input changed')
    (a.work/'completed.json').write_text(json.dumps({'scope': __doc__,
        'accepted_for_stellar_opacity': False, 'recipe': recipe,
        'prepared_sha256': digest(a.work/'prepared.json'), 'models': completed}, indent=2)+'\n')


if __name__ == '__main__':
    main()
