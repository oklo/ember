#!/usr/bin/env python3
"""Check depletion with fixed bulk input and a synthetic mutable-bulk caller.

MOLEQ keeps cbase separate from its mutable depleted gas arrays. The fixed-bulk
sequence models that interface. Feeding gas output back as bulk input is a
separate robustness experiment, not evidence of a defect in actual atmospheres.
This tests adapter state management, not full Fortran integration or grain EOS.
"""
import argparse
import json
from pathlib import Path
import subprocess

from generate_nongrey_grid import composition
from prepare_nongrey_sources import digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['model', 'adapter', 'library', 'work', 'output']:
        p.add_argument(name, type=Path)
    a = p.parse_args()
    if a.work.exists():
        raise FileExistsError('use a fresh adapter audit directory')
    saved = json.loads((a.model/'provenance.json').read_text())
    depletion = saved['prepared']['depletion']; fc = depletion['fastchem']
    if digest(a.library) != depletion['library_sha256']:
        raise ValueError('FastChem library differs from the pinned source model')
    source = Path(fc['source'])
    probe = Path(__file__).with_name('condensate_adapter_probe.cpp')
    dependencies = {str(path.resolve()): digest(path) for path in
                    [a.model/'provenance.json', a.model/'condensate-abundances.dat',
                     a.adapter, a.library, probe, Path(__file__)]}
    for name, expected in fc['chemistry_data_sha256'].items():
        path = source/'input/logK'/name
        if digest(path) != expected:
            raise ValueError('FastChem chemistry changed')
        dependencies[str(path.resolve())] = expected
    a.work.mkdir(parents=True); a.work = a.work.resolve()
    (a.work/'adapter.cpp').write_bytes(a.adapter.read_bytes())
    (a.work/'probe.cpp').write_bytes(probe.read_bytes())
    command = ['c++', '-std=c++17', '-O3', '-I'+str(source/'fastchem_src'),
               str(a.work/'probe.cpp'), str(a.work/'adapter.cpp'), str(a.library.resolve()),
               '-o', str(a.work/'probe')]
    with (a.work/'build.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    (a.work/'ember-condensates.cfg').write_text('equilibrium\n'+str(source)+'\n'+
                                               str((a.model/'condensate-abundances.dat').resolve())+'\n')
    bulk, _ = composition(saved['XH'], saved['X3'], saved['specification']['metals'])
    bulk = bulk[:92]
    prefix = ' '.join(format(v, '.17g') for v in bulk)+'\n'
    states = [(7000., 4000.), (1400., 4000.), (1769.41749, 4061.8898), (2400., 4000.)]

    def calculate(name, sequence, mutable=False):
        text = prefix+''.join(f'{t:.17g} {pressure:.17g} {int(i == 0)}\n'
                              for i, (t, pressure) in enumerate(sequence))
        (a.work/(name+'.input')).write_text(text)
        command = [str(a.work/'probe')]+(['--mutable-bulk'] if mutable else [])
        run = subprocess.run(command, input=text, text=True,
                             capture_output=True, cwd=a.work, check=True)
        (a.work/(name+'.txt')).write_text(run.stdout)
        (a.work/(name+'.log')).write_text(run.stderr)
        rows = [list(map(float, line.split())) for line in run.stdout.splitlines()]
        if len(rows) != len(sequence) or any(len(r) != 92 for r in rows):
            raise ValueError('invalid adapter probe output')
        return rows

    reference = {state: calculate('reference-'+str(i), [state])[0]
                 for i, state in enumerate(states)}
    order = [0, 1, 0, 3, 1, 3, 2, 0, 1, 2, 3, 0]
    sequence = [states[i] for i in order]
    comparisons = {}
    for mode, mutable in [('fixed_bulk', False), ('synthetic_mutable_bulk', True)]:
        rows = calculate(mode+'-sequence', sequence, mutable)
        defects = [max(abs(x-y)/b for x, y, b in zip(row, reference[state], bulk, strict=True)
                       if b > 1e-90) for row, state in zip(rows, sequence, strict=True)]
        comparisons[mode] = {'maximum_element_change_in_bulk_abundance_units': max(defects),
                             'sequence_element_defects': defects, 'passed': max(defects) < 1e-6}
    cold_depletion = max(abs(x/b-1) for x, b in zip(reference[states[1]], bulk, strict=True) if b > 1e-90)
    if cold_depletion < .01:
        raise ValueError('regression has no significant condensation')
    passed = comparisons['fixed_bulk']['passed']
    if any(digest(path) != expected for path, expected in dependencies.items()):
        raise ValueError('adapter audit inputs changed')
    report = {'scope': __doc__, 'input_files_sha256': dependencies, 'build_command': command,
              'executable_sha256': digest(a.work/'probe'), 'sequence': sequence,
              'comparisons': comparisons, 'maximum_cold_depletion': cold_depletion,
              'tolerance': 1e-6, 'passed': passed,
              'outputs_sha256': {f.name: digest(f) for f in a.work.iterdir() if f.is_file()}}
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'comparisons': comparisons, 'passed': passed}))
    if not passed:
        raise SystemExit('depletion depends on previous state with fixed bulk input')


if __name__ == '__main__':
    main()
