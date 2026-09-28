#!/usr/bin/env python3
"""Solve one prepared TLUSTY source column with an adaptive Newton/damped schedule.

usage: solve_column.py COLUMN_DIR --executable EXE --executable-sha256 SHA
                       [--attempt-cpu 600] [--total-cpu 1500] [--strategy hybrid|damped] [--chmax 1e-9]

COLUMN_DIR contains prepared TLUSTY inputs (fort.5, tas, ember-masses.dat, fort.15,
opacity.bin, data, optional fort.8), specification.json, provenance.json and opacity.sha256.
Only the numerical controls ORELAX (Newton relaxation) and, optionally, CHMAX
(temperature convergence) vary between attempts. Every attempt runs in COLUMN_DIR/attempts/NN on a copy of the
inputs. The final accepted attempt is promoted into COLUMN_DIR with its own tas, initial structure (fort.8) and outputs.
specification.json records the numerical controls it used; the original is kept as specification.prepared.json.
Its result.json has the same fields as the existing runners (accepted, CPU_seconds, diagnostics, retained_sha256),
so extend_mixed_atmosphere.load_column revalidates it with the unchanged source_inputs/source_state.

Schedule ("hybrid"):
  1. undamped probe (ORELAX=1) from the prepared guess; stop on growth or stall (at most 30 iterations);
  2. otherwise damped (the prepared ORELAX, default 0.3) from the best fully written iterate so far;
     hand over to undamped once max correction < 1e-2 and falling for 3 iterations;
  3. alternate at most 5 phases; the last phase always runs damped to TLUSTY's own end.
A monitor stop never claims that no solution exists. A column that fails every phase is reported unaccepted, with
per-attempt histories, so a different initial guess (continuation) can be tried.
"""
import argparse, gzip, hashlib, json, math, os, re, resource, shutil, signal, subprocess, sys, time
from pathlib import Path
def _scripts_dir():
    """Locate the Ember scripts directory: EMBER_ROOT, this file's own directory, or an ancestor containing it."""
    if os.environ.get('EMBER_ROOT'):
        return Path(os.environ['EMBER_ROOT']).resolve() / 'scripts'
    here = Path(__file__).resolve().parent
    for p in [here, *here.parents]:
        for c in (p, p / 'scripts'):
            if (c / 'generate_nongrey_grid.py').exists() and (c / 'import_nongrey_grid.py').exists():
                return c
    raise SystemExit('cannot locate Ember scripts/ (set EMBER_ROOT)')
sys.path.insert(0, str(_scripts_dir()))
from generate_nongrey_grid import temperatures, sequence, temperature_convergence_text, resample_initial_structure
from import_nongrey_grid import source_inputs, source_state

INPUTS = ['fort.5', 'tas', 'ember-masses.dat', 'fort.15', 'fort.8']

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def complete_history(text, nd):
    it = {}
    for line in text.splitlines():
        w = line.replace('D', 'E').split()
        if len(w) == 9 and w[0].isdigit() and w[1].isdigit():
            try: it.setdefault(int(w[0]), []).append(abs(float(w[6])))
            except ValueError: pass
    return [max(v) for k, v in sorted(it.items()) if len(v) == nd]

def valid_structure(text, nd):
    """A snapshot is reused only if it is a complete, finite, positive TLUSTY structure of the right depth count."""
    try:
        w = text.replace('D', 'E').split()
        if int(w[0]) != nd or int(w[1]) != -4 or len(w) != 2 + 5 * nd: return False
        v = [float(x) for x in w[2:]]
        if not all(math.isfinite(x) and x > 0 for x in v): return False
        resample_initial_structure(text, nd)
        return True
    except Exception:
        return False

def set_param(tas, key, value):
    pat = re.compile(r'(\b' + key + r'\s*=\s*)([0-9.eEdD+-]+)')
    if pat.search(tas): return pat.sub(lambda m: m.group(1) + value, tas, count=1)
    return tas.rstrip('\n') + f'\n{key}={value}\n'

def newton_monitor(h):
    if len(h) < 3: return None
    if h[-1] > 30 * min(h) and h[-1] > 1e-3: return 'undamped growth'
    if len(h) >= 8 and min(h[-5:]) > min(h[:-5]) / 2: return 'undamped stall'
    if len(h) >= 30: return 'undamped iteration budget'
    return None

def damped_switch(h, switch=1e-2):
    if len(h) >= 4 and h[-1] < switch and h[-1] < h[-2] < h[-3] < h[-4]: return 'switch to undamped'
    return None

def stop_and_reap(pid):
    """Kill a source process group and reap it; the process may already have exited after the last poll."""
    try: os.killpg(pid, signal.SIGKILL)
    except ProcessLookupError: pass
    return os.wait4(pid, 0)

def run(att, exe, cap, nd, rule):
    env = dict(os.environ, OMP_NUM_THREADS='1')
    def limit(): resource.setrlimit(resource.RLIMIT_CPU, (int(cap), int(cap) + 10))
    log = (att / 'run.log').open('wb'); stdin = (att / 'fort.5').open('rb'); t0 = time.monotonic()
    p = subprocess.Popen([exe], cwd=att, stdin=stdin, stdout=log, stderr=subprocess.STDOUT, env=env,
                         start_new_session=True, preexec_fn=limit)
    stdin.close(); stop = None; seen = 0; best = math.inf
    while True:
        pid, status, usage = os.wait4(p.pid, os.WNOHANG)
        if pid: break
        if (att / 'fort.9').exists():
            h = complete_history((att / 'fort.9').read_text(errors='ignore'), nd)
            if len(h) > seen:
                seen = len(h)
                if h[-1] < best and (att / 'fort.7').exists():
                    text = (att / 'fort.7').read_text(errors='ignore')
                    if valid_structure(text, nd):          # fully written snapshot only
                        best = h[-1]; (att / 'best.fort.7').write_text(text)
                if rule: stop = rule(h)
                if stop:
                    pid, status, usage = stop_and_reap(p.pid); break
        time.sleep(0.5)
    p.returncode = os.waitstatus_to_exitcode(status)
    log.close()
    h = complete_history((att / 'fort.9').read_text(errors='ignore'), nd) if (att / 'fort.9').exists() else []
    return dict(returncode=os.waitstatus_to_exitcode(status), CPU_seconds=usage.ru_utime + usage.ru_stime,
                wall_seconds=time.monotonic() - t0, monitor_stop=stop, iterations=len(h),
                history=[float('%.3g' % x) for x in h])

def main():
    a = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    a.add_argument('column', type=Path); a.add_argument('--executable', required=True)
    a.add_argument('--executable-sha256', required=True)
    a.add_argument('--attempt-cpu', type=float, default=600); a.add_argument('--total-cpu', type=float, default=1500)
    a.add_argument('--strategy', choices=['hybrid', 'damped'], default='hybrid')
    a.add_argument('--chmax', type=float, default=None, help='optional source stop for every attempt (e.g. 1e-9)')
    a.add_argument('--max-phases', type=int, default=5)
    args = a.parse_args()
    col = args.column.resolve()
    args.executable = str(Path(args.executable).resolve())          # attempts run with cwd inside the column
    if sha(args.executable) != args.executable_sha256: raise SystemExit('executable checksum mismatch')
    if (col / 'result.json').exists(): raise SystemExit('column already has a result; use a fresh prepared column')
    if (col / 'attempts').exists():
        raise SystemExit('attempts/ already exists (an interrupted run?); inspect it and use a fresh prepared column')
    for name in ('attempt_cpu', 'total_cpu'):
        v = getattr(args, name)
        if not (math.isfinite(v) and v > 0): raise SystemExit(f'--{name.replace("_", "-")} must be finite and positive')
    if args.max_phases < 1: raise SystemExit('--max-phases must be >= 1')
    if args.chmax is not None and not (math.isfinite(args.chmax) and 0 < args.chmax <= 1e-6):
        raise SystemExit('--chmax must lie in (0, 1e-6]')
    prepared = json.loads((col / 'provenance.json').read_text())
    if prepared['executables']['tlusty'] != args.executable_sha256:
        raise SystemExit('executable differs from the prepared source physics')
    if sha(col / 'opacity.bin') != (col / 'opacity.sha256').read_text().strip():
        raise SystemExit('opacity differs from the prepared source table')
    spec0 = json.loads((col / 'specification.json').read_text())
    x, y, t, g = (spec0[k][0] for k in ['hydrogen', 'helium3', 'teff_K', 'log_g'])
    tas0 = (col / 'tas').read_text(); nd = int(re.search(r'\bND\s*=\s*(\d+)', tas0).group(1))
    m = re.search(r'\bORELAX\s*=\s*([0-9.eE+-]+)', tas0); orelax0 = float(m.group(1)) if m else 1.0
    if args.strategy == 'hybrid' and orelax0 == 1.0: orelax0 = 0.3
    opacity = {'temperature_K': temperatures(spec0), 'density_g_cm3': sequence(spec0['log_density'])}
    guess = (col / 'fort.8').read_text() if (col / 'fort.8').exists() else None
    attempts, total, kind = [], 0.0, ('newton' if args.strategy == 'hybrid' else 'damped')
    budget_exhausted = False
    (col / 'attempts').mkdir()
    accepted = None
    for phase in range(args.max_phases if args.strategy == 'hybrid' else 1):
        budget = min(args.attempt_cpu, args.total_cpu - total)
        if budget < 30:
            budget_exhausted = True; break
        if args.strategy == 'hybrid' and phase == args.max_phases - 1:
            kind = 'damped'                               # the final phase is always the damped fallback
        att = col / 'attempts' / f'{phase:02d}_{kind}'
        att.mkdir(parents=True, exist_ok=False)
        for f in INPUTS:
            if (col / f).exists(): shutil.copyfile(col / f, att / f)
        for f in ['opacity.bin', 'data']: os.symlink(os.path.realpath(col / f), att / f)
        spec = dict(spec0); orelax = 1.0 if kind == 'newton' else orelax0
        tas = set_param(tas0, 'ORELAX', '1' if orelax == 1 else repr(orelax)); spec['newton_relaxation'] = orelax
        if args.chmax is not None:
            tas = set_param(tas, 'CHMAX', temperature_convergence_text(args.chmax)); spec['temperature_convergence'] = args.chmax
        (att / 'tas').write_text(tas)
        if guess is not None:
            (att / 'fort.8').write_text(guess)
            f5 = (att / 'fort.5').read_text()
            if '\nT T\n' in f5: (att / 'fort.5').write_text(f5.replace('\nT T\n', '\nT F\n', 1))
        last = args.strategy == 'damped' or phase == args.max_phases - 1
        rule = None if last else (newton_monitor if kind == 'newton' else damped_switch)
        r = run(att, args.executable, budget, nd, rule); total += r['CPU_seconds']
        r.update(phase=phase, kind=kind, ORELAX=orelax, CHMAX=args.chmax, accepted=False)
        try:
            if r['returncode']: raise RuntimeError(f"source process exit {r['returncode']}")
            if r['monitor_stop']: raise RuntimeError('monitor stop: ' + r['monitor_stop'])
            logtext = (att / 'run.log').read_text(errors='ignore')
            inputs = {k: (att / f).read_text() for k, f in [('atmosphere_input', 'fort.5'), ('parameters', 'tas'),
                      ('element_masses', 'ember-masses.dat')]}
            if (att / 'fort.8').exists(): inputs['initial_structure'] = (att / 'fort.8').read_text()
            source_inputs(inputs, spec, x, y, t, g, logtext)
            r['diagnostics'] = source_state(logtext, (att / 'fort.9').read_text(errors='ignore'), t, g, opacity, spec['tau'])
            r['accepted'] = True
        except Exception as e:
            r['failure'] = str(e)[:300]
        attempts.append(r)
        if r['accepted']: accepted = (att, spec, r); break
        if (att / 'best.fort.7').exists(): guess = (att / 'best.fort.7').read_text()
        stop = r.get('monitor_stop') or ''
        if kind == 'newton': kind = 'damped'
        elif stop.startswith('switch'): kind = 'newton'
        else: break
    for att in (col / 'attempts').iterdir():        # compact every attempt's logs
        for f in ['fort.9', 'run.log']:
            p = att / f
            if p.exists() and p.stat().st_size > 200000:
                (att / (f + '.gz')).write_bytes(gzip.compress(p.read_bytes(), compresslevel=6, mtime=0)); p.unlink()
        for p in att.glob('fort.*'):
            if p.name not in ('fort.5', 'fort.7', 'fort.8', 'fort.9', 'fort.9.gz', 'fort.15'): p.unlink()
    result = dict(name=col.name, strategy=args.strategy, accepted=accepted is not None,
                  CPU_seconds=total, wall_seconds=sum(r['wall_seconds'] for r in attempts),
                  selected_for_evolution=False, executable_sha256=args.executable_sha256,
                  attempts=[{k: v for k, v in r.items() if k != 'diagnostics'} for r in attempts])
    if accepted:
        att, spec, r = accepted
        if not (col / 'specification.prepared.json').exists():
            shutil.copyfile(col / 'specification.json', col / 'specification.prepared.json')
        (col / 'specification.json').write_text(json.dumps(spec, indent=2) + '\n')
        for f in ['tas', 'fort.5', 'fort.8', 'fort.7', 'fort.9', 'fort.9.gz', 'run.log', 'run.log.gz']:
            if (att / f).exists(): shutil.copyfile(att / f, col / f)
        result['diagnostics'] = r['diagnostics']
        result['final_attempt'] = str(att.relative_to(col))
    else:
        result['failure'] = ('CPU budget exhausted before another attempt (< 30 s left)' if budget_exhausted and attempts
                             else attempts[-1].get('failure') if attempts else 'no attempt within budget (< 30 CPU s)')
    result['retained_sha256'] = {p.name: sha(p) for p in sorted(col.iterdir()) if p.is_file() and not p.is_symlink()
                                 and p.name != 'result.json'}
    (col / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['name', 'accepted', 'CPU_seconds']} | {'attempts': len(attempts)}))

if __name__ == '__main__':
    main()
