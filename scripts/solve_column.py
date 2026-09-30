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
After a flux failure, one retry may remove isolated temperature spikes from
the optically thin starting guess. The mesh and acceptance checks stay fixed.
An otherwise converged, flux-rejected column may also restart up to twice from
its complete final profile with an undamped Newton step. The same source checks
and total CPU limit still apply, including when the initial strategy is damped.
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
from generate_nongrey_grid import composition, temperatures, sequence, temperature_convergence_text, resample_initial_structure, continuation_structure
from import_nongrey_grid import source_inputs, source_state

INPUTS = ['fort.5', 'tas', 'ember-masses.dat', 'fort.15', 'fort.8']
RETAINED_SOURCE_FILES = [*INPUTS, 'fort.7', 'fort.9', 'fort.9.gz', 'run.log', 'run.log.gz',
                         'specification.json', 'specification.prepared.json', 'provenance.json',
                         'physics.json', 'opacity.sha256']

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

def upper_layer_guess(text, log, nd):
    """Prepare a new guess from an otherwise converged, flux-rejected column.

    Isolated hot roots in optically thin layers can escape the source's energy
    equation while its independent convective-flux diagnostic rejects them.
    Interpolate those guess points between their neighbours, then solve again.
    No modified state is accepted without the usual full source checks.
    """
    if not valid_structure(text, nd) or 'FINAL MODEL ATMOSPHERE' not in log:
        return None
    profile = []
    for line in log.rsplit('FINAL MODEL ATMOSPHERE', 1)[1].splitlines():
        words = line.replace('D', 'E').split()
        if len(words) == 11 and words[0].isdigit():
            try: profile.append(list(map(float, words)))
            except ValueError: return None
    if len(profile) != nd or any(r[0] != i+1 or not all(math.isfinite(v) for v in r)
                                 for i, r in enumerate(profile)):
        return None
    words = text.replace('D', 'E').split()
    values = list(map(float, words[2:])); mass = values[:nd]
    state = [values[nd+4*i:nd+4*i+4] for i in range(nd)]
    if any(not math.isclose(r[1], m, rel_tol=1e-6) or
           not math.isclose(r[3], v[0], rel_tol=1e-6)
           for r, m, v in zip(profile, mass, state)):
        return None
    bad = [i for i, r in enumerate(profile) if abs(r[10]-1) > .002]
    if not bad or len(bad) > 3:
        return None
    if any(i == 0 or i == nd-1 or i-1 in bad or i+1 in bad or
           not 0 < profile[i][2] < 1e-5 or abs(profile[i][8]-1) > .002 or
           profile[i][9] <= 0 or state[i][0] <= 1.005*max(state[i-1][0], state[i+1][0])
           for i in bad):
        return None
    for i in bad:
        w = math.log(mass[i]/mass[i-1])/math.log(mass[i+1]/mass[i-1])
        state[i] = [math.exp((1-w)*math.log(a)+w*math.log(b))
                    for a, b in zip(state[i-1], state[i+1])]
    guess = (f'{nd} -4\n' + '\n'.join(format(v, '.17g') for v in mass) + '\n'
             + '\n'.join(' '.join(format(v, '.17g') for v in row) for row in state) + '\n')
    return (guess, [i+1 for i in bad]) if valid_structure(guess, nd) else None

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

def recover_damped_budget(result, attempt, nd):
    """Reuse an improving damped solve after its per-attempt CPU limit.

    The total CPU and phase limits still apply. Only a complete saved guess
    is reused; the interrupted calculation is never accepted as a source.
    """
    h = result.get('history', [])
    saved = attempt / 'best.fort.7'
    return (result['returncode'] == -signal.SIGXCPU and len(h) >= 4
            and 0 < min(h[-3:]) < min(1e-2, .25*h[0])
            and saved.is_file() and valid_structure(saved.read_text(), nd))

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

def donor_guess(donor, col, spec, nd):
    """Initial guess from an accepted neighbouring column (a guess only: the solve and acceptance are unchanged).

    T scales with the Teff ratio and column mass with 1/g (continuation_structure keeps n k T and g m fixed), then the
    profile is resampled to this column's depth count. The donor must be accepted, unchanged since its result, and
    have the same source physics (assemble_nongrey_grid.physical_identity); its composition may differ.
    """
    from assemble_nongrey_grid import physical_identity
    record = json.loads((donor / 'result.json').read_text())
    if not record.get('accepted'): raise SystemExit('donor column is not accepted')
    for name, checksum in record.get('retained_sha256', {}).items():
        if name in ('fort.7', 'specification.json') and sha(donor / name) != checksum:
            raise SystemExit(f'donor output changed since its result: {name}')
    dspec = json.loads((donor / 'specification.json').read_text())
    if physical_identity(dspec, json.loads((donor / 'provenance.json').read_text())) != \
            physical_identity(spec, json.loads((col / 'provenance.json').read_text())):
        raise SystemExit('donor source physics differs from this column')
    dt, dg = dspec['teff_K'][0], dspec['log_g'][0]
    text = continuation_structure((donor / 'fort.7').read_text(), dt, dg, spec['teff_K'][0], spec['log_g'][0])
    text = resample_initial_structure(text, nd, allow_coarsen=True)
    return text, dict(directory=str(donor), fort7_sha256=sha(donor / 'fort.7'), result_sha256=sha(donor / 'result.json'),
                      teff_K=dt, log_g=dg, hydrogen=dspec['hydrogen'][0], helium3=dspec['helium3'][0])

def main():
    solver_sha256 = sha(Path(__file__))
    a = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    a.add_argument('column', type=Path); a.add_argument('--executable', required=True)
    a.add_argument('--executable-sha256', required=True)
    a.add_argument('--attempt-cpu', type=float, default=600); a.add_argument('--total-cpu', type=float, default=1500)
    a.add_argument('--strategy', choices=['hybrid', 'damped'], default='hybrid')
    a.add_argument('--chmax', type=float, default=None, help='optional source stop for every attempt (e.g. 1e-9)')
    a.add_argument('--max-phases', type=int, default=5)
    a.add_argument('--initial-from', type=Path, default=None,
                   help='accepted donor column: scale its final structure to this Teff/log g as the initial guess')
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
    # Check mixture metadata before spending CPU on the source solve.
    composition(x, y, spec0['metals'])
    tas0 = (col / 'tas').read_text(); nd = int(re.search(r'\bND\s*=\s*(\d+)', tas0).group(1))
    m = re.search(r'\bORELAX\s*=\s*([0-9.eE+-]+)', tas0); orelax0 = float(m.group(1)) if m else 1.0
    if args.strategy == 'hybrid' and orelax0 == 1.0: orelax0 = 0.3
    opacity = {'temperature_K': temperatures(spec0), 'density_g_cm3': sequence(spec0['log_density'])}
    guess = (col / 'fort.8').read_text() if (col / 'fort.8').exists() else None
    donor = None
    if args.initial_from is not None:
        guess, donor = donor_guess(args.initial_from.resolve(), col, spec0, nd)
    attempts, total, kind = [], 0.0, ('newton' if args.strategy == 'hybrid' else 'damped')
    budget_exhausted = False
    (col / 'attempts').mkdir()
    accepted = None; upper_guess_retried = False; flux_restarts = 0
    phase_limit = args.max_phases if args.strategy == 'hybrid' else 3
    for phase in range(phase_limit):
        budget = min(args.attempt_cpu, args.total_cpu - total)
        if budget < 30:
            budget_exhausted = True; break
        if args.strategy == 'hybrid' and phase == args.max_phases - 1 and not flux_restarts:
            kind = 'damped'  # Final fallback, except for a requested undamped flux restart.
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
            convergence = att / 'fort.9'
            changes = convergence.read_text(errors='ignore') if convergence.exists() else ''
            r['diagnostics'] = source_state(logtext, changes, t, g, opacity, spec['tau'])
            r['accepted'] = True
        except Exception as e:
            r['failure'] = str(e)[:300]
        attempts.append(r)
        if r['accepted']: accepted = (att, spec, r); break
        if (not upper_guess_retried and args.strategy == 'hybrid' and phase+1 < args.max_phases
                and r.get('failure', '').startswith('unconverged source:')
                and r['history'] and r['history'][-1] <= 1e-6 and (att / 'fort.7').exists()):
            recovered = upper_layer_guess((att / 'fort.7').read_text(), logtext, nd)
            if recovered:
                guess, depths = recovered; upper_guess_retried = True
                r['next_guess'] = dict(method='interpolate isolated upper-layer spikes', depths=depths,
                                      parent_sha256=sha(att / 'fort.7'))
                kind = 'newton'
                continue
        if (r.get('failure', '').startswith('unconverged source:')
                and r['returncode'] == 0 and not r['monitor_stop']
                and r['history'] and r['history'][-1] <= 1e-6):
            saved = att / 'fort.7'
            if (flux_restarts < 2 and phase+1 < phase_limit and saved.exists()
                    and valid_structure(saved.read_text(), nd)):
                guess = saved.read_text(); flux_restarts += 1; kind = 'newton'
                r['next_guess'] = dict(method='undamped restart of complete flux-rejected profile; unchanged physics',
                                      parent_sha256=sha(saved))
                continue
            break
        if args.strategy == 'damped':
            break  # Extra phases are reserved for flux restarts.
        if (att / 'best.fort.7').exists(): guess = (att / 'best.fort.7').read_text()
        stop = r.get('monitor_stop') or ''
        if kind == 'newton': kind = 'damped'
        elif stop.startswith('switch'): kind = 'newton'
        elif recover_damped_budget(r, att, nd):
            r['next_guess'] = dict(method='Newton restart after improving damped CPU-limited attempt',
                                  parent_sha256=sha(att / 'best.fort.7'))
            kind = 'newton'
        else: break
    for att in (col / 'attempts').iterdir():        # compact every attempt's logs
        for f in ['fort.9', 'run.log']:
            p = att / f
            if p.exists() and p.stat().st_size > 200000:
                (att / (f + '.gz')).write_bytes(gzip.compress(p.read_bytes(), compresslevel=6, mtime=0)); p.unlink()
        for p in att.glob('fort.*'):
            if p.name not in ('fort.5', 'fort.7', 'fort.8', 'fort.9', 'fort.9.gz', 'fort.15'): p.unlink()
    result = dict(name=col.name, strategy=args.strategy, accepted=accepted is not None, initial_from=donor,
                  CPU_seconds=total, wall_seconds=sum(r['wall_seconds'] for r in attempts),
                  solver_sha256=solver_sha256,
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
    # A caller may redirect our stdout into this directory and keep writing after
    # result.json is closed. Hash the scientific record, not caller-owned logs.
    result['retained_sha256'] = {name: sha(col / name) for name in RETAINED_SOURCE_FILES
                                 if (col / name).is_file()}
    (col / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['name', 'accepted', 'CPU_seconds']} | {'attempts': len(attempts)}))

if __name__ == '__main__':
    main()
