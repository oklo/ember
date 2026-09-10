#!/usr/bin/env python3
"""Run independent, resumable atmosphere continuations from an explicit plan.

Only completed canonical runs count as accepted models. This controller
does not assemble or install a grid, or start a stellar evolution run.
An optional cancellation.json in the work directory, containing its absolute
work path, saved plan_sha256 and reason, prevents subsequent source launches.
Already running source calculations finish and retain their outputs.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import subprocess
import sys
import time

from prepare_nongrey_sources import digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('plan', type=Path)
    p.add_argument('work', type=Path)
    p.add_argument('--jobs', type=int, default=2)
    p.add_argument('--dependency-timeout', type=float, default=14400,
                   help='seconds to wait for an explicitly named canonical seed')
    a = p.parse_args()
    if not 1 <= a.jobs <= 8 or not 0 < a.dependency_timeout <= 86400:
        raise ValueError('invalid worker count or dependency timeout')
    plan = json.loads(a.plan.read_text())
    seen = set()
    for job in plan['requests']:
        name = job['name']
        if not name or Path(name).name != name or name in ['.', '..'] or name in seen:
            raise ValueError('unsafe or duplicate work label')
        seen.add(name)
    a.work.mkdir(parents=True, exist_ok=True)
    saved = a.work/'plan.json'
    if saved.exists() and json.loads(saved.read_text()) != plan:
        raise ValueError('changed plan requires a new work directory')
    saved.write_text(json.dumps(plan, indent=2)+'\n')

    def calculate(job):
        seed = Path(job['initial'])
        deadline = time.monotonic()+a.dependency_timeout
        while not (seed/'validated.json').exists():
            if not job.get('wait_for_seed'):
                raise FileNotFoundError(f'canonical seed is not ready: {seed}')
            if time.monotonic() > deadline:
                raise TimeoutError(f'canonical seed did not finish: {seed}')
            time.sleep(10)
        while not Path(job['opacity']).exists():
            if not job.get('wait_for_opacity'):
                raise FileNotFoundError('source opacity table is not ready: '+job['opacity'])
            if time.monotonic() > deadline:
                raise TimeoutError('source opacity table did not finish: '+job['opacity'])
            time.sleep(10)
        command = [sys.executable, '-B', str(Path(__file__).with_name('run_nongrey_continuation.py')),
                   plan['prepared'], job['specification'], job['opacity'], str(seed), str(a.work/job['name']),
                   '--initializer', plan['initializer'], '--hydrogen', str(job['XH']),
                   '--helium3', str(job['X3']), '--teff', str(job['teff_K']), '--logg', str(job['log_g']),
                   '--initializer-depths', str(plan['initializer_depths']),
                   '--initializer-frequencies', str(plan['initializer_frequencies'])]
        if job.get('composition_continuation'):
            command.append('--composition-continuation')
        if 'initial_bottom_tau' in job:
            command += ['--initial-bottom-tau', str(job['initial_bottom_tau'])]
        if 'convective_iterations' in job:
            command += ['--convective-iterations', str(job['convective_iterations'])]
        with (a.work/(job['name']+'.log')).open('a') as log:
            print('command:', json.dumps(command), file=log, flush=True)
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        validation = a.work/job['name']/'final/validated.json'
        record = json.loads(validation.read_text())
        if any(record[k] != job[k] for k in ['XH', 'X3', 'teff_K', 'log_g']):
            raise ValueError('final model label differs from request')
        return {'name': job['name'], 'status': 'canonical_validated',
                'validation_sha256': digest(validation), 'diagnostics': record['diagnostics']}

    results = {}
    with ThreadPoolExecutor(a.jobs) as pool:
        futures = {pool.submit(calculate, job): job for job in plan['requests']}
        for future in as_completed(futures):
            job = futures[future]
            try:
                result = future.result()
            except Exception as error:
                result = {'name': job['name'], 'status': 'rejected', 'error': str(error)}
            results[job['name']] = result
            print(json.dumps(result), flush=True)
            report = {'scope': __doc__, 'plan_sha256': digest(a.plan),
                      'complete': len(results) == len(plan['requests']),
                      'models': [results[j['name']] for j in plan['requests'] if j['name'] in results]}
            temporary = a.work/'manifest.pending'
            temporary.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
            temporary.replace(a.work/'manifest.json')
    if any(r['status'] != 'canonical_validated' for r in results.values()):
        raise SystemExit('some source models were rejected; completed models remain reusable')


if __name__ == '__main__':
    main()
