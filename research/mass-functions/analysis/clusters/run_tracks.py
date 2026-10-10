"""Fresh independent Ember tracks; two simultaneous two-thread jobs at most."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess, time, json, hashlib, resource, datetime

ROOT = Path(__file__).resolve().parent
BIN = Path('@EMBER_ROOT@/out/luhman16-independent-oct8-v1/build/apps/ember-evolve')

def run(case):
    d = ROOT / 'tracks' / case
    cfg = d / 'configuration.txt'
    scope=json.loads((d/'scope.json').read_text()) if (d/'scope.json').exists() else {}
    binary = Path(scope.get('executable', BIN))
    cmd = [str(binary), '--lifetime', str(cfg), str(scope.get('target_years',800000000)), str(d/'evolution'),
           '--max-steps', str(scope.get('maximum_steps',1000)),
           '--cpu-seconds', str(scope.get('maximum_cpu_seconds',600))]
    receipt = dict(command=cmd, started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   source_commit=scope.get('source_description','087351c'),
                   binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                   configuration_sha256=hashlib.sha256(cfg.read_bytes()).hexdigest())
    start = time.monotonic()
    with (d/'stdout.log').open('w') as stdout, (d/'stderr.log').open('w') as stderr:
        p = subprocess.Popen(cmd, stdout=stdout, stderr=stderr)
        receipt['pid'] = p.pid
        (d/'launch.json').write_text(json.dumps(receipt, indent=2)+'\n')
        rc = p.wait()
    receipt.update(returncode=rc, wall_seconds=time.monotonic()-start,
                   ended_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    report = d/'evolution/report.json'
    if report.exists(): receipt['report'] = json.loads(report.read_text())
    (d/'completion.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(case, rc, round(receipt['wall_seconds'],1), flush=True)
    return receipt

if __name__ == '__main__':
    import sys
    cases = sys.argv[1:] or ['control_m0.25_z0.02']
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(run,cases))
