"""Persistent, cached native EOS adapter for conservative material transport."""
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np

from conservative_material_transport import ThermodynamicPoint, positive_solve
from diffusion_burgers import KB, MU


RGAS = KB/MU


class NativeMaterialThermodynamics:
    def __init__(self, root, work, density, energy_scale):
        self.root, self.work = Path(root), Path(work)
        self.density = np.asarray(density, dtype=float)
        self.energy_scale = energy_scale
        self.cache = {}
        self.input_sha256 = {}
        retained = json.loads((self.root/'docs/results/smooth_eos_refined_retained_sources_v1.json').read_text())
        self.probe = Path('/private/tmp/ember-smooth-eos-four-plane-build-run-v1/probe')
        self.family = Path('/private/tmp/ember-smooth-eos-refined-family-v1/freeeos300_gs98_z020.dat')
        for p in (self.probe, self.family):
            sha = hashlib.sha256(p.read_bytes()).hexdigest()
            assert sha == retained['inputs_sha256'][str(p)]
            self.input_sha256[str(p)] = sha
        old_input = Path('/tmp/ember-smooth-eos-refined-physical-run-v1/retained/retained.input')
        old_output = old_input.with_suffix('.output')
        for p in (old_input, old_output):
            sha = hashlib.sha256(p.read_bytes()).hexdigest()
            assert sha == retained['artifacts_sha256'][str(p)]
            self.input_sha256[str(p)] = sha
        for line, reply in zip(old_input.read_text().splitlines(), old_output.read_text().splitlines(), strict=True):
            x, y, t, rho, chemical = map(float, line.split())
            result = json.loads(reply)
            if chemical == 1 and result['ok']:
                self.cache[x, y, float(np.log(t)), rho] = (t, np.array(result['values']))
        self.process = None
        self.new_queries = 0
        self.reused_queries = 0
        self.maximum_maxwell_error = 0.
        self.replies = (self.work/'native_queries.jsonl').open('x')
        self.stderr = (self.work/'native.stderr').open('x')

    def close(self):
        if self.process is not None:
            self.process.stdin.close()
            code = self.process.wait(timeout=30)
            if code:
                raise RuntimeError(f'native EOS process failed: {code}')
        self.replies.close(); self.stderr.close()

    def __call__(self, i, q, derivatives=True):
        x, y, lt = map(float, q)
        if x < 0 or y < 0 or x+y >= .98:
            raise ValueError('composition outside native material simplex')
        key = (x, y, lt, float(self.density[i]))
        if key in self.cache:
            t, v = self.cache[key]
            self.reused_queries += 1
        else:
            if self.process is None:
                self.process = subprocess.Popen([str(self.probe), str(self.family), 'smooth'],
                    stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=self.stderr,
                    text=True, bufsize=1)
            t = float(np.exp(lt))
            if not np.isfinite(t) or t <= 0:
                raise ValueError('invalid trial temperature')
            row = f'{x:.17g} {y:.17g} {t:.17g} {self.density[i]:.17g} {int(derivatives)}\n'
            # Always request chemical derivatives so a cached reply is complete.
            row = row.rsplit(' ', 1)[0]+' 1\n'
            self.process.stdin.write(row); self.process.stdin.flush()
            reply = self.process.stdout.readline()
            if not reply:
                raise RuntimeError('native EOS process returned no reply')
            result = json.loads(reply)
            self.new_queries += 1
            self.replies.write(json.dumps(dict(input=row.strip(), result=result))+'\n')
            self.replies.flush()
            if not result['ok']:
                raise ValueError(result['error'])
            v = np.asarray(result['values'])
            self.cache[key] = (t, v)
        u = np.array([x, y, v[1]/self.energy_scale])
        if not derivatives:
            return ThermodynamicPoint(u, None, None, None, v[2]/RGAS)
        cv = v[3];ec = v[19:21]
        error = float(np.max(abs(ec+t*v[28:30]))/(RGAS*t+np.max(abs(ec))))
        self.maximum_maxwell_error = max(error, self.maximum_maxwell_error)
        if error > 1e-9 or cv <= 0:
            raise ValueError('EOS chemical/energy derivative identity failed')
        h = np.empty((3, 3))
        h[:2, :2] = v[24:28].reshape(2, 2)+np.outer(ec, ec)/(cv*t*t)
        h[:2, 2] = -self.energy_scale*ec/(cv*t*t)
        h[2, :2] = h[:2, 2]
        h[2, 2] = self.energy_scale**2/(cv*t*t)
        h /= RGAS
        capacity = positive_solve(h, np.eye(3))
        conversion = np.eye(3)
        conversion[2, :2] = -ec/(cv*t)
        conversion[2, 2] = self.energy_scale/(cv*t)
        potential = np.r_[v[22:24], -self.energy_scale/t]/RGAS
        return ThermodynamicPoint(u, potential, capacity, conversion, v[2]/RGAS)
