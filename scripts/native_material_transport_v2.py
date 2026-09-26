"""Native material adapter with thermal-only old states at zero abundance."""
import json
import subprocess
import numpy as np

from conservative_material_transport import ThermodynamicPoint
from native_material_transport import NativeMaterialThermodynamics as InteriorThermodynamics, RGAS


class NativeMaterialThermodynamics(InteriorThermodynamics):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.thermal_cache = {}
        self.thermal_only_queries = 0

    def __call__(self, i, q, derivatives=True):
        x, y, lt = map(float, q)
        key = (x, y, lt, float(self.density[i]))
        if derivatives or key in self.cache:
            return super().__call__(i, q, derivatives=derivatives)
        if x < 0 or y < 0 or x+y >= .98:
            raise ValueError('composition outside native material simplex')
        if key in self.thermal_cache:
            self.reused_queries += 1
            v = self.thermal_cache[key]
        else:
            t = float(np.exp(lt))
            if not np.isfinite(t) or t <= 0:
                raise ValueError('invalid trial temperature')
            if self.process is None:
                self.process = subprocess.Popen([str(self.probe), str(self.family), 'smooth'],
                    stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=self.stderr,
                    text=True, bufsize=1)
            row = f'{x:.17g} {y:.17g} {t:.17g} {self.density[i]:.17g} 0\n'
            self.process.stdin.write(row);self.process.stdin.flush()
            reply = self.process.stdout.readline()
            if not reply:
                raise RuntimeError('native EOS process returned no thermal reply')
            result = json.loads(reply)
            self.new_queries += 1;self.thermal_only_queries += 1
            self.replies.write(json.dumps(dict(input=row.strip(), result=result))+'\n')
            self.replies.flush()
            if not result['ok']:
                raise ValueError(result['error'])
            v = np.asarray(result['values'])
            self.thermal_cache[key] = v
        return ThermodynamicPoint(np.array([x, y, v[1]/self.energy_scale]),
                                   None, None, None, v[2]/RGAS)
