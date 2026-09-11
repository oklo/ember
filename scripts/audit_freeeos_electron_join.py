#!/usr/bin/env python3
"""Compare FreeEOS free-energy changes with integrated source derivatives.

The density interval crosses the electron-fit switch near degeneracy parameter
4 at two mixtures; a third mixture is a control. Option 1 uses the fitted
electrons and includes radiation. Option 223 integrates the electrons and
omits radiation. Each is tested against its own thermodynamic derivative;
their raw energies must not be compared without accounting for radiation.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess

from metal_eos_composition import mixture


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('probe', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--jobs', type=int, default=4)
    p.add_argument('--source-directory', type=Path)
    p.add_argument('--source-library', type=Path)
    a = p.parse_args()
    if not 1 <= a.jobs <= 8:
        raise ValueError('invalid workers')
    T0 = 11750953.02469058
    rho0 = 4004.2197096372915
    q0 = math.log10(rho0)-1.5*(math.log10(T0)-6)
    cases = [(x,y,direction,option) for x,y in [(0.0375,0.12),(0.0375,0),(0.075,0)]
             for direction,option in [('temperature',1),('density',1),('density',223)]]
    inputs = {str(a.probe.resolve()): digest(a.probe), str(Path(__file__).resolve()): digest(__file__)}
    if a.source_directory:
        for name in ['mod_free_eos.f90','fermi_dirac.f90','fermi_dirac_ct.f90']:
            path = a.source_directory/'src'/name
            inputs[str(path.resolve())] = digest(path)
    if a.source_library:
        inputs[str(a.source_library.resolve())] = digest(a.source_library)

    def evaluate(job):
        case,i = job
        x,y,direction,option = case
        lo,hi = (7.0625,7.075) if direction=='temperature' else (1.9875,2.)
        z = lo+(hi-lo)*i/128
        t,q = (z,q0) if direction=='temperature' else (math.log10(T0),z)
        T,rho = 10**t,10**(q+1.5*(t-6))
        m = mixture(x,y)
        scale = m['source_mass_scale']
        request = ' '.join(map(str,m['eps']))+f'\n3 {option} -2\n'+f'{math.log(scale*rho):.17g} {math.log(T):.17g}\n'
        result = subprocess.run([str(a.probe.resolve())],input=request,text=True,capture_output=True,check=True)
        row = list(map(float,result.stdout.split()))
        if len(row)!=22 or row[0]!=0 or not all(math.isfinite(v) for v in row):
            raise ValueError('invalid source response')
        E,S,P = row[5]*scale,row[6]*scale,row[4]
        derivative = -E/T+1.5*P/rho/T if direction=='temperature' else P/rho/T
        return {'case':case,'log_coordinate':z,'input':request,'source':row,
                'phi':E/T-S,'derivative':derivative}

    with ThreadPoolExecutor(a.jobs) as pool:
        rows = list(pool.map(evaluate,[(case,i) for case in cases for i in range(129)]))
    results = []
    for k,case in enumerate(cases):
        part = rows[k*129:(k+1)*129]
        change = part[-1]['phi']-part[0]['phi']
        span = (part[-1]['log_coordinate']-part[0]['log_coordinate'])*math.log(10)
        quadrature = {}
        for stride in [2,1]:
            subset = part[::stride]
            n = len(subset)-1
            integral = span/(3*n)*sum((1 if i in [0,n] else 4 if i%2 else 2)*r['derivative']
                                      for i,r in enumerate(subset))
            quadrature[str(n)] = {'potential_change':change,'integrated_derivative':integral,
                                  'defect':change-integral,'relative_defect':(change-integral)/integral}
        results.append({'mixture':case[:2],'direction':case[2],'source_options':[3,case[3],-2],
                        'quadrature':quadrature})
    for path,expected in inputs.items():
        if digest(path)!=expected:
            raise ValueError('source input changed during audit')
    a.output.parent.mkdir(parents=True,exist_ok=True)
    raw = a.output.with_suffix('.source.json.gz')
    raw.write_bytes(gzip.compress(json.dumps(rows,allow_nan=False).encode(),mtime=0))
    report = {'scope':__doc__,'source_queries':len(rows),'input_sha256':inputs,
              'raw_source_sha256':digest(raw),'results':results}
    a.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(results,indent=2))


if __name__=='__main__':
    main()
