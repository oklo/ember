#!/usr/bin/env python3
"""Measure interpolation sensitivity in the printed dilute opacity spectra.

Compare integration of inverse opacity with interpolation of opacity itself.
Agreement with a source mean does not prove that its atomic frequency grid
resolves every line; this diagnostic preserves that distinction.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from audit_tops_spectral_means import digest, source
from audit_tops_electron_dispersion import RW


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('comparison',type=Path)
    p.add_argument('output',type=Path)
    a=p.parse_args()
    if a.output.exists():
        raise FileExistsError('keep completed interpolation diagnostics')
    previous=json.loads(a.comparison.read_text())
    inputs=dict(previous['input_sha256'])
    for path in (a.comparison,Path(__file__)):
        inputs[str(path.resolve())]=digest(path)
    for name,h in inputs.items():
        if digest(Path(name))!=h:
            raise ValueError('changed source input')
    cache,rows={},[]
    for row in previous['records']:
        root=row['full_spectrum_source']
        if root not in cache:
            cache[root]=source(Path(root))[0]
        key=row['temperature_keV'],row['density_atomic_g_cm3']
        data=cache[root][key]
        u=data[:,0]/key[0];keep=u<=700;u=u[keep];k=data[keep,1]
        w=u**4*np.exp(-u)/(-np.expm1(-u))**2
        values={'weighted_inverse_trapezoid':float(RW/np.trapezoid(w/k,u))}
        nodes,weights=leggauss(16)
        x=(u[1:,None]+u[:-1,None])/2+(u[1:,None]-u[:-1,None])*nodes/2
        wr=x**4*np.exp(-x)/(-np.expm1(-x))**2
        linear=(x-u[:-1,None])/(u[1:,None]-u[:-1,None])
        logarithmic=np.log(x/u[:-1,None])/np.log(u[1:]/u[:-1])[:,None]
        recipes={'linear_inverse':(1-linear)/k[:-1,None]+linear/k[1:,None],
            'log_energy_linear_inverse':(1-logarithmic)/k[:-1,None]+logarithmic/k[1:,None],
            'linear_opacity':1/((1-linear)*k[:-1,None]+linear*k[1:,None])}
        for name,inverse in recipes.items():
            values[name]=float(RW/np.sum((u[1:]-u[:-1])[:,None]/2*weights*wr*inverse))
        rows.append({**row,'alternative_integrals':values,
            'group_vs_linear_inverse_relative_difference':row['group_mean']/values['linear_inverse']-1,
            'source_grey_vs_linear_inverse_relative_difference':row['source_uncut_grey_mean']/values['linear_inverse']-1})
    for name,h in inputs.items():
        if digest(Path(name))!=h:
            raise ValueError('input changed during diagnostic')
    result={'scope':__doc__,'accepted_for_stellar_opacity':False,'states':len(rows),
        'maximum_group_vs_linear_inverse_relative_difference':max(abs(r['group_vs_linear_inverse_relative_difference']) for r in rows),
        'maximum_source_grey_vs_linear_inverse_relative_difference':max(abs(r['source_grey_vs_linear_inverse_relative_difference']) for r in rows),
        'interpretation':'The log-opacity interpolation used in the earlier monochromatic integration introduces a larger difference from the source means than inverse-opacity quadrature does. This is evidence about numerical averaging, not an independent demonstration of atomic frequency resolution.',
        'records':rows,'input_sha256':inputs}
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('records','input_sha256')}),flush=True)


if __name__=='__main__':
    main()
