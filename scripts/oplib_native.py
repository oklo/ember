"""Read native MESA OPLIB values and reject missing interpolation stencils."""
import hashlib
from pathlib import Path

import numpy as np

from compare_mesa_oplib_means_v2 import interpolate, stencil


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_native(entry):
    assert digest(entry['path']) == entry['sha256']
    lines = Path(entry['path']).read_text().splitlines()
    assert '25 Elements' in lines[0] and 'Grevesse & Sauval (1998)' in lines[0]
    metadata = list(map(float, lines[2].split()))
    assert len(metadata) == 10 and metadata[:2] == [1, 37]
    assert metadata[4:] == [39, -8, 1.5, 74, 3.764, 9.065]
    rr = list(map(float, lines[5].split()))
    a = np.loadtxt(lines[6:])
    assert a.shape == (74, 40) and np.isfinite(a).all()
    assert rr == [-8 + .25*i for i in range(39)] and np.all(np.diff(a[:,0]) > 0)
    mask = a[:,1:] >= 9.99
    assert np.all(np.abs(a[:,1:][~mask]) < 8)
    missing = [dict(logT=float(a[i,0]),logR=rr[j],source_value=float(a[i,j+1]))
               for i,j in np.argwhere(mask)]
    return dict(X=metadata[2],Z=metadata[3],logT=a[:,0].tolist(),logR=rr,
                logk=a[:,1:].tolist(),mask=mask,missing=missing,source=entry)


def evaluate(table,t,r,order):
    ii,_ = stencil(table['logT'],t,order)
    jj,_ = stencil(table['logR'],r,order)
    if any(table['mask'][i,j] for i in ii for j in jj):
        raise ValueError('interpolation touches an explicitly missing native value')
    return interpolate(table,t,r,order)
