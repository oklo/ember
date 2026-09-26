"""A bounded extension that preserves every old collision-table query exactly."""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import BPoly, CubicSpline


class ExtendedIonTable:
    def __init__(self, specification):
        data=json.loads(Path(specification).read_text())
        p=Path(data['base_table'])
        if hashlib.sha256(p.read_bytes()).hexdigest()!=data['base_table_sha256']:
            raise ValueError('base collision table changed')
        base=json.loads(p.read_text())
        x=np.asarray(base['log10_strength']);y=np.log(base['dimensionless_integrals'])
        self.base=CubicSpline(x,y,axis=0,extrapolate=False)
        xn=np.asarray(data['additional_log10_strength']);yn=np.log(data['additional_dimensionless_integrals'])
        if (xn.ndim!=1 or len(xn)<3 or yn.shape!=(len(xn),4) or
                not np.all(np.isfinite(xn)) or not np.all(np.isfinite(yn)) or
                xn[0]<=x[-1] or np.any(np.diff(xn)<=0)):
            raise ValueError('ordered finite new collision data required')
        support=CubicSpline(np.r_[x[-6:],xn],np.vstack([y[-6:],yn]),axis=0,extrapolate=False)
        xe=np.r_[x[-1],xn]
        derivatives=[np.array([self.base(x[-1],order) for order in range(3)])]
        derivatives.extend(np.array([support(value,order) for order in range(3)]) for value in xn)
        self.extension=BPoly.from_derivatives(xe,derivatives,extrapolate=False)
        self.minimum,self.join,self.maximum=float(x[0]),float(x[-1]),float(xn[-1])

    def log_moments(self, log10_strength, derivative=0):
        x=np.asarray(log10_strength,dtype=float)
        flat=x.ravel();result=np.empty((len(flat),4))
        old=flat<=self.join
        result[old]=self.base(flat[old],derivative)
        result[~old]=self.extension(flat[~old],derivative)
        return result.reshape(x.shape+(4,))

    def moments(self, strength):
        strength=np.asarray(strength,dtype=float)
        if np.any(strength<=0) or not np.all(np.isfinite(strength)):
            raise ValueError('positive finite collision strengths required')
        return np.exp(self.log_moments(np.log10(strength)))
