"""Thermal moments using finite, split speed integrals.

The scattering calculation is the unchanged, independently orbit-checked v1
implementation. Integrate in thermal speed rather than a global generalized
Laguerre rule, resolving the low-speed contribution explicitly.
"""
import math
import numpy as np
from yukawa_scattering_control import cross_sections,gauss


def collision_integrals(g, *, thermal_order=32, angle_order=128,
                        radius_order=256, thermal_max=8.):
    if not math.isfinite(g) or g<=0:
        raise ValueError('positive collision strength required')
    if thermal_max<6:
        raise ValueError('thermal tail insufficient')
    split=sorted({0.,min(math.sqrt(g),thermal_max),1.,2.,4.,6.,thermal_max})
    node,weight=gauss(thermal_order)
    result=np.zeros(4)
    for lo,hi in zip(split,split[1:]):
        v=lo+(hi-lo)*node
        phi=np.array([cross_sections(x/math.sqrt(g),angle_order=angle_order,
                                     radius_order=radius_order) for x in v])
        base=(hi-lo)*weight*np.exp(-v*v)*v**5/(g*g)
        result += [base@phi[:,0],(base*v*v)@phi[:,0],
                   (base*v**4)@phi[:,0],(base*v*v)@phi[:,1]]
    if not np.all(np.isfinite(result)) or np.any(result<=0):
        raise ValueError('invalid collision integral')
    return dict(zip(((1,1),(1,2),(1,3),(2,2)),map(float,result)))
