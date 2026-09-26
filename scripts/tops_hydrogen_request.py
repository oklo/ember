"""Request another H/He ratio while preserving a verified metal mixture."""
from decimal import Decimal
import math


def hydrogen_request(plan, baseline, request):
    expected = {k:baseline[k] for k in ('X','Z','metals')}
    if 'hydrogen_atomic_mass_fraction' not in plan:
        return dict(request), expected
    x = plan['hydrogen_atomic_mass_fraction']
    if isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x):
        raise ValueError('hydrogen fraction must be a finite number')
    dx,dz = Decimal(str(x)),Decimal(str(baseline['Z']))
    if dx < 0 or dx+dz >= 1 or dx != dx.quantize(Decimal('1e-10')):
        raise ValueError('unrepresentable hydrogen or absent helium')
    if request.get('fractype') != 'mass' or request.get('isotope') != 'noniso':
        raise ValueError('hydrogen override requires elemental atomic mass fractions')
    parts = request['mixture'].split()
    if len(parts)%2 or parts[1:4:2] != ['h','he']:
        raise ValueError('unexpected source mixture layout')
    metals = {parts[i+1].capitalize():float(parts[i]) for i in range(4,len(parts),2)}
    if len(metals) != (len(parts)-4)//2 or metals != baseline['metals']:
        raise ValueError('source metal inventory differs')
    if (abs(float(parts[0])-baseline['X']) > 5e-11 or
            abs(float(parts[2])-(1-baseline['X']-baseline['Z'])) > 5e-11 or
            abs(sum(metals.values())-baseline['Z']) > 2e-9):
        raise ValueError('source atomic fractions do not close')
    result = dict(request)
    result['mixture'] = f'{dx:.10f} h {Decimal(1)-dx-dz:.10f} he ' + ' '.join(parts[4:])
    return result, {**expected,'X':float(dx)}
