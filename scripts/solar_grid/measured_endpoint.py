"""Read the selected atmosphere after either a single-region or detached-region solve."""
import math
from pathlib import Path
import sys
from release import digest, read, save


def collect(root,mass):
    case=root/f'cool-m{mass:02}-1000'
    selected=read(root/'recovery'/f'{case.name}-native/selected-column.json')
    receipt=Path(selected['receipt'])
    column=read(receipt)
    spec=read(case/'spec.json')
    actual=column['spec']
    numerical={'initial_regions','resume_region_top','resume_guess'}
    def physical(s):
        return {k:v for k,v in s.items() if k not in numerical}
    assert receipt.parent==Path(selected['directory']) and actual==read(receipt.parent/'spec.json')
    assert column['quality']['usable'] and physical(actual)==physical(spec)
    assert spec['mass_Msun']==mass/100 and spec['teff']==1000
    layer=column['fresh']['at_tau']['3']
    point=[layer['T'],layer['P_bar']*1e6,layer['depth_cm_mu_He3corr']]
    assert all(math.isfinite(v) and v>0 for v in point)
    save(case/'measured-endpoint.json',dict(passes=True,endpoint=point,receipt=str(receipt),
        receipt_sha256=digest(receipt),scope='Measured 1000 K endpoint only.'))


if __name__=='__main__':
    collect(Path(sys.argv[1]).resolve(),int(sys.argv[2]))
