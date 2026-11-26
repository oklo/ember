"""Infer initial convective regions from a resolved radiative gap on the pressure grid."""
import json
from pathlib import Path


def split_regions(top, inward, unstable, levels):
    if not inward and unstable and unstable==sorted(set(unstable)):
        groups=[[unstable[0]]]
        for layer in unstable[1:]:
            if layer!=groups[-1][-1]+1:groups.append([])
            groups[-1].append(layer)
        if len(groups)==1 and groups[0][-1]<top-1:
            upper,deep=groups[0],top
        elif len(groups)==2 and groups[-1][-1]==top-1:
            upper,deep=groups[0],groups[1][0]
        else:return None
        if 3<=upper[0]<=upper[-1]<=deep-2 and deep<=levels-2:
            return [[upper[0],upper[-1]],[deep,levels-2]]
        return None
    if not inward or inward != list(range(inward[0], inward[-1]+1)):
        return None
    if unstable and unstable != list(range(unstable[0],unstable[-1]+1)):
        return None
    if inward[0] > top+1:
        if unstable and unstable[-1] != top-1:
            return None
        upper=[unstable[0] if unstable else top,inward[0]-2]
    elif inward[0]==top+1 and unstable and unstable[-1]<top-1:
        upper=[unstable[0],unstable[-1]]
    else:
        return None
    deep=[inward[-1],levels-2]
    if 3<=upper[0]<=upper[1] and upper[1]+2<=deep[0]<=deep[1]:
        return [upper,deep]
    return None


def neighbor_regions(receipt, target):
    if not receipt['quality']['usable']:return None
    native=receipt['spec']
    if (abs(native['teff']-target['teff'])>100 or abs(native['logg']-target['logg'])>.15
            or native['levels']!=target['levels'] or native['mass_Msun']!=target['mass_Msun']):
        return None
    raw=json.loads(Path(receipt['column']).read_bytes())
    if len(raw['pressure_bar'])!=target['levels']:return None
    cvz=[int(x) for x in str(raw['cvz_locs']).strip('[]').replace(',',' ').split()]
    if len(cvz)<6:return None
    t,b,u,v=cvz[1],cvz[2],cvz[4],cvz[5]
    if 3<=t<=b and b+2<=u<=v==target['levels']-2:
        return [[t,b],[u,v]]
    return None
