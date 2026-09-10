"""Audit actual cool layers and the explicit fully vaporized hot continuation.

FastChem is independently queried at actual T/P only through 6000 K. For
hotter native atmosphere layers, test the 6000 K join at each actual pressure;
never relabel those chemistry queries as calculations at the hotter state.
"""
import numpy as np
from audit_condensates import closure

HOT_JOIN_TEMPERATURE=6000.


def audit_profile(chemistry,profile,numbers,masses,hot_join=None):
    cool=[i for i,r in enumerate(profile) if r[3]<=HOT_JOIN_TEMPERATURE]
    hot=[i for i,r in enumerate(profile) if r[3]>HOT_JOIN_TEMPERATURE]
    if not cool:raise ValueError('this atmosphere audit requires at least one layer in the tested chemistry domain')
    audit=closure(chemistry,[profile[i] for i in cool],numbers,masses)
    condensed=np.zeros(len(profile),dtype=bool)
    condensed[cool]=[max(r['condensed'])/r['total_element_density']>1e-14 for r in chemistry['rows']]
    if hot:
        if hot_join is None:raise ValueError('hot atmosphere requires independent vaporized-join chemistry')
        if hot_join['mode']!='equilibrium':raise ValueError('hot join requires equilibrium chemistry')
        queries=[]
        for i in hot:
            row=profile[i].copy();row[3]=HOT_JOIN_TEMPERATURE;queries.append(row)
        joined=closure(hot_join,queries,numbers,masses)
        maximum=max(max(r['condensed'])/r['total_element_density'] for r in hot_join['rows'])
        if maximum>1e-14:raise ValueError('condensates at 6000 K invalidate the fully vaporized continuation')
        audit['vaporized_continuation']={
            'layers':len(hot),'actual_temperature_range_K':[min(profile[i][3] for i in hot),max(profile[i][3] for i in hot)],
            'pressure_range_bar':[min(profile[i][6]/1e6 for i in hot),max(profile[i][6]/1e6 for i in hot)],
            'independently_queried_temperature_K':HOT_JOIN_TEMPERATURE,
            'maximum_join_condensate_particles_per_nucleus':maximum,
            'join_closure':{k:joined[k] for k in ['element_closure','pressure_closure','nuclei_closure']},
            'scope':'Actual T/P chemistry is audited through6000K. Hotter layers use the native fully vaporized prescription; independent equilibrium is checked at6000K and each actual pressure, not extrapolated to the hotter temperature.'}
    elif hot_join is not None:raise ValueError('unexpected hot-join chemistry for a wholly cool atmosphere')
    return audit,condensed
