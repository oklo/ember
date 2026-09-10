"""GS98 element numbers for the baryonic FreeEOS source family.

Uses exactly the atmosphere's fixed GS98 element pattern and representative
integer masses. FreeEOS has no potassium entry; that trace omission is
recorded explicitly, rather than folding all metals into helium.
"""
import hashlib
import json
from pathlib import Path
from generate_nongrey_grid import composition
from stellar_composition import interior_composition

ROOT=Path(__file__).resolve().parents[1]
ELEMENTS=['H','He','C','N','O','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','Ca','Ti','Cr','Mn','Fe','Ni']


def mixture(hydrogen,helium3):
    carried=interior_composition(hydrogen)
    carried[1]=helium3;carried[2]-=helium3
    if min(carried)<0:raise ValueError('invalid baryonic composition')
    abundance,_=composition(hydrogen,helium3,carried[3:])
    path=ROOT/'data/atmosphere/sources/synple-elements.json'
    elements=json.loads(path.read_text())
    index={s.lower():i for i,s in enumerate(elements['symbol'])}
    number=[hydrogen*abundance[index[s.lower()]] for s in ELEMENTS]
    # Convert to a normalized atomic source gram, then convert every source
    # result back to the conserved baryonic gram. Helium uses source He4.
    weights=[elements['mass'][index[s.lower()]] for s in ELEMENTS]
    weights[0]=1.00782503;weights[1]=4.00260325
    scale=sum(n*w for n,w in zip(number,weights,strict=True))
    eps=[n/scale for n in number]
    k=index['k']
    return {'composition_basis':'baryon_mass','metal_mixture':'GS98',
            'hydrogen':hydrogen,'helium3':helium3,'composition':carried,
            'elements':ELEMENTS,'eps':eps,'source_mass_scale':scale,
            'potassium_number_per_baryon_mass':hydrogen*abundance[k],
            'potassium_baryonic_mass_fraction':hydrogen*abundance[k]*round(elements['mass'][k]),
            'element_metadata_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'metal_request_sha256':hashlib.sha256((ROOT/'data/opacity/sources/tops_gs98_x070_z020.request.json').read_bytes()).hexdigest(),
            'approximation':'Shared atmosphere GS98 element inventory; FreeEOS omits potassium; representative metal isotopes; source He4 electronic physics and analytic helium isotope entropy; inert carried metal slots retain their existing labels'}
