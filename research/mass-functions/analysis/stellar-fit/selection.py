"""Gaia parent-catalogue detection, distinct from cluster-member selection.

The empirical function/parameters are from GaiaUnlimited surveyTCG.py,
Copyright (c) 2022 GaiaUnlimited, BSD-3-Clause; full license in raw/.
M10 map: Cantat-Gaudin et al. (2023), Zenodo 8063930, CC-BY-4.0.
"""
from pathlib import Path
import ast,hashlib,json
import numpy as np
import h5py
from astropy_healpix import lonlat_to_healpix
from astropy import units as u

W=Path(__file__).resolve().parent
# Load only the two published numerical functions, avoiding astroquery and
# their automatic network/cache setup. Retain the downloaded source verbatim.
tree=ast.parse((W/'raw/surveyTCG.py').read_text())
funcs=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['sigmoid','m10_to_completeness']]
assert len(funcs)==2
scope={'np':np}
exec(compile(ast.Module(body=funcs,type_ignores=[]),'published_Gaia_completeness','exec'),scope)
with h5py.File(W/'raw/allsky_M10_hpx7.hdf5','r') as f: data=f['data'][:]
assert np.all(data[:,0]==7) and len(data)==12*128**2
M10=np.empty(len(data));M10[data[:,1].astype(int)]=data[:,2]

def m10(ra,dec):
    index=lonlat_to_healpix(np.asarray(ra)*u.deg,np.asarray(dec)*u.deg,128,order='nested')
    return M10[index]

def completeness(g,local_m10):
    g=np.asarray(g);p=scope['m10_to_completeness'](g,np.broadcast_to(local_m10,g.shape))
    return np.clip(p,0,1)

if __name__=='__main__':
    g=np.linspace(12,25,101);c=completeness(g,20.9)
    assert np.isfinite(c).all() and ((c>=0)&(c<=1)).all() and np.all(np.diff(c)<=0)
    assert c[0]>.999 and c[-1]<.001
    receipt=dict(source='https://github.com/gaia-unlimited/gaiaunlimited/blob/main/src/gaiaunlimited/selectionfunctions/surveyTCG.py',
        data='https://zenodo.org/records/8063930',
        sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (W/'raw').iterdir() if p.is_file()},
        scope='Gaia DR3 source-catalogue inclusion, NOT astrometry-quality or cluster-membership recovery.',
        monotonicity_and_probability_checks=True)
    (W/'selection_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    for cluster,ra,dec in [('hyades',66.75,15.87),('pleiades',56.75,24.12)]:
        m=float(m10(ra,dec));print(cluster,'M10',m,'G=19,20,21,22',completeness(np.array([19.,20.,21.,22.]),m))
