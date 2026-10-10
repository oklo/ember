"""Fill a photometry gap for the matched-aperture check; no mass changes."""
from pathlib import Path
from contextlib import closing
import hashlib,io,json,time,urllib.parse,urllib.request,sqlite3
import pandas as pd
W=Path(__file__).resolve().parent;V1=W.parent/'cluster-dissolution-oct9-v1'
d=pd.read_csv(V1/'hyades_spatial_selection.csv',dtype={'gaia_dr3_id':str})
with closing(sqlite3.connect(W.parent/'cluster-mass-functions-oct9-v1/cluster_catalogues.sqlite')) as con:
    exact=pd.read_sql('SELECT object_id,gaia_dr3_id FROM objects',con).set_index('object_id').gaia_dr3_id
d['gaia_dr3_id']=d.object_id.map(exact)
ids=sorted(set(d.loc[((d.p_mass_reference_star>.001)|(d.p_mass_BD>.001))&d.Jmag.isna(),'gaia_dr3_id'].dropna()))
url='https://gea.esac.esa.int/tap-server/tap/sync'
parts=[];receipt=[]
for i in range(0,len(ids),100):
    batch=ids[i:i+100];assert all(x.isdigit() for x in batch)
    query='SELECT g.source_id, t.designation, t.j_m, t.ks_m, t.j_msigcom, t.ks_msigcom FROM gaiadr3.tmass_psc_xsc_best_neighbour AS g JOIN gaiadr1.tmass_original_valid AS t ON g.original_ext_source_id=t.designation WHERE g.source_id IN ('+','.join(batch)+')'
    path=W/f'hyades_2mass_raw_{hashlib.sha256(query.encode()).hexdigest()[:12]}.csv'
    if not path.exists():
        req=urllib.request.Request(url,data=urllib.parse.urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=query)).encode(),headers={'User-Agent':'Ember-cluster-research'})
        with urllib.request.urlopen(req,timeout=60) as response:data=response.read()
        assert b'j_m' in data and b'VOTABLE' not in data,data[:300]
        path.write_bytes(data)
    table=pd.read_csv(path,dtype={'source_id':str});parts.append(table)
    receipt.append(dict(query=query,url=url,file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    print(i,len(table),flush=True)
pd.concat(parts,ignore_index=True).to_csv(W/'hyades_2mass.csv',index=False)
(W/'hyades_2mass_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
