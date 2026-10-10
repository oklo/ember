"""Fetch actual Gaia eligibility for retained candidates; exact integer IDs."""
from pathlib import Path
import sqlite3,json,urllib.request,urllib.parse,hashlib,time
import pandas as pd
W=Path(__file__).resolve().parent;C=W.parent/'cluster-mass-functions-oct9-v1'
URL='https://gea.esac.esa.int/tap-server/tap/sync'
with sqlite3.connect(C/'cluster_catalogues.sqlite') as con:
    o=pd.read_sql('SELECT object_id,gaia_dr3_id FROM objects WHERE gaia_dr3_id IS NOT NULL',con)
a=pd.read_csv(W/'candidate_selection_audit.csv')
ids=sorted(set(o[o.object_id.isin(a.object_id)].gaia_dr3_id))
assert all(x.isdigit() for x in ids)
frames=[];receipts=[]
for i in range(0,len(ids),1000):
    path=W/'raw'/f'candidate_gaia_{i//1000}.csv'
    query='''SELECT source_id,ra,dec,parallax,parallax_error,pmra,pmdec,
pmra_error,pmdec_error,phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,
ruwe,astrometric_excess_noise,astrometric_excess_noise_sig,
astrometric_params_solved,visibility_periods_used
FROM gaiadr3.gaia_source WHERE source_id IN ('''+','.join(ids[i:i+1000])+')'
    if not path.exists():
        t=time.monotonic()
        req=urllib.request.Request(URL,data=urllib.parse.urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=query)).encode())
        with urllib.request.urlopen(req,timeout=180) as r:data=r.read()
        assert data.startswith(b'source_id,')
        path.write_bytes(data)
        path.with_suffix('.json').write_text(json.dumps(dict(query=query,url=URL,sha256=hashlib.sha256(data).hexdigest(),wall_seconds=time.monotonic()-t),indent=2)+'\n')
    frames.append(pd.read_csv(path,dtype={'source_id':str}))
    print(i,len(frames[-1]),flush=True)
allrows=pd.concat(frames,ignore_index=True)
assert allrows.source_id.is_unique
allrows.to_csv(W/'candidate_gaia.csv',index=False)
(W/'candidate_gaia_receipt.json').write_text(json.dumps(dict(requested=len(ids),returned=len(allrows),missing=sorted(set(ids)-set(allrows.source_id))),indent=2)+'\n')
