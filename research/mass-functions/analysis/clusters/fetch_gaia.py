"""Fetch DR3 photometry/errors for selected catalogue objects in small batches."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json, hashlib, time, urllib.request, urllib.parse
import pandas as pd

W=Path(__file__).resolve().parent
frame=pd.read_csv(W/'catalogues/combined_objects.csv',dtype={'gaia_dr3_id':str},low_memory=False)
ids=sorted(set(frame.loc[frame.selected_by_any & frame.gaia_dr3_id.notna(),'gaia_dr3_id']))
dest=W/'raw/gaia_dr3';dest.mkdir(exist_ok=True)
columns='source_id,ra,dec,parallax,parallax_error,pmra,pmdec,phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,phot_g_mean_flux_over_error,phot_bp_mean_flux_over_error,phot_rp_mean_flux_over_error,ruwe,phot_bp_rp_excess_factor'
def run(block):
    query=f'SELECT {columns} FROM gaiadr3.gaia_source WHERE source_id IN ('+','.join(block)+')'
    key=hashlib.sha256(query.encode()).hexdigest()[:16];p=dest/f'{key}.csv'
    if not p.exists():
        url='https://gea.esac.esa.int/tap-server/tap/sync?'+urllib.parse.urlencode(
            dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=query))
        for attempt in range(3):
            try:
                payload=urllib.request.urlopen(url,timeout=90).read()
                assert payload.startswith(b'source_id,')
                p.write_bytes(payload);break
            except Exception:
                if attempt==2:raise
                time.sleep(2)
        p.with_suffix('.json').write_text(json.dumps(dict(url=url,query=query,
            sha256=hashlib.sha256(payload).hexdigest(),requested_ids=block),indent=2)+'\n')
    data=pd.read_csv(p,dtype={'source_id':str})
    assert set(data.source_id)==set(block),(p,len(data),len(block))
    print('Gaia batch',key,len(data),flush=True)
    return data
# Reuse previously verified batches even if adding a catalogue changes the
# sorting of the requested IDs. Only query the newly required sources.
previous=W/'catalogues/gaia_dr3_photometry.csv'
cached=pd.read_csv(previous,dtype={'source_id':str}) if previous.exists() else pd.DataFrame(columns=['source_id'])
cached=cached[cached.source_id.isin(ids)]
missing=sorted(set(ids)-set(cached.source_id))
with ThreadPoolExecutor(max_workers=3) as pool:
    tables=[cached]+list(pool.map(run,[missing[i:i+200] for i in range(0,len(missing),200)]))
out=pd.concat(tables,ignore_index=True)
assert out.source_id.nunique()==len(out)==len(ids)
out.to_csv(W/'catalogues/gaia_dr3_photometry.csv',index=False)
print('Complete:',len(out),'distinct Gaia sources',flush=True)
