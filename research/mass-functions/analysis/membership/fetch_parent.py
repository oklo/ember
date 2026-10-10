"""Fetch the published CP parent sample, retaining the real Galactic background.

The remote query uses only the published necessary cuts. The broad Z window
is narrowed for the Hyades in the local selection comparison. No magnitude cut
is imposed beyond Gaia's catalogue and the original required photometry.
"""
from pathlib import Path
import urllib.request,urllib.parse,json,hashlib,time
W=Path(__file__).resolve().parent
URL='https://gea.esac.esa.int/tap-server/tap/sync'
PARAMS={'hyades':(97.3,6.9,-17.),'pleiades':(89.8,-44.7,-54.)}
for cluster,(ac,dc,zc) in PARAMS.items():
    a=f'SIN(RADIANS({ac}-ra))'
    b=f'(COS(RADIANS(dec))*TAN(RADIANS({dc}))-SIN(RADIANS(dec))*COS(RADIANS({ac}-ra)))'
    vp=f'4.74047*(-({b})*pmra+({a})*pmdec)/(parallax*SQRT(POWER({a},2)+POWER({b},2)))'
    vmean=4.74047*((104.136**2+28.732**2)**.5)/21.2 if cluster=='hyades' else 4.74047*((19.955**2+45.457**2)**.5)/7.38
    vparallel=f'4.74047*(({a})*pmra+({b})*pmdec)/(parallax*SQRT(POWER({a},2)+POWER({b},2)))'
    query=f'''SELECT TOP 1000000 source_id,ra,dec,l,b,parallax,parallax_error,pmra,pmdec,pmra_error,pmdec_error,parallax_pmra_corr,parallax_pmdec_corr,pmra_pmdec_corr,phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,ruwe
FROM gaiadr3.gaia_source
WHERE parallax>2 AND pmra IS NOT NULL AND pmdec IS NOT NULL
AND 1000*parallax_error/POWER(parallax,2)<10
AND (astrometric_excess_noise<1 OR astrometric_excess_noise_sig<2) AND ruwe<1.4
AND phot_g_mean_mag IS NOT NULL AND phot_bp_mean_mag IS NOT NULL AND phot_rp_mean_mag IS NOT NULL
AND COS(RADIANS(l))<0
AND ABS(1000*SIN(RADIANS(b))/parallax-({zc}))<55
AND ABS({vp})<10.1
AND ABS(({vparallel})-({vmean}))<50.1'''
    path=W/'raw'/f'{cluster}_gaia_cp_parent.csv'
    (W/'raw'/f'{cluster}_gaia_cp_query.adql').write_text(query+'\n')
    if path.exists():continue
    start=time.monotonic()
    request=urllib.request.Request(URL,data=urllib.parse.urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=query)).encode(),headers={'User-Agent':'Ember-cluster-research'})
    try:
        with urllib.request.urlopen(request,timeout=240) as response:data=response.read()
        if not data.startswith(b'source_id,'):raise RuntimeError(data[:1000])
        lines=data.count(b'\n')-1
        if lines>=1000000:raise RuntimeError(f'Query hit row cap: {lines}')
        path.write_bytes(data)
        receipt=dict(cluster=cluster,url=URL,query=query,rows=lines,bytes=len(data),wall_seconds=time.monotonic()-start,sha256=hashlib.sha256(data).hexdigest(),scope='Gaia parent sample for CP selection; NOT selected cluster members.')
        path.with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n')
        print(cluster,lines,len(data),flush=True)
    except Exception as exc:
        print(cluster,type(exc).__name__,str(exc),flush=True)
        if hasattr(exc,'read'):(W/'raw'/f'{cluster}_query_error.txt').write_bytes(exc.read())
