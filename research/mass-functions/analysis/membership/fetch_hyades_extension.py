from pathlib import Path
import urllib.request,urllib.parse,hashlib,json,time
w=Path(__file__).resolve().parent
query=(w/'raw/hyades_gaia_cp_query.adql').read_text().replace('COS(RADIANS(l))<0','COS(RADIANS(l))>=0')
p=w/'raw/hyades_gaia_cp_positive_x.csv'
url='https://gea.esac.esa.int/tap-server/tap/sync'
if not p.exists():
 start=time.monotonic();req=urllib.request.Request(url,data=urllib.parse.urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=query)).encode())
 with urllib.request.urlopen(req,timeout=240) as response:data=response.read()
 assert data.startswith(b'source_id,') and data.count(b'\n')-1<1000000
 p.write_bytes(data);p.with_suffix('.json').write_text(json.dumps(dict(query=query,url=url,sha256=hashlib.sha256(data).hexdigest(),rows=data.count(b'\n')-1,wall_seconds=time.monotonic()-start,reason='The published Hyades CP catalogue contains members at positive X as well as negative X. Retain both hemispheres for reproduction; no hand-drawn rejection masks are claimed.'),indent=2)+'\n')
 print('Fetched',data.count(b'\n')-1,'additional parent sources',flush=True)
