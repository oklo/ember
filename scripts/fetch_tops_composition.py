#!/usr/bin/env python3
"""Archive TOPS mixtures sequentially, validating each result before the next.

The server prepares a calculation at /submit and serves it at /results.
No service requests are needed for normal builds or reproducible import.
"""
import argparse
from decimal import Decimal
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
import secrets

from import_tops_composition import read as read_source


def fraction_label(value, scale):
    """Exact decimal labels: .125 is x012p5, never the rounded x012.

    Old archives keep their original names and explicit request mixtures.
    Refinement must not overwrite a nearby composition with the same label.
    """
    number = Decimal(str(value))
    if not number.is_finite() or number < 0 or number != number.quantize(Decimal('1e-10')):
        raise ValueError('composition must be finite, nonnegative and representable in the source request')
    whole, _, fractional = format(number * scale, 'f').partition('.')
    fractional = fractional.rstrip('0')
    return f'{int(whole):03d}' + ('p' + fractional if fractional else '')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def cached_result(directory, label, x, z, metals):
    """Only a complete, pinned, revalidated result can skip a service call."""
    request = directory / (label + '.request.json')
    source = directory / (label + '.txt')
    receipt = directory / (label + '.receipt.json')
    if not any(p.exists() for p in [request, source, receipt]):
        return False
    if not all(p.is_file() for p in [request, source, receipt]):
        raise ValueError('incomplete existing source archive; preserve it and use a new output directory')
    record = json.loads(receipt.read_text())
    expected = {'X': x, 'Z': z, 'metals': metals, 'file': source.name,
                'request': request.name}
    if any(record.get(k) != v for k, v in expected.items()):
        raise ValueError('cached source composition differs')
    if digest(request) != record['request_sha256']:
        raise ValueError('cached request changed')
    read_source(source, record)
    return True


class Form(HTMLParser):
    def __init__(self):
        super().__init__(); self.active=False; self.done=False; self.data={}
    def handle_starttag(self, tag, attrs):
        a=dict(attrs)
        if tag=='form' and a.get('action')=='./results' and not self.done:
            self.active=True
        if tag=='input' and self.active and 'name' in a:
            self.data[a['name']]=a.get('value','')
    def handle_endtag(self, tag):
        if tag=='form' and self.active: self.active=False; self.done=True


class Text(HTMLParser):
    def __init__(self): super().__init__(); self.parts=[]
    def handle_data(self, data): self.parts.append(data)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('output',type=Path)
    ap.add_argument('hydrogen',type=float,nargs='+')
    ap.add_argument('--metallicity',type=float,default=.02)
    ap.add_argument('--insecure',action='store_true',help='request-specific workaround for local certificate-chain failure')
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    ctx=ssl._create_unverified_context() if a.insecure else ssl.create_default_context()
    baseline=json.loads(Path('data/opacity/sources/tops_gs98_x070_z020.request.json').read_text())
    if not 0<a.metallicity<.1: raise ValueError('invalid metallicity')
    metal_tokens=baseline['mixture'].split(' he ',1)[1].split()
    metals=' '.join(f'{float(metal_tokens[i])*a.metallicity/.02:.10f} {metal_tokens[i+1]}'
                    for i in range(0,len(metal_tokens),2))
    tokens=metals.split()
    metal_fractions={tokens[i+1].capitalize():float(tokens[i]) for i in range(0,len(tokens),2)}
    def post(endpoint, data):
        req=urllib.request.Request('https://aphysics2.lanl.gov/'+endpoint,
                                  data=urllib.parse.urlencode(data).encode())
        return urllib.request.urlopen(req,context=ctx,timeout=180).read().decode()
    for x in a.hydrogen:
        if not 0<=x<=1-a.metallicity: raise ValueError('invalid hydrogen abundance')
        label=f'tops_gs98_x{fraction_label(x,100)}_z{fraction_label(a.metallicity,1000)}'
        if cached_result(a.output,label,x,a.metallicity,metal_fractions):
            print(f'reused independently revalidated source archive X={x:g}',flush=True)
            continue
        # Results are keyed by user ID. Isolate every requested mixture;
        # reusing an ID can return a completed, stale calculation.
        p=baseline.copy(); p['userid']='e'+''.join(secrets.choice('abcdefghijklmnopqrstuvwxyz') for _ in range(2))
        p['mixture']=f'{x:.10f} h {1-x-a.metallicity:.10f} he '+metals
        # The submission gateway can time out after accepting a calculation.
        # Keep the same request/ID and use bounded backoff; completed local
        # receipts ensure a later plan restart does not duplicate earlier jobs.
        for attempt in range(3):
            try:
                submitted=post('submit',p)
                break
            except (TimeoutError,urllib.error.HTTPError) as e:
                if isinstance(e,urllib.error.HTTPError) and e.code!=504:raise
                if attempt==2:raise
                print(f'TOPS submission timed out for X={x:g}, retry {attempt+1}',flush=True)
                time.sleep(5*3**attempt)
        form=Form();form.feed(submitted)
        if not form.done or form.data.get('output')!='tabcol':
            raise ValueError('unexpected TOPS response form')
        # A calculation can outlive the gateway, or /results can initially
        # serve the previous calculation. Never archive that stale mixture.
        for attempt in range(6):
            try:
                html=post('results',form.data)
                parser=Text();parser.feed(html)
                text='\n'.join(v.strip() for v in ''.join(parser.parts).replace('\u00a0',' ').splitlines() if v.strip())+'\n'
                rows=text.split('No. Fraction Mass Fraction  At. No.  Chem. Sym.  Mat ID.\n')[1].split('Temperature grid')[0]
                elements={r.split()[3]:float(r.split()[1]) for r in rows.splitlines() if r.strip()}
                if len(elements)==21 and abs(elements['H']-x)<1e-6 and abs(elements['He']-(1-a.metallicity-x))<1e-6:
                    break
                print(f'rejecting stale TOPS result for X={x:g}, retry {attempt+1}',flush=True)
            except (TimeoutError,urllib.error.HTTPError) as e:
                if isinstance(e,urllib.error.HTTPError) and e.code!=504:raise
                print(f'TOPS request timed out for X={x:g}, retry {attempt+1}',flush=True)
            time.sleep(5)
        else:raise ValueError('TOPS did not return the requested mixture')
        request=a.output/(label+'.request.json');source=a.output/(label+'.txt')
        # Exclusive creation preserves prior responses even if two fetchers
        # accidentally target the same directory. A partial pair is rejected
        # on the next invocation, rather than silently being replaced.
        with request.open('x') as f:f.write(json.dumps(p,indent=2)+'\n')
        with source.open('x') as f:f.write(text)
        record={'X':x,'Z':a.metallicity,'metals':metal_fractions,
                'file':source.name,'sha256':digest(source),
                'request':request.name,'request_sha256':digest(request)}
        read_source(source,record)
        with (a.output/(label+'.receipt.json')).open('x') as f:
            f.write(json.dumps(record,indent=2)+'\n')
        print(f'archived verified X={x:g}',flush=True)


if __name__=='__main__':main()
