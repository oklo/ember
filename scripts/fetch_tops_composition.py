#!/usr/bin/env python3
"""Archive TOPS mixtures sequentially, validating each result before the next.

The server prepares a calculation at /submit and serves it at /results.
No service requests are needed for normal builds or reproducible import.
"""
import argparse
from html.parser import HTMLParser
import json
from pathlib import Path
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
import secrets


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
    def post(endpoint, data):
        req=urllib.request.Request('https://aphysics2.lanl.gov/'+endpoint,
                                  data=urllib.parse.urlencode(data).encode())
        return urllib.request.urlopen(req,context=ctx,timeout=180).read().decode()
    for x in a.hydrogen:
        if not 0<=x<=1-a.metallicity: raise ValueError('invalid hydrogen abundance')
        label=f'tops_gs98_x{round(x*100):03d}_z{round(a.metallicity*1000):03d}'
        # Results are keyed by user ID. Isolate every requested mixture;
        # reusing an ID can return a completed, stale calculation.
        p=baseline.copy(); p['userid']='e'+''.join(secrets.choice('abcdefghijklmnopqrstuvwxyz') for _ in range(2))
        p['mixture']=f'{x:.10f} h {1-x-a.metallicity:.10f} he '+metals
        submitted=post('submit',p)
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
        (a.output/(label+'.request.json')).write_text(json.dumps(p,indent=2)+'\n')
        (a.output/(label+'.txt')).write_text(text)
        print(f'archived verified X={x:g}',flush=True)


if __name__=='__main__':main()
