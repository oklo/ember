#!/usr/bin/env python3
"""Direct FreeEOS sensitivity to replacing the He4 metal proxy by GS98 metals.

This is an audit, not a new EOS grid. Potassium (4.56e-6 by mass) is absent
from FreeEOS's element list and omitted. Integer masses use the baryon basis.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('probe',type=Path);ap.add_argument('output',type=Path);a=ap.parse_args()
    mixture=json.loads(Path('data/opacity/sources/tops_gs98_x070_z020.request.json').read_text())['mixture'].split()
    metals={mixture[i+1]:float(mixture[i]) for i in range(4,len(mixture),2)}
    elements=['c','n','o','ne','na','mg','al','si','p','s','cl','ar','ca','ti','cr','mn','fe','ni']
    masses=[12,14,16,20,23,24,27,28,31,32,35,40,40,48,52,55,56,58]
    rows=[]
    for x,y3 in [(.7,0),(.6,.08),(.45,.08)]:
        for t,r in [(4100,1.4e-4),(6000,.001),(3e4,.03),(1e6,3),(4.6e6,390)]:
            states=[]
            for metal in [False,True]:
                eps=[x,y3/3+((.98 if metal else 1)-x-y3)/4]+([metals[e]/m for e,m in zip(elements,masses,strict=True)] if metal else [0]*18)
                request=' '.join(map(str,eps))+'\n3 1 -2\n'+f'{math.log(r)} {math.log(t)}\n'
                response=subprocess.run([str(a.probe.resolve())],input=request,capture_output=True,text=True,check=True,timeout=60)
                v=list(map(float,response.stdout.split()))
                if len(v)!=22 or v[0]!=0 or not all(math.isfinite(z) for z in v) or abs(v[2]/r-1)>1e-10 or abs(v[3]/t-1)>1e-10:raise ValueError('invalid source response')
                states.append(v)
            rows.append({'X':x,'Y3':y3,'T':t,'rho':r,'metal_to_proxy_P':states[1][4]/states[0][4],
                'metal_to_proxy_E':states[1][5]/states[0][5],'metal_to_proxy_cp':states[1][13]/states[0][13],'metal_to_proxy_ad':states[1][14]/states[0][14]})
    a.output.write_text(json.dumps({'description':'Direct FreeEOS GS98 metal sensitivity; K omitted (4.56e-6), integer masses, no assertion of actual EOS consistency with carried AAG21 lumped Zrest','rows':rows},indent=2)+'\n')
if __name__=='__main__':main()
