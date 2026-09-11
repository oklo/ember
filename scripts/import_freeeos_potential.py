#!/usr/bin/env python3
"""Build C2 Helmholtz interpolation jets from direct FreeEOS 3.0 evaluations.

P, E, S and their first responses must satisfy the local first law before
a node is accepted. Small source-fit discontinuities across joins are
regularized by this C2 potential representation, with the resulting
source-response differences and grid sensitivity audited in docs/FREEEOS.md.
Higher mixed derivatives complete the biquintic representation using
fourth-order centered differences of the source Hessian. The outer two
source rows on every edge are omitted, providing the complete stencil.
"""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
from stellar_composition import interior_composition


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('source', type=Path)
    ap.add_argument('output', type=Path)
    ap.add_argument('--hydrogen', type=float, default=.7)
    a = ap.parse_args()
    original = a.source.read_bytes()
    src = json.loads(gzip.decompress(original))
    if src['version'] != 'FreeEOS 3.0.0' or src['options'] not in ([3, 1, -2], [3, 223, -2]):
        raise ValueError('expected FreeEOS 3.0 EOS1 or numerical-electron source')
    numerical_electrons = src['options'] == [3, 223, -2]
    if not math.isfinite(a.hydrogen) or not 0 <= a.hydrogen <= .98:
        raise ValueError('invalid hydrogen abundance')
    baryonic = src.get('composition_basis') == 'baryon_mass'
    if numerical_electrons and not baryonic:
        raise ValueError('numerical-electron source requires the explicit baryonic mixture')
    helium = .3 if a.hydrogen == .7 else 1-a.hydrogen
    expected_eps=[a.hydrogen/1.00782503,helium/4.00260325]+[0.]*18
    if baryonic:
        from metal_eos_composition import mixture
        expected = mixture(src['hydrogen'],src['helium3'])
        for key,value in expected.items():
            if src.get(key) != value:
                raise ValueError(f'metal source metadata mismatch: {key}')
        expected_eps = expected['eps']
    if len(src['eps'])!=20 or any(abs(x-y)>1e-14 for x,y in zip(src['eps'],expected_eps)):
        raise ValueError('source composition differs from the specified H/He mixture')
    lt, lq, rows = src['logT'], src['logQ'], src['data']
    nt, nq = len(lt), len(lq)
    if len(rows) != nt*nq or nt < 5 or nq < 5:
        raise ValueError('bad grid dimensions')
    dx, dy = (lt[1]-lt[0])*math.log(10), (lq[1]-lq[0])*math.log(10)
    for axis in [lt, lq]:
        if any(abs((axis[i]-axis[i-1])/(axis[1]-axis[0])-1) > 1e-10 for i in range(1,len(axis))):
            raise ValueError('requires a uniform source grid')
    phi = [[0.]*9 for _ in rows]
    good = []
    # FreeEOS derives radiation constants from CODATA inputs, rather than a
    # rounded sigma_SB. Remove its radiation here; ember adds its own once.
    arad = 8*math.pi**5*1.380649e-16**4/(15*6.62607015e-27**3*2.99792458e10**3)
    for k,r in enumerate(rows):
        if not all(math.isfinite(v) for v in r):
            raise ValueError('nonfinite source value')
        info, _, rho, T, P, E, S, chir, chit, Er, Et, Sr, St, cp, ad, delta, *_ = r
        if info != 0 and baryonic:
            # Never interpolate through a failed source state. Placeholder
            # jets are written only under a mask expanded by the full stencil.
            good.append(False)
            continue
        if info != 0:
            raise ValueError('nonconverged source state')
        expected_t=lt[k//nq]*math.log(10)
        expected_r=(lq[k%nq]+1.5*(lt[k//nq]-6))*math.log(10)
        if abs(math.log(T)-expected_t)>1e-10 or abs(math.log(rho)-expected_r)>1e-9:
            raise ValueError('source failed to match requested temperature/density')
        defects = [rho*Er/P+chit-1, T*St/Et-1, rho*T*Sr/P+chit]
        if max(abs(v) for v in defects) > 1e-7:
            raise ValueError(f'inconsistent source state {k}: {defects}')
        # Option 223 retains the EOS1 material choices but computes the
        # electron integrals numerically and omits radiation. Ember adds
        # radiation once at runtime for either source treatment.
        pr = 0. if numerical_electrons else arad*T**4/3
        em, sm = E-3*pr/rho, S-4*pr/(rho*T)
        qt = -em/T
        qr = (P-pr)/(rho*T)
        qtt = em/T-(Et-12*pr/rho)/T
        qtr = (P*chit-P-3*pr)/(rho*T)
        qrr = (P*chir-P+pr)/(rho*T)
        # x=ln T, y=ln[rho/(T/1e6)^1.5]; phi=F_material/T.
        p=phi[k]
        p[0]=em/T-sm; p[1]=qr; p[2]=qrr
        p[3]=qt+1.5*qr; p[4]=qtr+1.5*qrr
        p[6]=qtt+3*qtr+2.25*qrr
        good.append(P>pr and chir>0 and Et>0 and cp>0 and ad>0 and delta>0)
    first=[1.,-8.,0.,8.,-1.]
    second=[-1.,16.,-30.,16.,-1.]
    def diff(k, col, stride, h, order):
        coeff=first if order==1 else second
        return sum(w*phi[k+(i-2)*stride][col] for i,w in enumerate(coeff))/(12*h**order)
    lines=['EMBER_HELMHOLTZ 1',
           (f'source "FreeEOS 3.0.0 EOS1; H=.7 He=.3; direct-source SHA256 {hashlib.sha256(original).hexdigest()}"' if a.hydrogen==.7 else
            f'source "FreeEOS 3.0.0 EOS1; H={a.hydrogen:g} He={helium:g}; direct-source SHA256 {hashlib.sha256(original).hexdigest()}"'),
           ('composition_proxy "Metals represented by helium; fixed X=.7, effective Y=.3; He3 and composition changes unsupported"' if a.hydrogen==.7 else
            f'composition_proxy "Metals represented by helium; fixed X={a.hydrogen:g}; composition family node"'),
           'composition '+' '.join(format(v,'.17g') for v in interior_composition(a.hydrogen)),
           f'log_t {nt-4} '+' '.join(format(v,'.17g') for v in lt[2:-2]),
           f'log_q {nq-4} '+' '.join(format(v,'.17g') for v in lq[2:-2]),'data']
    if baryonic:
        source_method = ('option 223; numerical electron integrals; radiation omitted'
                         if numerical_electrons else 'EOS1')
        lines[:4] = ['EMBER_HELMHOLTZ 2',
            f'source "FreeEOS 3.0.0 {source_method}; GS98; baryonic H={src["hydrogen"]:g} He3={src["helium3"]:g}; direct-source SHA256 {hashlib.sha256(original).hexdigest()}"',
            'composition_proxy '+json.dumps(src['approximation']),
            'basis baryon_mass',
            'metal_inventory gs98',
            'composition '+' '.join(format(v,'.17g') for v in src['composition'])]
    invalid=0
    for i in range(2,nt-2):
        for j in range(2,nq-2):
            k=i*nq+j
            p=phi[k][:]
            p[7]=diff(k,6,1,dy,1)
            p[5]=diff(k,2,nq,dx,1)
            p[8]=.5*(diff(k,6,1,dy,2)+diff(k,2,nq,dx,2))
            valid=all(good[k+d] and good[k+d*nq] for d in range(-2,3))
            invalid+=not valid
            lines.append(str(int(valid))+' '+' '.join(format(v,'.17g') for v in p))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text('\n'.join(lines)+'\n')
    print(f'{(nt-4)*(nq-4)} potential nodes; {invalid} excluded by source-stencil stability checks')


if __name__=='__main__':
    main()
