#!/usr/bin/env python3
"""Independent BGS reference: eigenvalues of the linearized fluid equations.

The native implementation solves the published cubic and its stationarity
condition. This reference instead diagonalizes the velocity/temperature/
composition matrix and maximizes its largest real eigenvalue over log(l^2).
No MESA source is used. Equations: Brown et al. 2013, arXiv:1212.1688.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
from decimal import Decimal, localcontext
import numpy as np
from scipy.optimize import minimize_scalar


def reference(pr, tau, ratio):
    critical_q = math.sqrt((1 / ratio - tau) / tau)

    def growth(logq):
        q = math.exp(logq)
        matrix = np.array([[-pr*q, pr, -pr], [-1., -q, 0.],
                           [-1./ratio, 0., -tau*q]])
        roots = np.linalg.eigvals(matrix)
        value = roots[np.argmax(roots.real)]
        assert abs(value.imag) < 1e-10 * max(1e-30, abs(value.real))
        return value.real

    # Search the entire unstable q interval on a log mesh, then refine its
    # maximum. This avoids assuming the native stationarity root is maximal.
    upper = math.log(critical_q) - 1e-7
    grid = np.linspace(min(-25., upper-25.), upper, 500)
    values = [growth(v) for v in grid]
    i = int(np.argmax(values))
    assert 0 < i < len(grid)-1
    result = minimize_scalar(lambda x: -growth(x), bounds=(grid[i-1], grid[i+1]),
                             method='bounded', options={'xatol': 2e-10})
    assert result.success
    # At small Pr the growth maximum is very flat: double-precision matrix
    # eigenvalues locate lambda accurately but can misplace q enough to alter
    # the predicted flux. Refine the same matrix determinant with decimal
    # arithmetic and derivative-free maximization; no native cubic or q
    # stationarity expression is used here.
    with localcontext() as context:
        context.prec = 60
        p, t, r = [Decimal(str(v)) for v in (pr, tau, ratio)]

        def determinant(lam, q):
            matrix = [[lam+p*q, -p, p], [Decimal(1), lam+q, Decimal(0)],
                      [1/r, Decimal(0), lam+t*q]]
            a, b, c = matrix[0]; d, e, f = matrix[1]; g, h, j = matrix[2]
            return a*(e*j-f*h)-b*(d*j-f*g)+c*(d*h-e*g)

        def decimal_growth(q):
            lam = Decimal(str(growth(math.log(float(q)))))
            assert lam > 0
            for _ in range(30):
                h = lam*Decimal('1e-22')
                derivative = (determinant(lam+h,q)-determinant(lam-h,q))/(2*h)
                step = determinant(lam,q)/derivative
                lam -= step
                if abs(step) < abs(lam)*Decimal('1e-45'):
                    return lam
            raise AssertionError('decimal eigenvalue refinement failed')

        left, right = [Decimal(str(math.exp(v))) for v in (grid[i-1], grid[i+1])]
        fraction = (Decimal(5).sqrt()-1)/2
        a = right-fraction*(right-left); b = left+fraction*(right-left)
        fa, fb = decimal_growth(a), decimal_growth(b)
        for _ in range(90):
            if fa > fb:
                right, b, fb = b, a, fa
                a = right-fraction*(right-left); fa = decimal_growth(a)
            else:
                left, a, fa = a, b, fb
                b = left+fraction*(right-left); fb = decimal_growth(b)
        q_decimal = (left+right)/2
        q, lam = float(q_decimal), float(decimal_growth(q_decimal))
    assert lam > 0
    return [pr, tau, ratio, lam, q,
            49*lam*lam/(q*(lam+q)), 49*lam*lam/(tau*q*(lam+tau*q))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('output', type=Path)
    args = ap.parse_args()
    cases = [(4e-6, 2e-6, 1700.)]
    for pr in [1e-7, 1e-6, 1e-4, .01, .1]:
        for factor in [.1, .5, 1., 2.]:
            tau = pr*factor
            for ratio in [1.001, 2., 10., .5/tau, .99/tau]:
                if 1 < ratio < 1/tau:
                    cases.append((pr, tau, ratio))
    rows = [reference(*p) for p in cases]
    assert abs((1+rows[0][-1])/1294-1) < .005, 'published section 6.2 example'
    with args.output.open('x') as f:
        f.write('# Pr tau R0 lambda l_squared Nu_T_minus_1 Nu_mu_minus_1\n')
        for row in rows:
            f.write(' '.join(format(v, '.17g') for v in row)+'\n')
    receipt = dict(cases=len(rows), method=__doc__, source=__file__,
                   source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   output_sha256=hashlib.sha256(args.output.read_bytes()).hexdigest(),
                   published_example_nu_mu=1+rows[0][-1])
    args.output.with_suffix('.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
