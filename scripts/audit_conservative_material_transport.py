#!/usr/bin/env python3
"""Independent thermodynamic and conservative controls for material transport."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq
from scipy.special import xlogy

from conservative_material_transport import ThermodynamicPoint, implicit_step, positive_solve, solve_chain


class AnalyticMixture:
    """Classical ideal ions plus an optional cold degenerate electron energy.

    Units set Rgas=1. The analytic potential, not an opacity or stellar EOS,
    supplies an independent nonlinear control with composition-energy coupling.
    """
    def __init__(self, masses, charges, binding, density, degeneracy=0., energy_scale=1.):
        self.a = np.asarray(masses, dtype=float)
        self.z = np.asarray(charges, dtype=float)
        self.binding = np.asarray(binding, dtype=float)
        self.rho = np.asarray(density, dtype=float)
        self.degeneracy = degeneracy
        self.energy_scale = energy_scale

    def composition(self, q):
        return np.r_[q[:-1], 1-sum(q[:-1])]

    def energy_terms(self, i, x):
        cv = float(x@(1.5/self.a))
        ye = float(x@(self.z/self.a))
        factor = self.degeneracy*self.rho[i]**(2/3)
        cold = float(x@self.binding+factor*ye**(5/3))
        return cv, ye, factor, cold

    def primitive(self, i, u):
        x = np.r_[u[:-1], 1-sum(u[:-1])]
        cv, ye, factor, cold = self.energy_terms(i, x)
        t = (u[-1]*self.energy_scale-cold)/cv
        return np.r_[u[:-1], np.log(t)]

    def __call__(self, i, q, derivatives=True):
        x = self.composition(q)
        if np.any(x < 0):
            raise ValueError('negative analytic abundance')
        t = float(np.exp(q[-1]))
        if not np.isfinite(t) or t <= 0:
            raise ValueError('invalid analytic temperature')
        cv, ye, factor, cold = self.energy_terms(i, x)
        energy = cv*t+cold
        entropy = cv*np.log(t)-float(np.sum(xlogy(x/self.a, self.rho[i]*x/self.a)))
        u = np.r_[x[:-1], energy/self.energy_scale]
        if not derivatives:
            return ThermodynamicPoint(u, None, None, None, entropy)
        if np.any(x <= 0):
            raise ValueError('entropy derivatives require an interior state')
        dy = self.z[:-1]/self.a[:-1]-self.z[-1]/self.a[-1]
        dcv = 1.5/self.a[:-1]-1.5/self.a[-1]
        dcold = self.binding[:-1]-self.binding[-1]+5/3*factor*ye**(2/3)*dy
        ec = dcv*t+dcold
        g = dcv*(1-np.log(t))+(np.log(self.rho[i]*x[:-1]/self.a[:-1])+1)/self.a[:-1]
        g -= (np.log(self.rho[i]*x[-1]/self.a[-1])+1)/self.a[-1]
        g += dcold/t
        hphi = np.diag(1/(self.a[:-1]*x[:-1]))+1/(self.a[-1]*x[-1])
        hphi += 10/9*factor*ye**(-1/3)*np.outer(dy, dy)/t
        size = len(x)
        h = np.empty((size, size))
        h[:-1, :-1] = hphi+np.outer(ec, ec)/(cv*t*t)
        h[:-1, -1] = -self.energy_scale*ec/(cv*t*t)
        h[-1, :-1] = h[:-1, -1]
        h[-1, -1] = self.energy_scale**2/(cv*t*t)
        capacity = positive_solve(h, np.eye(size))
        conversion = np.eye(size)
        conversion[-1, :-1] = -ec/(cv*t)
        conversion[-1, -1] = self.energy_scale/(cv*t)
        return ThermodynamicPoint(u, np.r_[g, -self.energy_scale/t], capacity, conversion, entropy)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    assert not args.output.exists()
    start = time.process_time()
    checks, runs = [], []
    def check(name, error, tolerance=2e-10, **meta):
        checks.append(dict(name=name, error=float(error), tolerance=tolerance,
            passed=bool(np.isfinite(error) and abs(error) <= tolerance), **meta))
    rng = np.random.default_rng(2046)
    # Independently assemble the entire matrix for noncommuting block controls.
    for stiffness in (1e-7, 1., 1e5):
        n, size = 7, 3
        a = rng.normal(size=(n, size, size))
        capacity = a@np.swapaxes(a, -1, -2)+np.eye(size)
        b = rng.normal(size=(n-1, size, size))
        k = stiffness*(b@np.swapaxes(b, -1, -2)+np.eye(size))
        rhs = rng.normal(size=(n, size))
        full = np.zeros((n*size, n*size))
        for i in range(n):
            section = slice(i*size, (i+1)*size)
            full[section, section] += capacity[i]
            if i+1 < n:
                other = slice((i+1)*size, (i+2)*size)
                full[section, section] += k[i]
                full[other, other] += k[i]
                full[section, other] -= k[i]
                full[other, section] -= k[i]
        answer = solve_chain(capacity, k, rhs)
        reference = np.linalg.solve(full, rhs.ravel()).reshape(n, size)
        check('independent dense block solution', np.max(abs(answer-reference))/np.max(abs(reference)),
              2e-9, stiffness=stiffness)
    for stiffness in (0., 1., 1e6, 1e12):
        vector = np.array([.3, -.8, 1.1])
        rhs = np.einsum('nij,j->ni', capacity, vector)
        answer = solve_chain(capacity, stiffness*np.broadcast_to(np.eye(3), (6, 3, 3)), rhs)
        check('stiff uniform mode', np.max(abs(answer-vector)), 2e-10, stiffness=stiffness)
    for degenerate in (0., 3.):
        eos = AnalyticMixture([1, 3, 4], [1, 2, 2], [1., -.4, .2], [1., 4.], degenerate)
        for i in (0, 1):
            q = np.array([.2, .07, np.log(1.3)])
            p = eos(i, q)
            h = positive_solve(p.capacity, np.eye(3))
            for j in range(3):
                du = np.zeros(3);du[j] = 2e-6
                upper = eos(i, eos.primitive(i, p.conserved+du))
                lower = eos(i, eos.primitive(i, p.conserved-du))
                derivative = (upper.potential-lower.potential)/(2*du[j])
                check('entropy Hessian from finite differences', np.max(abs(derivative-h[:, j]))/np.max(abs(h[:, j])),
                      2e-8, degeneracy=degenerate, density=eos.rho[i], column=j)
    # A single constituent gives an independent scalar two-cell heat solution.
    eos = AnalyticMixture([4], [2], [0.], [1., 1.])
    old = np.log([[2.], [1.]])
    mass = np.array([.3, .7]);dt = .4;conductance = np.array([[[.2]]])
    answer = implicit_step(old_primitives=old, initial_guess=old, mass=mass,
        conductance=conductance, dt=dt, thermodynamics=eos)
    cv = 1.5/4;mean = float(mass@np.exp(old[:, 0]))
    def heat_equation(t0):
        t1 = (mean-mass[0]*t0)/mass[1]
        return mass[0]*cv*(t0-2)+dt*.2*(1/t1-1/t0)
    t0 = brentq(heat_equation, mean, 2., xtol=1e-14)
    check('independent two-cell heat root', abs(np.exp(answer['primitives'][0, -1])/t0-1))
    for degenerate in (0., 2.):
        n = 12
        mass = np.linspace(1., 2., n);mass /= sum(mass)
        eos = AnalyticMixture([1, 3, 4], [1, 2, 2], [1., -.4, .2], np.linspace(.7, 1.3, n), degenerate)
        x = np.empty((n, 3));x[:n//2] = [.65, 0., .35];x[n//2:] = [0., .15, .85]
        old = np.c_[x[:, :2], np.log(np.linspace(.8, 1.4, n))]
        matrix = np.array([[1., .2, .1], [.2, .6, -.15], [.1, -.15, .9]])
        k = np.broadcast_to(matrix, (n-1, 3, 3)).copy()
        for dt in (1e-4, .05, 100.):
            answers = []
            for initial_mixing in (.01, .1):
                interior = (1-initial_mixing)*x+initial_mixing*(mass@x)
                guess = np.c_[interior[:, :2], old[:, -1]]
                result = implicit_step(old_primitives=old, initial_guess=guess, mass=mass,
                    conductance=k, dt=dt, thermodynamics=eos)
                answers.append(result)
                meta = dict(degeneracy=degenerate, dt=dt, initial_mixing=initial_mixing)
                check('conserved species and full energy', np.max(abs(result['integrated_conservation_error'])), **meta)
                check('backward-Euler entropy inequality', min(result['entropy_change']-result['backward_euler_entropy_bound'], 0.), **meta)
                final_x = np.c_[result['primitives'][:, :-1], 1-result['primitives'][:, :-1].sum(axis=1)]
                check('positive incoming abundance without clipping', min(float(final_x.min()), 0.), **meta)
                assert np.all(final_x[x == 0] > 0)
                check('cell mass fraction sum', np.max(abs(final_x.sum(axis=1)-1)), **meta)
                runs.append(dict(**meta, iterations=result['iterations'], residual=result['residual'],
                    minimum_final_fraction=float(final_x.min()), entropy_change=result['entropy_change'],
                    backward_euler_entropy_bound=result['backward_euler_entropy_bound'],
                    final_primitives=result['primitives'].tolist()))
            check('independence of interior Newton guess', np.max(abs(answers[0]['conserved']-answers[1]['conserved'])),
                  1e-9, degeneracy=degenerate, dt=dt)
    passed = all(c['passed'] for c in checks)
    root = Path(__file__).resolve().parents[1]
    files = [Path(__file__), root/'scripts/conservative_material_transport.py']
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(), outcome='passed' if passed else 'failed',
        accepted_for_stellar_evolution=False, checks=checks, runs=runs, cpu_seconds=time.process_time()-start,
        input_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        limitations=['Analytic thermodynamics and prescribed constant Onsager conductances; no stellar opacity or EOS accepted.',
                     'Face coefficients are frozen during the implicit step; changing coefficients require outer iteration or timestep control.',
                     'No burning, convective mixing, hydrostatic work or stellar structure update.'])
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False);stream.write('\n')
    print(json.dumps(dict(outcome=report['outcome'], checks=len(checks),
        failed=[c for c in checks if not c['passed']], cpu_seconds=report['cpu_seconds']), indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
