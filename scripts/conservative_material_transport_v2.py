"""Conservative material transport with state-dependent face coefficients."""
import numpy as np
from conservative_material_transport import ThermodynamicPoint, positive_solve, solve_chain


def implicit_step(*, old_primitives, initial_guess, mass, conductance, dt,
                  thermodynamics, tolerance=2e-11, max_iterations=60):
    """Return a converged backward-Euler update or raise without modifying input.

Conductance may be a fixed array or a function of the trial primitives.
    A symmetric block solve freezes its derivative during each Newton direction;
    the line search evaluates the complete changing-coefficient residual.
    An interior initial guess is a Newton iterate, never a floor on the solution.
The original, possibly zero, abundances remain in the conserved right side.
The caller removes a species absent from the entire connected domain.
"""
    old = np.asarray(old_primitives, dtype=float)
    q = np.array(initial_guess, dtype=float, copy=True)
    mass = np.asarray(mass, dtype=float)
    n, size = old.shape
    if (q.shape != old.shape or mass.shape != (n,) or
            not np.isfinite(mass).all() or np.any(mass <= 0) or
            not np.isfinite(dt) or dt <= 0):
        raise ValueError('positive mass, time and symmetric finite conductances required')
    def matrix_for(candidate):
        k = np.asarray(conductance(candidate) if callable(conductance) else conductance, dtype=float)
        if (k.shape != (n-1, size, size) or not np.isfinite(k).all() or
                not np.array_equal(k, np.swapaxes(k, -1, -2))):
            raise ValueError('finite symmetric face conductances required')
        for matrix in k:
            if np.any(matrix != 0):
                np.linalg.cholesky(matrix)
        return k
    old_points = [thermodynamics(i, v, derivatives=False) for i, v in enumerate(old)]
    old_u = np.array([p.conserved for p in old_points])
    old_entropy = float(mass@np.array([p.entropy for p in old_points]))
    scale = np.maximum(1., np.max(abs(old_u), axis=0))
    def evaluate(candidate):
        if (not np.isfinite(candidate).all() or np.any(candidate[:, :-1] <= 0) or
                np.any(candidate[:, :-1].sum(axis=1) >= 1)):
            raise ValueError('Newton iterate leaves the composition interior')
        points = [thermodynamics(i, v, derivatives=True) for i, v in enumerate(candidate)]
        u = np.array([p.conserved for p in points])
        w = np.array([p.potential for p in points])
        k = matrix_for(candidate)
        flow = -np.einsum('nij,nj->ni', k, np.diff(w, axis=0))
        residual = mass[:, None]*(u-old_u)
        residual[:-1] += dt*flow
        residual[1:] -= dt*flow
        norm = float(np.max(abs(residual)/mass[:, None]/scale))
        return points, u, w, flow, residual, norm, k
    history = []
    for iteration in range(max_iterations):
        points, u, w, flow, residual, norm, k = evaluate(q)
        history.append(norm)
        if norm <= tolerance:
            entropy = float(mass@np.array([p.entropy for p in points]))
            dissipation = float(-np.sum(flow*np.diff(w, axis=0)))
            return dict(primitives=q, conserved=u, face_flux=flow, iterations=iteration,
                residual=norm, residual_history=history,
                integrated_conservation_error=mass@(u-old_u),
                entropy_change=entropy-old_entropy,
                backward_euler_entropy_bound=dt*dissipation)
        capacity = np.array([p.capacity for p in points])
        delta_w = solve_chain(mass[:, None, None]*capacity, dt*k, -residual)
        delta_u = np.einsum('nij,nj->ni', capacity, delta_w)
        delta_q = np.array([p.primitive_from_conserved@du for p, du in zip(points, delta_u, strict=True)])
        damping = 1.
        accepted = False
        for _ in range(50):
            candidate = q+damping*delta_q
            try:
                trial = evaluate(candidate)
            except (ValueError, FloatingPointError, OverflowError):
                damping *= .5
                continue
            if trial[5] < norm*(1-1e-4*damping) or trial[5] <= tolerance:
                q = candidate
                accepted = True
                break
            damping *= .5
        if not accepted:
            raise RuntimeError(f'implicit transport line search failed at residual {norm}')
    raise RuntimeError(f'implicit transport did not converge: {history[-1]}')
