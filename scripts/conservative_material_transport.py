"""Implicit finite-volume species and energy transport on a fixed mass mesh.

The caller supplies thermodynamics and full material Onsager conductances,
including face area/distance and a common energy/entropy normalization.
Unknown primitives are independent mass fractions and log temperature.
Every face has one shared flux; exterior fluxes vanish. Conductances are held
fixed during a step. There is no burning, gravity work or structural update.
"""
from dataclasses import dataclass
import numpy as np


@dataclass
class ThermodynamicPoint:
    conserved: np.ndarray       # independent X and E/e0
    potential: np.ndarray       # (mu_exchange/T, -e0/T)/s0
    capacity: np.ndarray        # d(conserved)/d(potential), positive
    primitive_from_conserved: np.ndarray  # d(X,log T)/d(X,E/e0)
    entropy: float              # specific entropy / s0


def positive_solve(a, b):
    a = (a+a.T)/2
    scale = np.sqrt(np.diag(a))
    chol = np.linalg.cholesky(a/scale[:, None]/scale[None, :])
    return np.linalg.solve(chol.T, np.linalg.solve(chol, b/scale[:, None]))/scale[:, None]


def solve_chain(capacities, couplings, rhs):
    """Block elimination without subtracting two large diffusion diagonals.

The Schur contribution is the parallel sum of E and K. Whiten with E to
compute lambda/(1+lambda), preserving the small mass term for stiff K.
"""
    count, size = rhs.shape
    effective = np.array(capacities, copy=True)
    right = np.array(rhs, copy=True)
    transfers = []
    inverses = []
    for i in range(count-1):
        e, k = effective[i], couplings[i]
        chol = np.linalg.cholesky((e+e.T)/2)
        normalized = np.linalg.solve(chol, k)
        normalized = np.linalg.solve(chol, normalized.T).T
        value, vectors = np.linalg.eigh((normalized+normalized.T)/2)
        if value.min() < -1e-12*max(value.max(), 1.):
            raise ValueError('negative face conductance')
        # Tiny negative eigenvalues are not replaced: all supported nonzero
        # conductances must be positive definite, or exactly zero matrices.
        if np.any(k != 0) and value.min() <= 0:
            raise ValueError('face conductance is not resolved as positive')
        fraction = value/(1+value)
        parallel = (vectors*fraction)@vectors.T
        transfer = np.linalg.solve(chol.T, (chol@parallel).T).T
        effective[i+1] += chol@parallel@chol.T
        right[i+1] += transfer@right[i]
        transfers.append(transfer.T)
        inv = positive_solve(e+k, np.eye(size))
        inverses.append(inv)
    answer = np.empty_like(rhs)
    answer[-1] = positive_solve(effective[-1], right[-1, :, None])[:, 0]
    for i in range(count-2, -1, -1):
        # (E+K)^-1 K = [K(E+K)^-1]^T.
        answer[i] = inverses[i]@right[i]+transfers[i]@answer[i+1]
    return answer


def implicit_step(*, old_primitives, initial_guess, mass, conductance, dt,
                  thermodynamics, tolerance=2e-11, max_iterations=60):
    """Return a converged backward-Euler update or raise without modifying input.

An interior initial guess is a Newton iterate, never a floor on the solution.
The original, possibly zero, abundances remain in the conserved right side.
The caller removes a species absent from the entire connected domain.
"""
    old = np.asarray(old_primitives, dtype=float)
    q = np.array(initial_guess, dtype=float, copy=True)
    mass, k = np.asarray(mass, dtype=float), np.asarray(conductance, dtype=float)
    n, size = old.shape
    if (q.shape != old.shape or mass.shape != (n,) or k.shape != (n-1, size, size) or
            not np.isfinite(mass).all() or np.any(mass <= 0) or
            not np.isfinite(k).all() or not np.isfinite(dt) or dt <= 0 or
            not np.array_equal(k, np.swapaxes(k, -1, -2))):
        raise ValueError('positive mass, time and symmetric finite conductances required')
    for matrix in k:
        if np.any(matrix != 0):
            np.linalg.cholesky(matrix)
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
        flow = -np.einsum('nij,nj->ni', k, np.diff(w, axis=0))
        residual = mass[:, None]*(u-old_u)
        residual[:-1] += dt*flow
        residual[1:] -= dt*flow
        norm = float(np.max(abs(residual)/mass[:, None]/scale))
        return points, u, w, flow, residual, norm
    history = []
    for iteration in range(max_iterations):
        points, u, w, flow, residual, norm = evaluate(q)
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
            if trial[-1] < norm*(1-1e-4*damping) or trial[-1] <= tolerance:
                q = candidate
                accepted = True
                break
            damping *= .5
        if not accepted:
            raise RuntimeError(f'implicit transport line search failed at residual {norm}')
    raise RuntimeError(f'implicit transport did not converge: {history[-1]}')
