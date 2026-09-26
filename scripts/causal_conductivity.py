"""A constrained low-frequency conductivity continuation for diagnostics.

The positive, index-one absorption is proportional to the real conductivity.
Below a chosen join it is replaced by a Drude response continuous with the
measured spectrum; its collision rate satisfies a prescribed conductivity sum.
Above the join it retains the measured spectrum. The real dielectric response
comes from the same conductivity through a principal-value integral.

This is an explicit approximation to test, not a recovery of the source's
unmodified free-free calculation or an accepted stellar-opacity prescription.
Frequency u is dimensionless; D=c*rho/omega_scale, so epsilon_2=D*kappa/u.
"""
import math

import numpy as np
from scipy.optimize import brentq


def real_refractive_index(real, imaginary):
    magnitude = np.hypot(real, imaginary)
    # The alternative expression avoids cancellation for negative epsilon_1.
    denominator = 2*np.where(real < 0, magnitude-real, 1.)
    square = np.where(real >= 0, (magnitude+real)/2, imaginary**2/denominator)
    return np.sqrt(square)


class ContinuedConductivity:
    def __init__(self, u, opacity, D, plasma_u, join, refinement=1, tail_power=3):
        u, opacity = np.asarray(u, dtype=float), np.asarray(opacity, dtype=float)
        if (u.ndim != 1 or opacity.shape != u.shape or len(u) < 3
                or not np.isfinite(u).all() or not np.isfinite(opacity).all()
                or np.any(u <= 0) or np.any(opacity <= 0) or np.any(np.diff(u) <= 0)
                or not u[0] < join < u[-1] or D <= 0 or plasma_u <= 0
                or refinement not in (1, 2, 4) or tail_power not in (2, 3)):
            raise ValueError('invalid conductivity continuation inputs')
        # Refinement converges the analytic piecewise-linear transform toward
        # the log-linear spectral interpolation used by the other audits.
        kept = u[u > join]
        initial = np.r_[join, kept]
        grid = np.exp((np.log(initial[:-1])[:, None]
                       + np.arange(refinement)/refinement
                       * np.diff(np.log(initial))[:, None]).ravel())
        self.x = np.r_[grid, initial[-1]]
        self.x[0] = join
        self.k = np.exp(np.interp(np.log(self.x), np.log(u), np.log(opacity)))
        self.D, self.plasma_u, self.join = float(D), float(plasma_u), float(join)
        self.tail_power = tail_power
        self.m = np.diff(self.k)/np.diff(self.x)
        self.b = self.k[:-1]-self.m*self.x[:-1]
        self.high_area = float(np.sum(np.diff(self.x)*(self.k[:-1]+self.k[1:])/2))
        self.tail_area = float(self.k[-1]*self.x[-1]/(tail_power-1))
        self.target_area = math.pi*plasma_u**2/(2*D)
        self.required_low_area = self.target_area-self.high_area-self.tail_area
        ratio = self.required_low_area/(self.k[0]*join)
        if not ratio > 1:
            raise ValueError('no positive continuous Drude response can satisfy the requested sum')

        def scaled_area(log_y):
            y = math.exp(log_y)
            # For large y avoid loss of the positive difference from one.
            if y > 100:
                z = 1/y**2
                return 1+2*z/3-2*z*z/15+2*z**3/35
            return (y+1/y)*math.atan(1/y)

        log_y = brentq(lambda v: scaled_area(v)-ratio, -40, 40, xtol=1e-13)
        self.collision_u = join*math.exp(log_y)
        self.amplitude = float(self.k[0]*(join**2+self.collision_u**2))
        self.low_area = self.amplitude/self.collision_u*math.atan(join/self.collision_u)
        self.sum_relative_error = (self.low_area+self.high_area+self.tail_area)/self.target_area-1
        if abs(self.sum_relative_error) > 1e-11:
            raise ValueError('conductivity sum does not close')

    def absorption(self, u):
        u = np.asarray(u, dtype=float)
        if np.any(u <= 0) or not np.isfinite(u).all():
            raise ValueError('positive finite evaluation frequency required')
        value = np.interp(u, self.x, self.k)
        value = np.where(u < self.join, self.amplitude/(u*u+self.collision_u**2), value)
        return np.where(u > self.x[-1], self.k[-1]*(self.x[-1]/u)**self.tail_power, value)

    def _tail_integral(self, u):
        """PV integral of kappa_tail(x)/(x^2-u^2), without the subtraction."""
        b = self.x[-1]
        y = u/b
        if self.tail_power == 2:
            # I = kappa(b)/b * [atanh(y)/y-1]/y^2, with a real PV logarithm.
            with np.errstate(divide='ignore', invalid='ignore'):
                value = (.5*np.log(np.abs((1+y)/(1-y)))/y-1)/y**2
            series = sum(y**(2*j)/(2*j+3) for j in range(10))
        else:
            z = y*y
            with np.errstate(divide='ignore', invalid='ignore'):
                value = (-z-np.log(np.abs(1-z)))/(2*z*z)
            series = .5*sum(z**j/(j+2) for j in range(10))
        return self.k[-1]/b*np.where(y < .05, series, value)

    def dielectric(self, u, chunk=32):
        u = np.asarray(u, dtype=float)
        shape = u.shape
        flat = u.reshape(-1)
        opacity = self.absorption(flat)
        if np.any(flat == self.x[-1]):
            raise ValueError('evaluate away from the artificial upper tail join')
        result = np.empty(flat.size)
        for start in range(0, len(flat), chunk):
            v = flat[start:start+chunk]
            kv = opacity[start:start+chunk]
            # Subtract kappa(u). PV integral of 1/(x^2-u^2) over [0,infinity)
            # is zero, and the subtraction removes the internal singularity.
            kk = self.amplitude/(self.collision_u**2+v*v)
            with np.errstate(divide='ignore', invalid='ignore'):
                J = np.log(np.abs((self.join-v)/(self.join+v)))/(2*v)
                low = np.where(v <= self.join, 0., (kk-kv)*J)
            low -= kk*math.atan(self.join/self.collision_u)/self.collision_u
            left, right = self.x[:-1][None, :], self.x[1:][None, :]
            q = v[:, None]
            log_plus = np.log1p((right-left)/(left+q))
            inside = (q >= left) & (q <= right)
            with np.errstate(divide='ignore', invalid='ignore'):
                log_ratio = np.log(np.abs(right-q))-np.log(np.abs(left-q))-log_plus
                residual = self.m[None, :]*q+self.b[None, :]-kv[:, None]
                remainder = np.where(inside, 0., residual*log_ratio/(2*q))
            high = np.sum(self.m[None, :]*log_plus+remainder, axis=1)
            J_tail = np.log(np.abs((self.x[-1]+v)/(self.x[-1]-v)))/(2*v)
            tail = self._tail_integral(v)-kv*J_tail
            result[start:start+chunk] = 1+2*self.D/math.pi*(low+high+tail)
        imaginary = self.D*opacity/flat
        if not np.isfinite(result).all() or not np.isfinite(imaginary).all():
            raise ValueError('nonfinite causal dielectric response')
        return result.reshape(shape), imaginary.reshape(shape)

    def refractive_index(self, u):
        real, imaginary = self.dielectric(u)
        return real_refractive_index(real, imaginary)
