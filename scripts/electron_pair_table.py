"""Positive interpolation of the nonrelativistic Born/Pauli pair matrices.

The momentum mode remains exactly null. Higher energy modes use the same
continuous, positive-leading-coefficient polynomial convention as the source
integrals. Positivity is structural; accuracy needs independent source checks.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RectBivariateSpline


class ElectronPairTable:
    def __init__(self, report, *, interpolation_degree=3):
        path = Path(report)
        data = json.loads(path.read_text())
        if data['outcome'] != 'completed_unvalidated_table' or data['failed_sources']:
            raise ValueError('complete source grid required')
        if interpolation_degree not in (1, 3) or data['mode_count'] != 10:
            raise ValueError('linear or cubic interpolation of ten modes required')
        self.report_sha256 = hashlib.sha256(path.read_bytes()).hexdigest()
        self.eta = np.asarray(data['eta_grid'], dtype=float)
        self.b = np.asarray(data['b_thermal_grid'], dtype=float)
        for grid in (self.eta, self.b):
            if (grid.ndim != 1 or len(grid) <= interpolation_degree or
                    not np.isfinite(grid).all() or np.any(np.diff(grid) <= 0)):
                raise ValueError('finite increasing grid required')
        if np.any(self.b <= 0):
            raise ValueError('positive screening grid required')
        self.log_b = np.log(self.b)
        self.lower = np.tril_indices(9)
        self.diagonal = self.lower[0] == self.lower[1]
        values = np.empty((len(self.eta), len(self.b), len(self.lower[0])))
        seen = set()
        self.sources = []
        for entry in data['completed_sources']:
            source_path = Path(entry['path'])
            raw = source_path.read_bytes()
            if hashlib.sha256(raw).hexdigest() != entry['sha256']:
                raise ValueError('source identity changed: ' + str(source_path))
            source = json.loads(raw)
            i = np.flatnonzero(self.eta == source['eta'])
            j = np.flatnonzero(self.b == source['b_thermal'])
            if len(i) != 1 or len(j) != 1 or (int(i[0]), int(j[0])) in seen:
                raise ValueError('source coordinates absent or repeated')
            i, j = int(i[0]), int(j[0])
            seen.add((i, j))
            if (source['degree'] != 9 or source['statistics'] != 'fermi' or
                    source['orders'] != [80, 40, 40, 48, 48] or
                    source['tail'] != 50 or source['normalization_tail'] != 80):
                raise ValueError('source integral conventions differ')
            matrix = np.asarray(source['collision_matrix'], dtype=float)
            if (matrix.shape != (10, 10) or not np.isfinite(matrix).all() or
                    np.any(matrix[0] != 0) or np.any(matrix[:, 0] != 0) or
                    not np.array_equal(matrix, matrix.T)):
                raise ValueError('symmetric source with exact momentum null mode required')
            chol = np.linalg.cholesky(matrix[1:, 1:])
            encoded = (chol / np.diag(chol)[:, None])[self.lower]
            encoded[self.diagonal] = np.log(np.diag(chol))
            values[i, j] = encoded
            self.sources.append(source)
        if len(seen) != len(self.eta) * len(self.b):
            raise ValueError('source grid is incomplete')
        self.splines = [RectBivariateSpline(self.eta, self.log_b, values[:, :, k],
                         kx=interpolation_degree, ky=interpolation_degree, s=0)
                        for k in range(values.shape[-1])]
        self.interpolation_degree = interpolation_degree

    def matrices(self, eta, b_thermal):
        eta, b = np.broadcast_arrays(np.asarray(eta, dtype=float),
                                     np.asarray(b_thermal, dtype=float))
        if (not np.isfinite(eta).all() or not np.isfinite(b).all() or
                np.any(eta < self.eta[0]) or np.any(eta > self.eta[-1]) or
                np.any(b < self.b[0]) or np.any(b > self.b[-1])):
            raise ValueError('electron pair query outside the supplied grid')
        encoded = np.stack([s.ev(eta.ravel(), np.log(b.ravel()))
                            for s in self.splines], axis=-1)
        chol = np.zeros((eta.size, 9, 9))
        chol[:, self.lower[0], self.lower[1]] = encoded
        diagonal = np.exp(encoded[:, self.diagonal])
        chol *= diagonal[:, :, None]
        chol[:, np.arange(9), np.arange(9)] = diagonal
        result = np.zeros((eta.size, 10, 10))
        result[:, 1:, 1:] = chol @ np.swapaxes(chol, -1, -2)
        if not np.isfinite(result).all():
            raise ValueError('nonfinite interpolated matrix')
        return result.reshape(eta.shape + (10, 10))
