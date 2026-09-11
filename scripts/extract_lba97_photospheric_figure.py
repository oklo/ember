#!/usr/bin/env python3
"""Recover Figure 6 geometry and an approximate opacity scale from LBA97.

The figure has no colorbar. Fit its gray levels to AJR83 Table 2 at cell
centers, using 3.5 <= log T < 3.7; withhold cooler and hotter cells. This
calibrates a diagram, not a replacement opacity table for stellar evolution.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pymupdf
from scipy.interpolate import RegularGridInterpolator

# AJR83, ApJ 272, 773, Table 2 (p. 776), read from the published page.
# Columns log rho = -8, -7, -6, -5, -4, -3, -2. Missing cells stay missing.
TEMPERATURE = [3.30, 3.35, 3.40, 3.45, 3.50, 3.55, 3.60, 3.65, 3.70, 3.80, 3.90, 4.00]
OPACITY = [
    [-2.093, -1.919, -1.858, -1.828, None, None, None],
    [-1.821, -1.764, -1.734, -1.709, None, None, None],
    [-1.729, -1.619, -1.573, -1.509, None, None, None],
    [-1.957, -1.495, -1.343, -1.227, None, None, None],
    [-2.410, -1.530, -1.122, -.908, -.585, None, None],
    [-2.389, -1.657, -1.025, -.613, -.290, None, None],
    [-2.035, -1.415, -.838, -.343, .065, None, None],
    [-1.854, -1.114, -.492, -.002, .459, .913, None],
    [-1.769, -.988, -.243, .341, .817, 1.261, None],
    [-.832, -.362, .181, .782, 1.339, 1.780, None],
    [.426, .746, 1.126, 1.574, 2.101, 2.510, None],
    [1.799, 2.026, 2.270, 2.578, 3.008, 3.537, 4.016],
]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['lba97_pdf', 'ajr83_pdf', 'output']:
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    if digest(args.lba97_pdf) != '3f5298af22f71bdb50936a35c807b5f6081a11e97a89780e77421793351f8e4e':
        raise ValueError('unexpected LBA97 source PDF')
    document = pymupdf.open(args.lba97_pdf)
    page = document[8]
    paths = page.get_drawings()
    # Axis corners in the original PDF, independently fixed by labeled ticks.
    left, right, top, bottom = 125.765, 518.424, 406.191, 683.362
    xcoord = lambda x: -8 + 6*(x-left)/(right-left)
    ycoord = lambda y: 3 + (bottom-y)/(bottom-top)
    cells, circles = [], []
    for index, path in enumerate(paths):
        rect = path['rect']
        if path['type'] == 'f':
            if len(path['items']) != 1 or path['items'][0][0] != 're':
                raise ValueError('unexpected filled geometry')
            cells.append({'logrho': [xcoord(rect.x0), xcoord(rect.x1)],
                          'logT': [ycoord(rect.y1), ycoord(rect.y0)],
                          'gray': float(np.mean(path['fill']))})
        elif len(path['items']) in [4, 10, 15]:
            if path['items'][0][0] == 'l':
                points = np.array([list(item[1]) for item in path['items']])
                x, y, c = np.linalg.lstsq(np.c_[2*points, np.ones(len(points))],
                                        np.sum(points**2, axis=1), rcond=None)[0]
                radius = np.sqrt(c+x*x+y*y)
                if max(abs(np.linalg.norm(points-[x,y], axis=1)-radius)) > .002:
                    raise ValueError('marker is not a circular polygon')
                diameter = 2*radius
            elif all(item[0] == 'c' for item in path['items']):
                x, y = (rect.x0+rect.x1)/2, (rect.y0+rect.y1)/2
                diameter = rect.width
            else:
                raise ValueError('unexpected marker geometry')
            circles.append([index, xcoord(x), ycoord(y), diameter])
    if len(cells) != 156 or len(circles) != 4662:
        raise ValueError('unexpected Figure 6 geometry')

    centers = np.array([[np.mean(c['logT']), np.mean(c['logrho'])] for c in cells])
    # Suppress only PDF coordinate rounding, well below its plotted precision.
    estimate = RegularGridInterpolator((TEMPERATURE, np.arange(-8, -1)),
        np.array(OPACITY, dtype=float), bounds_error=False, fill_value=np.nan)(centers.round(4))
    gray = np.array([c['gray'] for c in cells])
    available = np.isfinite(estimate) & (gray > .51) & (gray < .94)
    fit = available & (centers[:,0] >= 3.5) & (centers[:,0] < 3.7)
    withheld = available & ~fit
    slope, intercept = np.polyfit(estimate[fit], gray[fit], 1)
    inferred = (gray-intercept)/slope
    residual = inferred[withheld]-estimate[withheld]
    if slope >= 0 or sum(fit) != 15 or sum(withheld) != 19 or max(abs(residual)) > .15:
        raise ValueError('opacity calibration failed its withheld-cell check')
    for cell, value in zip(cells, inferred):
        cell['estimated_log10_kappa_cm2_g'] = float(value)
    comparisons = []
    for i in np.flatnonzero(available):
        comparisons.append({'cell': int(i), 'coordinates_logT_logrho': centers[i].tolist(),
                            'AJR83_interpolated_log_kappa': float(estimate[i]),
                            'figure_estimated_log_kappa': float(inferred[i]),
                            'used_in_fit': bool(fit[i])})

    markers = np.array(circles)
    breaks = [0] + [i+1 for i, (a,b) in enumerate(zip(markers, markers[1:]))
                   if np.hypot(a[1]-b[1], 10*(a[2]-b[2])) > .5] + [len(markers)]
    if len(breaks) != 7:
        raise ValueError('expected six published mass sequences')
    # Published order and temperature maxima identify 0.06,.08,.10,.12,.16,.25 Msun.
    track = markers[breaks[2]:breaks[3]]
    hottest = np.argmax(track[:,2])
    zams = np.argmin(track[:hottest,2])
    if abs(track[zams,2]-np.log10(2228)) > .001:
        raise ValueError('Figure 6 minimum disagrees with the published ZAMS temperature')
    radius_zams = 10**(-3.38/2)*(5772/2228)**2
    radii = track[:,3]/track[zams,3]*radius_zams
    result = {'scope': __doc__, 'source': 'Laughlin, Bodenheimer & Adams 1997, Figure 6',
              'source_doi': '10.1086/304125', 'page': 428,
              'calibration_source': 'Alexander, Johnson & Rypma 1983, Table 2, p. 776',
              'input_sha256': {str(p): digest(p) for p in [args.lba97_pdf, args.ajr83_pdf, Path(__file__)]},
              'calibration': {'gray_equals': 'intercept + slope * log10(kappa / (cm^2/g))',
                  'slope': float(slope), 'intercept': float(intercept), 'fit_cells': int(sum(fit)),
                  'withheld_cells': int(sum(withheld)),
                  'withheld_rms_log10_kappa': float(np.sqrt(np.mean(residual**2))),
                  'withheld_max_abs_log10_kappa': float(max(abs(residual))),
                  'limitations': 'Approximate diagram calibration. Cell-center assignment, coarse shading and source interpolation contribute uncertainty. The checks do not bound grain opacity errors or recover the original numerical table. Extremal gray levels are outside the fit.',
                  'comparisons': comparisons},
              'cells': cells,
              'track': {'mass_Msun': .1, 'selection_path_range': [int(track[0,0]), int(track[-1,0])],
                  'radius_normalization': {'published_ZAMS_T_K': 2228, 'published_ZAMS_log10_L_Lsun': -3.38,
                     'adopted_solar_T_K': 5772, 'inferred_ZAMS_R_Rsun': radius_zams,
                     'ZAMS_marker_index': int(zams), 'note': 'Relative diameters are recovered from vector markers; absolute radii are anchored to the published ZAMS luminosity and temperature.'},
                  'log10_density_g_cm3': track[:,1].tolist(), 'log10_temperature_K': track[:,2].tolist(),
                  'R_Rsun': radii.tolist(), 'original_diameter_pdf_points': track[:,3].tolist()}}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result['calibration'] | {'comparisons': len(comparisons)}))


if __name__ == '__main__':
    main()
