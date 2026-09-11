"""Validate explicit absent cells in a hot-only raw EOS density extension."""
import math


def absent_source_rows(raw):
    ts, qs, data = raw['logT'], raw['logQ'], raw['data']
    if len(data) != len(ts)*len(qs):
        raise ValueError('source coverage dimensions differ')
    declaration = raw.get('source_coverage')
    if declaration is None:
        if any(row is None for row in data):
            raise ValueError('absent source row without a coverage declaration')
        return [False]*len(data)
    if (raw.get('composition_basis') != 'baryon_mass'
            or set(declaration) != {'kind', 'original_logQ_max', 'minimum_added_logT'}
            or declaration['kind'] != 'hot_density_extension'):
        raise ValueError('unknown source coverage declaration')
    boundary, minimum = declaration['original_logQ_max'], declaration['minimum_added_logT']
    if (not all(math.isfinite(v) for v in [boundary, minimum])
            or boundary not in qs[:-1] or minimum not in ts[1:]):
        raise ValueError('source coverage boundary must be an interior grid coordinate')
    missing = [t < minimum and q > boundary for t in ts for q in qs]
    if any(expected != (row is None) for expected, row in zip(missing, data, strict=True)):
        raise ValueError('source rows do not match declared absent coverage')
    return missing
