"""Validate explicit absent cells in a hot-only raw EOS density extension."""
import math


def inconsistent_source_rows(source):
    """Validate explicit exclusions without changing a returned source flag.

    Each excluded finite source row must actually fail the unchanged first-law
    criterion. Importers still reject any undeclared inconsistency.
    """
    excluded = set()
    for record in source.get('source_consistency_exclusions', []):
        k = record['index']
        if (type(k) is not int or not 0 <= k < len(source['data']) or k in excluded
                or record['criterion'] != 1e-7 or record['reason'] != 'first_law_defect'):
            raise ValueError('invalid source consistency exclusion')
        row = source['data'][k]
        if row is None or len(row) != 22 or row[0] != 0 or not all(map(math.isfinite, row)):
            raise ValueError('consistency exclusion requires a finite converged source row')
        rho, T, P = row[2:5]; chit, Er, Et, Sr, St = row[8:13]
        defect = max(abs(rho*Er/P+chit-1), abs(T*St/Et-1), abs(rho*T*Sr/P+chit))
        if defect <= 1e-7 or defect != record['maximum_absolute_defect']:
            raise ValueError('source consistency exclusion does not match its defect')
        excluded.add(k)
    return excluded


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
