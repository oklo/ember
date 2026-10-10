"""Build the compact SQLite view of the published CSV tables (standard library)."""
from pathlib import Path
import csv
import json
import math
import os
import re
import sqlite3
import tempfile

ROOT = Path(__file__).resolve().parent
DESCRIPTIONS = {
    'clusters': 'Candidate mass inference; uncorrected survey selection.',
    'field': 'Field population and model sensitivities; see field report.',
    'dynamics': 'Dynamical histories and observed aperture diagnostics.',
    'ensemble': 'Complete 17-model comparison; birth ratios are assumed.',
    'stellar-fit': 'Stellar fits and brown-dwarf predictions; incomplete selection.',
    'membership': 'Candidate retention and conditional membership recovery.',
}


def quote(name):
    return '"' + name.replace('"', '""') + '"'


def table_name(path):
    return re.sub(r'[^a-zA-Z0-9_]', '_', path.parent.name + '__' + path.stem)


def read_table(path):
    with path.open(newline='') as f:
        reader = csv.reader(f)
        columns = next(reader)
        rows = list(reader)
    if len(columns) != len(set(columns)) or any(len(r) != len(columns) for r in rows):
        raise ValueError(f'Invalid CSV shape: {path}')
    types = []
    for i, column in enumerate(columns):
        values = [r[i] for r in rows if r[i] not in ('', 'nan', 'NaN')]
        numeric = bool(values) and not re.search(r'(^|_)(id|source_id|gaia_id)$', column.lower())
        try:
            numeric = numeric and all(math.isfinite(float(v)) for v in values)
        except ValueError:
            numeric = False
        types.append('REAL' if numeric else 'TEXT')
    converted = [tuple(None if v in ('', 'nan', 'NaN') else float(v) if kind == 'REAL' else v
                       for v, kind in zip(row, types)) for row in rows]
    return columns, types, converted


def build():
    fd, filename = tempfile.mkstemp(prefix='.results-', suffix='.sqlite', dir=ROOT)
    os.close(fd)
    try:
        with sqlite3.connect(filename) as db:
            db.execute('CREATE TABLE datasets (table_name TEXT PRIMARY KEY, csv_path TEXT, '
                       'row_count INTEGER, description TEXT, columns_json TEXT)')
            for path in sorted((ROOT / 'tables').glob('*/*.csv')):
                columns, types, rows = read_table(path)
                name = table_name(path)
                db.execute(f'CREATE TABLE {quote(name)} (' + ', '.join(
                    f'{quote(c)} {t}' for c, t in zip(columns, types)) + ')')
                db.executemany(f'INSERT INTO {quote(name)} VALUES (' +
                               ','.join('?' for _ in columns) + ')', rows)
                db.execute('INSERT INTO datasets VALUES (?,?,?,?,?)', (
                    name, str(path.relative_to(ROOT)), len(rows),
                    DESCRIPTIONS[path.parent.name], json.dumps(columns)))
            assert db.execute('PRAGMA integrity_check').fetchone()[0] == 'ok'
        os.replace(filename, ROOT / 'results.sqlite')
    finally:
        Path(filename).unlink(missing_ok=True)


if __name__ == '__main__':
    build()
    print('Rebuilt results.sqlite from the published CSV tables.')
