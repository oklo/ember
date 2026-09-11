"""Read the latest plotted history and verify its recorded identity.

The earlier matched-mesh transition keeps its separate CSV and profiles.
"""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def current_history(from_archive=False):
    record = json.loads((HERE / 'evolution_latest_provenance.json').read_text())
    path = HERE / record['history_csv']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == record['history_csv_sha256']
    with path.open(newline='') as stream:
        reader = csv.reader(stream)
        columns = next(reader)
        rows = [[None if value == '' else float(value) for value in row]
                for row in reader]
    assert len(rows) == record['states']
    assert rows[0][0] == 0 and rows[-1][0] == record['endpoint']['age_yr']
    assert all(a[0] < b[0] for a, b in zip(rows, rows[1:]))
    for row in rows:
        for key, value in zip(columns, row, strict=True):
            assert (key == 'last_halfstep_gravothermal_Lsun' if value is None
                    else math.isfinite(value))
    if from_archive:
        packed = (HERE / record['raw_archive']).read_bytes()
        assert hashlib.sha256(packed).hexdigest() == record['raw_archive_sha256']
        raw = gzip.decompress(packed)
        assert hashlib.sha256(raw).hexdigest() == record['raw_sha256']
        data = json.loads(raw)
        assert data['mass_Msun'] == .1 and data['points'] == 512
        assert data['converged'] == record['requested_age_reached']
        assert data['columns'] == columns and data['history'] == rows
    return columns, rows
