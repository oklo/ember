"""Check the published bundle without running any stellar or cluster models."""
from pathlib import Path
import ast
import hashlib
import json
import re
import sqlite3
from urllib.parse import unquote

from build_database import ROOT, quote, read_table


def main():
    manifest = json.loads((ROOT / 'manifest.json').read_text())
    for record in manifest['files']:
        path = ROOT / record['path']
        assert path.is_file(), f'Missing: {path}'
        data = path.read_bytes()
        assert len(data) == record['bytes'], f'Size changed: {path}'
        assert hashlib.sha256(data).hexdigest() == record['sha256'], f'Hash changed: {path}'
    expected = {r['path'] for r in manifest['files']} | {'manifest.json'}
    actual = {str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file()
              and 'local' not in p.relative_to(ROOT).parts
              and '__pycache__' not in p.parts and p.suffix != '.pyc'}
    assert expected == actual, f'Unlisted or missing files: {expected ^ actual}'
    for relative in sorted(actual):
        path = ROOT / relative
        if path.suffix == '.py':
            ast.parse(path.read_text(), filename=str(path))
        if path.suffix == '.md':
            for target in re.findall(r'\]\(([^)]+)\)', path.read_text()):
                if target.startswith(('https:', 'http:', '#', 'mailto:')):
                    continue
                target = unquote(target.split('#')[0])
                assert (path.parent / target).exists(), f'Broken link: {relative}: {target}'
    with sqlite3.connect(f'file:{ROOT / "results.sqlite"}?mode=ro', uri=True) as db:
        assert db.execute('PRAGMA integrity_check').fetchone()[0] == 'ok'
        datasets = db.execute('SELECT table_name,csv_path,row_count FROM datasets').fetchall()
        csvs = {str(p.relative_to(ROOT)) for p in (ROOT / 'tables').glob('*/*.csv')}
        assert csvs == {r[1] for r in datasets}, 'Database/CSV table lists differ'
        total = 0
        for name, relative, count in datasets:
            columns, _, rows = read_table(ROOT / relative)
            actual_rows = db.execute('SELECT ' + ','.join(quote(c) for c in columns) +
                                     f' FROM {quote(name)} ORDER BY rowid').fetchall()
            assert rows == actual_rows and len(rows) == count, f'Database mismatch: {name}'
            total += count
    print(f'Verified {len(manifest["files"])} files, links, Python syntax and '
          f'{len(datasets)} database tables ({total} rows).')


if __name__ == '__main__':
    main()
