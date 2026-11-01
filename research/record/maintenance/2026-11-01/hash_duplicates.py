"""Read-only checksum pass over equal-sized inventory files (>= 1 MiB)."""
from concurrent.futures import ThreadPoolExecutor
import hashlib
from pathlib import Path
import sqlite3
import time

ROOT = Path(__file__).resolve().parents[4]
DB = Path(__file__).with_name('inventory.sqlite')


def checksum(row):
    name, size, mtime = row
    path = ROOT / name
    before = path.lstat()
    if path.is_symlink() or (before.st_size, before.st_mtime_ns) != (size, mtime):
        return name, None, 'changed since inventory'
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    after = path.lstat()
    if (after.st_size, after.st_mtime_ns, after.st_ino) != (size, mtime, before.st_ino):
        return name, None, 'changed while hashing'
    return name, digest, None


if __name__ == '__main__':
    connection = sqlite3.connect(DB)
    connection.execute('create table if not exists hashes(path text primary key, sha256 text, error text)')
    rows = connection.execute('''select path,bytes,mtime from files where kind='file'
        and bytes in (select bytes from files where kind='file' and bytes>=1048576
                      group by bytes having count(*)>1)
        and path not in (select path from hashes) order by bytes desc''').fetchall()
    print({'files': len(rows), 'GiB_to_read': sum(r[1] for r in rows)/2**30}, flush=True)
    last = time.monotonic()
    with ThreadPoolExecutor(max_workers=4) as pool:
        for index, result in enumerate(pool.map(checksum, rows), 1):
            connection.execute('insert into hashes values(?,?,?)', result)
            if time.monotonic() - last > 20:
                connection.commit()
                print({'hashed': index, 'total': len(rows)}, flush=True)
                last = time.monotonic()
    connection.commit()
    print('Hash pass complete.', flush=True)
