#!/usr/bin/env python3
"""Compare a continuous Hayashi calculation with an interrupted restart."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import tempfile


def history(directory):
    rows = [json.loads(line) for line in (directory / "history.jsonl").read_text().splitlines()]
    for row in rows:
        row.pop("cpu_seconds")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("configuration", type=Path)
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    config = args.configuration.resolve(strict=True)
    with tempfile.TemporaryDirectory(prefix="ember-lifetime-restart-") as temporary:
        work = Path(temporary)

        def run(name, extra=(), expected=0, maximum_steps=20):
            command = [str(executable), "--lifetime", str(config), "10000", str(work / name),
                       "--max-steps", str(maximum_steps), "--cpu-seconds", "120", *extra]
            result = subprocess.run(command, capture_output=True, text=True, timeout=180)
            if result.returncode != expected:
                raise AssertionError(f"{name}: exit {result.returncode}\n{result.stdout}\n{result.stderr}")
            return result

        with ThreadPoolExecutor(max_workers=2) as pool:
            full = pool.submit(run, "full")
            prefix = pool.submit(run, "prefix", maximum_steps=3)
            full.result()
            prefix.result()
        run("resumed", ("--restart", str(work / "prefix/final.checkpoint")))
        a = work / "full/final.checkpoint"
        b = work / "resumed/final.checkpoint"
        assert a.read_bytes() == b.read_bytes(), "restart changes the physical state or next-step decision"
        assert history(work / "full") == history(work / "prefix") + history(work / "resumed")[1:]
        assert history(work / "full")[-1]["years"] == 10000

        # A real restart must reject a damaged model, even if all input-table
        # identities are unchanged. Keep the checkpoint header intact.
        lines = (work / "prefix/final.checkpoint").read_text().splitlines()
        identity_count = int(lines[8])
        first_zone = 10 + identity_count
        row = lines[first_zone].split()
        row[0] = "-1"  # negative enclosed mass
        lines[first_zone] = " ".join(row)
        corrupt = work / "corrupt.checkpoint"
        corrupt.write_text("\n".join(lines) + "\n")
        rejected = run("damaged", ("--restart", str(corrupt)), expected=1)
        assert "checkpoint" in rejected.stderr
        assert not (work / "damaged/final.checkpoint").exists()
        print("Lifetime restart is exact; a damaged checkpoint is rejected.")


if __name__ == "__main__":
    main()
