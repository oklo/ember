#!/usr/bin/env python3
"""Compare a continuous Hayashi calculation with an interrupted restart."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import shlex
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

        def run(name, extra=(), expected=0, maximum_steps=20, configuration=config):
            command = [str(executable), "--lifetime", str(configuration), "10000", str(work / name),
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

        # Rename the configuration, all family manifests and every payload;
        # preserve only their content and logical roles. Hard links avoid
        # copying the large material dataset and are never modified below.
        relocated = work / "relocated"
        relocated.mkdir()
        settings = dict(shlex.split(line) for line in config.read_text().splitlines()
                        if line.strip() and not line.startswith("#"))
        linked = {}

        def payload(source):
            source = source.resolve(strict=True)
            if source not in linked:
                target = relocated / f"input-{len(linked)}.dat"
                os.link(source, target)
                linked[source] = target
            return linked[source].name

        for role in ("eos", "opacity_low", "opacity_warm", "opacity_bridge", "opacity_hot"):
            source = config.parent / settings[role]
            lines = source.read_text().splitlines()
            first_plane = 4 if role == "eos" else 1
            for i in range(first_plane, len(lines)):
                tokens = shlex.split(lines[i])
                tokens[-1] = json.dumps(payload(source.parent / tokens[-1]))
                lines[i] = " ".join(tokens)
            name = f"renamed-{role}.dat"
            (relocated / name).write_text("\n".join(lines) + "\n")
            settings[role] = name
        for role in ("conduction", "atmosphere", "collisions", "composition"):
            settings[role] = payload(config.parent / settings[role])
        settings.update(zone_threads="1", maximum_step_years="100000000", initial_step_years="10",
                        mass_Msun="0.100")

        def configuration_file(name, changes=None):
            p = relocated / name
            values = dict(settings, **(changes or {}))
            p.write_text("# Moved input dataset; different execution controls.\n" +
                         "".join(f"{k} {json.dumps(v)}\n" for k, v in sorted(values.items())))
            return p

        moved = configuration_file("renamed-configuration.txt")
        run("moved", ("--restart", str(work / "prefix/final.checkpoint")), configuration=moved)
        assert (work / "moved/final.checkpoint").read_bytes() == a.read_bytes()
        assert history(work / "moved") == history(work / "resumed")

        changed = configuration_file("changed-physics.txt", {"energy_tolerance": ".006"})
        rejected = run("physics-change", ("--restart", str(work / "prefix/final.checkpoint")),
                       expected=1, configuration=changed)
        assert "checkpoint" in rejected.stderr and "configuration.energy_tolerance" in rejected.stderr

        # A corrupted table must be rejected by identity, before its parser is
        # constructed. Do not edit a hard link to the retained source dataset.
        (relocated / "broken-conduction.dat").write_text("not a conduction table\n")
        changed = configuration_file("changed-table.txt", {"conduction": "broken-conduction.dat"})
        rejected = run("table-change", ("--restart", str(work / "prefix/final.checkpoint")),
                       expected=1, configuration=changed)
        assert "checkpoint executable or input tables differ: conduction" in rejected.stderr

        # Manifest coordinates matter even when all plane contents agree.
        family = relocated / settings["opacity_low"]
        lines = family.read_text().splitlines()
        tokens = shlex.split(lines[2])
        lines[2] = ".0041 " + json.dumps(tokens[1])
        (relocated / "changed-opacity.dat").write_text("\n".join(lines) + "\n")
        changed = configuration_file("changed-axis.txt", {"opacity_low": "changed-opacity.dat"})
        rejected = run("axis-change", ("--restart", str(work / "prefix/final.checkpoint")),
                       expected=1, configuration=changed)
        assert "checkpoint" in rejected.stderr and "opacity_low.plane.1.metallicity" in rejected.stderr

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
        print("Lifetime restart and relocated/thread-changed replay are exact; changed physics, "
              "table content, family coordinates and damaged checkpoints are rejected.")


if __name__ == "__main__":
    main()
