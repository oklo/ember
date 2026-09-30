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

        # An explicit restart interval changes the controller, not the saved star.
        restart = ("--restart", str(work / "prefix/final.checkpoint"))
        run("reset-step", (*restart, "--restart-step-years", "100"))
        execution = json.loads((work / "reset-step/execution.json").read_text())
        assert execution["restart_step_years"] == 100
        assert execution["effective_initial_step_years"] == 100
        prefix_lines = (work / "prefix/final.checkpoint").read_text().splitlines()
        reset_lines = (work / "reset-step/seed.checkpoint").read_text().splitlines()
        differences = [(a, b) for a, b in zip(prefix_lines, reset_lines) if a != b]
        assert len(differences) == 1
        saved_fields, reset_fields = (line.split() for line in differences[0])
        assert saved_fields[:3] == reset_fields[:3] and saved_fields[4:] == reset_fields[4:]
        assert float(reset_fields[3]) == 100 * 31557600
        endpoint = history(work / "reset-step")[-1]
        assert endpoint["years"] == 10000
        for key in ("luminosity_Lsun", "radius_Rsun", "Teff_K", "H_mass_g", "He3_mass_g"):
            assert abs(endpoint[key] / history(work / "full")[-1][key] - 1) < 1e-4
        rejected = run("step-without-restart", ("--restart-step-years", "100"), expected=1)
        assert "requires --restart" in rejected.stderr
        for index, value in enumerate(("0", "-1", "nan", "inf")):
            run(f"invalid-restart-step-{index}", (*restart, "--restart-step-years", value), expected=1)

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

        # A conductivity prescription is physical input, not a runtime control.
        changed = configuration_file("ionized-conduction.txt", {"conduction_envelope": "ionized"})
        rejected = run("conduction-change", restart, expected=1, configuration=changed)
        assert "checkpoint input selection differs" in rejected.stderr, rejected.stderr
        invalid = configuration_file("invalid-conduction.txt", {"conduction_envelope": "unknown"})
        rejected = run("invalid-conduction", expected=1, configuration=invalid)
        assert "unknown conduction envelope selection" in rejected.stderr

        changed = configuration_file("changed-transport.txt", {
            "transport": "screened_core", "screened_heat_upper_T_K": "3000000",
            "maximum_relative_mixing_gradient": ".01"})
        rejected = run("transport-change", ("--restart", str(work / "prefix/final.checkpoint")),
                       expected=1, configuration=changed)
        assert "checkpoint" in rejected.stderr

        changed = configuration_file("changed-atmosphere-approximation.txt", {
            "atmosphere_metals": "bounded_fixed_Z", "atmosphere_maximum_delta_Z": "1e-8"})
        rejected = run("atmosphere-approximation-change", ("--restart", str(work / "prefix/final.checkpoint")),
                       expected=1, configuration=changed)
        assert "checkpoint" in rejected.stderr and "atmosphere." in rejected.stderr

        changed = configuration_file("changed-atmosphere-overlap.txt", {
            "atmosphere_main_sequence": settings["atmosphere"],
            "atmosphere_join_log_g_low": "4.95", "atmosphere_join_log_g_high": "5.1",
            "atmosphere_join_hydrogen_low": ".6985", "atmosphere_join_hydrogen_high": ".6995"})
        rejected = run("atmosphere-overlap-change", ("--restart", str(work / "prefix/final.checkpoint")),
                       expected=1, configuration=changed)
        assert "checkpoint" in rejected.stderr

        # The later metal-dependent boundary must be selectable from the PMS
        # without querying its unsupported early-temperature/He3 cells. Its
        # response payloads and interval coordinates remain restart identities.
        response_dir = Path(__file__).resolve().parents[1] / "data/atmosphere/lifetime"
        chain_rows = ["EMBER_METAL_ATMOSPHERE_CHAIN 1 2"]
        for z, name in ((".02", "metal_z020_z010.dat"), (".01", "metal_z010_z005.dat")):
            chain_rows.append(z + " " + json.dumps(payload(response_dir / name)))
        chain_path = relocated / "metal-chain.dat"
        chain_path.write_text("\n".join(chain_rows) + "\n")
        chain_settings = {
            "atmosphere_metals": "bounded_fixed_Z", "atmosphere_maximum_delta_Z": "2e-5",
            "atmosphere_main_sequence": settings["atmosphere"],
            "atmosphere_join_log_g_low": "4.95", "atmosphere_join_log_g_high": "5.1",
            "atmosphere_join_hydrogen_low": ".6985", "atmosphere_join_hydrogen_high": ".6995",
            "atmosphere_metal_chain": chain_path.name,
            "atmosphere_metal_join_low": ".01998", "atmosphere_metal_join_high": ".01999",
            "zone_threads": "2", "initial_step_years": "1000"}
        # Use exactly the original invocation's timestep for the comparison.
        initial_settings = dict(shlex.split(line) for line in config.read_text().splitlines()
                                if line.strip() and not line.startswith("#"))
        chain_settings["initial_step_years"] = initial_settings["initial_step_years"]
        chain_config = configuration_file("with-metal-chain.txt", chain_settings)
        run("with-metal-chain", maximum_steps=3, configuration=chain_config)
        assert history(work / "with-metal-chain") == history(work / "prefix")
        restart_chain = ("--restart", str(work / "with-metal-chain/final.checkpoint"))
        renamed_chain = relocated / "renamed-metal-chain.dat"
        renamed_chain.write_bytes(chain_path.read_bytes())
        renamed_config = configuration_file("renamed-chain.txt", dict(chain_settings,
                                           atmosphere_metal_chain=renamed_chain.name))
        run("renamed-chain", restart_chain, configuration=renamed_config)
        assert history(work / "renamed-chain") == history(work / "resumed")
        broken_rows = chain_rows[:]
        # Quoted manifest filenames use the C++ quoted-string convention.
        broken_rows[2] = ".009 " + json.dumps(shlex.split(chain_rows[2])[1])
        (relocated / "changed-chain.dat").write_text("\n".join(broken_rows) + "\n")
        bad_chain = configuration_file("changed-chain.txt", dict(chain_settings,
                                       atmosphere_metal_chain="changed-chain.dat"))
        rejected = run("chain-axis-change", restart_chain, expected=1, configuration=bad_chain)
        assert "checkpoint" in rejected.stderr and "atmosphere.metal_chain" in rejected.stderr

        # Recursive hydrogen intervals must leave the early branch untouched
        # and bind every nested table and coordinate, independent of filenames.
        reference = payload(response_dir / "main_sequence_reference.dat")
        inner = relocated / "hydrogen-inner.dat"
        inner.write_text('EMBER_HYDROGEN_ATMOSPHERE_INTERVAL 1\nreference ' + json.dumps(reference) +
                         '\nchain "metal-chain.dat"\nreference_Z .02\nhydrogen .75 .8\n')
        outer = relocated / "hydrogen-outer.dat"
        outer.write_text('EMBER_HYDROGEN_ATMOSPHERE_INTERVAL 2\nlower_interval "hydrogen-inner.dat"\n'
                         'lower_metal_chain "metal-chain.dat"\nlower_reference_Z .02\nreference ' +
                         json.dumps(reference) + '\nchain "metal-chain.dat"\nreference_Z .02\nhydrogen .82 .85\n')
        interval_settings = dict(chain_settings, atmosphere_hydrogen_interval=outer.name)
        interval_config = configuration_file("hydrogen-interval.txt", interval_settings)
        run("hydrogen-interval", maximum_steps=3, configuration=interval_config)
        assert history(work / "hydrogen-interval") == history(work / "prefix")
        interval_restart = ("--restart", str(work / "hydrogen-interval/final.checkpoint"))
        renamed_outer = relocated / "renamed-interval.dat"
        renamed_outer.write_bytes(outer.read_bytes())
        moved_interval = configuration_file("moved-interval.txt", dict(interval_settings,
                                            atmosphere_hydrogen_interval=renamed_outer.name))
        run("moved-interval", interval_restart, configuration=moved_interval)
        assert history(work / "moved-interval") == history(work / "resumed")
        inner.write_text(inner.read_text().replace("hydrogen .75 .8", "hydrogen .751 .8"))
        rejected = run("changed-inner-interval", interval_restart, expected=1, configuration=moved_interval)
        assert "checkpoint" in rejected.stderr and "hydrogen_interval.lower.hydrogen_low" in rejected.stderr
        inner.write_text(inner.read_text().replace("hydrogen .751 .8", "hydrogen .75 .8"))

        # Very hydrogen-rich sources must also be selectable before they are
        # needed, with all coordinates and leaf bytes bound into the restart.
        envelope_dir = response_dir.parent / "lifetime_hydrogen_envelope"
        envelope_rows = (envelope_dir / "envelope.dat").read_text().splitlines()
        for i in range(1, 5):
            field, filename = shlex.split(envelope_rows[i])
            envelope_rows[i] = field + " " + json.dumps(payload(envelope_dir / filename))
        envelope = relocated / "envelope.dat"
        envelope.write_text("\n".join(envelope_rows) + "\n")
        envelope_settings = dict(interval_settings, atmosphere_hydrogen_envelope=envelope.name)
        envelope_config = configuration_file("envelope.txt", envelope_settings)
        run("envelope", maximum_steps=3, configuration=envelope_config)
        assert history(work / "envelope") == history(work / "prefix")
        envelope_restart = ("--restart", str(work / "envelope/final.checkpoint"))
        renamed_envelope = relocated / "renamed-envelope.dat"
        renamed_envelope.write_bytes(envelope.read_bytes())
        moved_envelope = configuration_file("moved-envelope.txt", dict(envelope_settings,
                                            atmosphere_hydrogen_envelope=renamed_envelope.name))
        run("moved-envelope", envelope_restart, configuration=moved_envelope)
        assert history(work / "moved-envelope") == history(work / "resumed")
        renamed_envelope.write_text(envelope.read_text().replace("gravity_high 5.8 5.9", "gravity_high 5.81 5.9"))
        rejected = run("changed-envelope-axis", envelope_restart, expected=1, configuration=moved_envelope)
        assert "checkpoint" in rejected.stderr and "hydrogen_envelope.gravity_high.low" in rejected.stderr
        changed_rows = envelope_rows[:]
        (relocated / "invalid-envelope-table.dat").write_text("not a physical atmosphere table\n")
        changed_rows[4] = 'high_gravity "invalid-envelope-table.dat"'
        renamed_envelope.write_text("\n".join(changed_rows) + "\n")
        rejected = run("changed-envelope-table", envelope_restart, expected=1, configuration=moved_envelope)
        assert "checkpoint" in rejected.stderr and "hydrogen_envelope.high_gravity" in rejected.stderr

        # A local nonlinear correction tolerance and a global inventory budget
        # measure different quantities; neither must be ordered against the other.
        independent = configuration_file("independent-tolerances.txt", {
            "coupling_abundance_tolerance": "1e-12", "inventory_abundance_tolerance": "1e-14",
            "initial_step_years": "1000"})
        run("independent-tolerances", maximum_steps=1, configuration=independent)
        audit_rows = [json.loads(line) for line in (work / "independent-tolerances/attempts.jsonl").read_text().splitlines()]
        assert audit_rows and audit_rows[-1]["accepted"]

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
