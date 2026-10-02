"""A relocated runtime package retains every selected envelope input."""
from pathlib import Path
import json
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import package_runtime_data as pack

with tempfile.TemporaryDirectory(prefix="ember-package-test-") as directory:
    root = Path(directory)
    child = root / "opacity-node.dat"
    child.write_text("synthetic opacity payload\n")
    cfg = {}
    for role in pack.INPUT_KEYS:
        path = root / (role + ".dat")
        if role == "eos":
            path.write_text("EMBER_VARIABLE_METAL_HELMHOLTZ_BINARY 1\nfixture\n")
        elif role.startswith("opacity_"):
            path.write_text('EMBER_OPACITY_MIXTURE 1 1 0 1\n0 "opacity-node.dat"\n')
        else:
            path.write_text("synthetic " + role + " payload\n")
        cfg[role] = path.name
    gibbs = ["envelope_gibbs_hydrogen", "envelope_gibbs_helium",
             "envelope_gibbs_hydrogen_warm", "envelope_gibbs_helium_warm"]
    for index, roles in enumerate((["envelope_source"], ["envelope_map"], gibbs)):
        settings = dict(cfg)
        for role in roles:
            path = root / (role + ".dat")
            path.write_text("distinct synthetic envelope payload: " + role + "\n")
            settings[role] = path.name
        config = root / f"configuration-{index}.txt"
        config.write_text("".join(f"{k} {json.dumps(v)}\n" for k, v in settings.items()))
        metadata = root / f"metadata-{index}.json"
        metadata.write_text(json.dumps({k: {"generator": "test fixture", "licence": "synthetic test data"}
                                      for k in [*settings, "configuration"]}))
        destination = root / f"package-{index}"
        pack.package(config, destination / "data", metadata, destination / "manifest.json", root)
        relocated = root / f"relocated-{index}"
        destination.rename(relocated)
        pack.verify(relocated / "manifest.json")
        copied = pack.settings(relocated / "data/configuration.txt")
        for role in roles:
            target = relocated / "data" / copied[role]
            assert not Path(copied[role]).is_absolute()
            assert target.read_bytes() == (root / settings[role]).read_bytes()
        target.write_bytes(target.read_bytes() + b"changed\n")
        try:
            pack.verify(relocated / "manifest.json")
        except ValueError:
            pass
        else:
            raise AssertionError("changed envelope input was not detected")
