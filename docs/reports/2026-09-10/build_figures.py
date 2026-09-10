#!/usr/bin/env python3
"""Rebuild figures from the published CSV, or recheck local raw histories."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import numpy as np

HERE = Path(__file__).resolve().parent
FILES = ["evolution-cold-remnant-x015-transition-512-3560gyr-v2.json.gz"]


def f77_history(work):
    """Extract accepted states, joining the two printed records by model ID."""
    audit = json.loads((HERE.parents[1] / "results" /
                       "f77_core_onset_comparison_v1.json").read_text())
    columns = ["run", "model", "age_yr", "central_X", "central_Y3",
               "L_Lsun", "R_Rsun", "Teff_K", "central_convective_flag"]
    rows, sources = [], []
    for name in ["ajr-20000", "ferguson-20000"]:
        path = work / (name + ".out")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        expected = next(value for key, value in audit["input_files_sha256"].items()
                        if key.endswith("/" + path.name))
        assert digest == expected, f"F77 output changed: {path}"
        lines = path.read_text().splitlines()
        run_rows = []
        for i, line in enumerate(lines):
            if not line.startswith("EMBER_CORE "):
                continue
            core = list(map(float, line.split()[1:]))
            surface = list(map(float, lines[i - 1].split()))
            assert len(core) == 10 and len(surface) == 11
            assert core[0] == surface[0]
            # Surface output has four significant figures; central quantities
            # come from the higher-precision diagnostic on that same state.
            assert abs(surface[1] / core[1] - 1) < 5e-5
            assert abs(surface[10] - core[2]) < max(5e-4 * abs(core[2]), 1e-25)
            run_rows.append([name, int(core[0]), core[1], core[2], core[3],
                             surface[4] * 3.86 / 3.828,
                             surface[5] * 6.96 / 6.957,
                             surface[6], int(core[4])])
        run = next(r for r in audit["runs"] if r["name"] == name)
        assert len(run_rows) == run["recorded_models"]
        assert all(a[2] < b[2] for a, b in zip(run_rows, run_rows[1:]))
        assert run_rows[-1][2] == run["final"]["age_yr"]
        rows.extend(run_rows)
        sources.append({"run": name, "output_sha256": digest,
                        "accepted_states": len(run_rows),
                        "returncode": run["receipt"]["returncode"],
                        "first_core_flag_model": run["central_flag_transitions"][0]["after"]["model"]})
    with (HERE / "f77_history.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)
    provenance = {
        "sources": sources,
        "source_audit": "../../results/f77_core_onset_comparison_v1.json",
        "join": "Accepted core diagnostics joined to the immediately preceding surface record; model IDs, ages and hydrogen fractions checked.",
        "solar_unit_conversion": {"F77_Lsun_erg_s": 3.86e33, "F77_Rsun_cm": 6.96e10,
                                  "figure_Lsun_erg_s": 3.828e33, "figure_Rsun_cm": 6.957e10},
        "precision": "L, R and Teff retain the information in the four-significant-figure F77 output. Unit conversion adds no information. Other columns retain diagnostic precision.",
        "figure_selection": "Age figure shows accepted F77 states from 1 Gyr onward, through the final accepted state. Hydrogen figure shows the common central-H interval, excluding XH above .699. No age shift or extrapolation.",
        "csv_sha256": hashlib.sha256((HERE / "f77_history.csv").read_bytes()).hexdigest(),
    }
    (HERE / "f77_history_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


def archived_history():
    rows = []
    columns = None
    manifest = json.loads((HERE / 'recovery_manifest.json').read_text())
    for name in FILES:
        entry = next(e for e in manifest['entries'] if e['archive'] == 'artifacts/' + name)
        raw = (HERE / 'artifacts' / name).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == entry['gzip_sha256']
        assert hashlib.sha256(gzip.decompress(raw)).hexdigest() == entry['sha256']
        with gzip.open(HERE / "artifacts" / name, "rt") as stream:
            data = json.load(stream)
        assert data["converged"] and data["mass_Msun"] == 0.1
        assert columns is None or columns == data["columns"]
        columns = data["columns"]
        segment = data["history"]
        if rows:
            assert rows[-1][0] == segment[0][0]
            # Initial restart rows reset step diagnostics, but retain the state.
            assert rows[-1][2:7] == segment[0][2:7]
            segment = segment[1:]
        rows.extend(segment)
    assert all(a[0] < b[0] for a, b in zip(rows, rows[1:]))
    assert rows[0][0] == 0 and rows[-1][0] == 3.56e12
    assert len(rows) == 2354
    return columns, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from-archives", action="store_true",
                        help="recheck local raw-history joins and regenerate the CSV")
    parser.add_argument("--f77-work", type=Path,
                        help="re-extract and verify F77 data from the recorded local run directory")
    args = parser.parse_args()
    if args.f77_work:
        f77_history(args.f77_work)
    if args.from_archives:
        columns, rows = archived_history()
        with (HERE / "evolution_history.csv").open("w", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(columns)
            writer.writerows(rows)
    else:
        with (HERE / "evolution_history.csv").open(newline="") as stream:
            reader = csv.reader(stream)
            columns = next(reader)
            # The initial/restored state has no preceding implicit half step.
            # Its gravothermal diagnostic is null in JSON and empty in CSV;
            # retain that missing value rather than manufacturing zero heat.
            rows = [[float('nan') if not value and name == 'last_halfstep_gravothermal_Lsun'
                     else float(value) for name, value in zip(columns,row,strict=True)]
                    for row in reader]
        assert all(a[0] < b[0] for a, b in zip(rows, rows[1:]))
        assert rows[0][0] == 0 and rows[-1][0] == 3.56e12
    values = dict(zip(columns, np.asarray(rows).T))
    with (HERE / "f77_history.csv").open(newline="") as stream:
        f77_rows = list(csv.DictReader(stream))
    provenance = json.loads((HERE / "f77_history_provenance.json").read_text())
    assert hashlib.sha256((HERE / "f77_history.csv").read_bytes()).hexdigest() == provenance["csv_sha256"]
    tracks = [("Ember", values, "#ad4d28", "-")]
    for name, label, color, style in [("ajr-20000", "F77: AJR opacity", "#267c86", "--"),
                                       ("ferguson-20000", "F77: Ferguson opacity", "#7160a1", "-.")]:
        selected = [row for row in f77_rows if row["run"] == name]
        data = {key: np.array([float(row[key]) for row in selected])
                for key in selected[0] if key != "run"}
        assert np.all(np.diff(data["age_yr"]) > 0)
        tracks.append((label, data, color, style))
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42})
    for mode in ["age", "hydrogen"]:
        fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.3))
        fig.subplots_adjust(left=.12, right=.975, bottom=.10, top=.90,
                            wspace=.34, hspace=.35)
        quantities = [("central_X", "Central hydrogen mass fraction", 1),
                      ("L_Lsun", "$L/L_\\odot$", 1),
                      ("Teff_K", "$T_{\\rm eff}$ (K)", 1),
                      ("R_Rsun", "$R/R_\\odot$", 1)] if mode == "age" else [
                      ("L_Lsun", "$L / (10^{-3} L_\\odot)$", 1e3),
                      ("Teff_K", "$T_{\\rm eff}$ (K)", 1),
                      ("R_Rsun", "$R/R_\\odot$", 1),
                      ("central_Y3", "Central $^3$He mass fraction", 1)]
        for label, data, color, style in tracks:
            # Omit only the initial F77 contraction from the age overview.
            mask = np.ones(len(data["age_yr"]), dtype=bool)
            if mode == "age" and label != "Ember":
                mask = data["age_yr"] >= 1e9
            x = data["age_yr"] / 1e12 if mode == "age" else data["central_X"]
            for ax, (key, ylabel, scale) in zip(axes.flat, quantities):
                y = scale * data[key][mask]
                if label == "Ember":
                    ax.plot(x[mask], y, color=color, lw=1.7, zorder=3)
                else:
                    # De-emphasize abrupt changes without deleting or smoothing
                    # data, or claiming they are necessarily numerical errors.
                    if key in ["central_X", "central_Y3"]:
                        abrupt = np.abs(np.diff(y)) > .01
                    else:
                        abrupt = np.abs(np.diff(np.log(y))) > np.log(1.05)
                    # Draw continuous runs as paths, not separate translucent
                    # segments: overlap at thousands of segment ends otherwise
                    # darkens the F77 curves in vector PDF renderers.
                    points = np.column_stack((x[mask], y))
                    starts = np.r_[0, np.flatnonzero(np.diff(abrupt)) + 1]
                    ends = np.r_[starts[1:], len(abrupt)]
                    for faint, alpha in [(False, .30), (True, .08)]:
                        parts = []
                        for start, end in zip(starts, ends):
                            if abrupt[start] == faint:
                                parts.extend([points[start:end+1], [[np.nan, np.nan]]])
                        if parts:
                            path = np.concatenate(parts)
                            ax.plot(path[:, 0], path[:, 1], color=color,
                                    lw=.85, alpha=alpha, zorder=1)
                if mode == "age":
                    marker = "x" if "Ferguson" in label else "o"
                    ax.plot(x[-1], scale * data[key][-1], marker=marker,
                            color=color, alpha=1 if label == "Ember" else .35,
                            ms=4, clip_on=False)
                ax.set_ylabel(ylabel)
        for ax, (key, _, _) in zip(axes.flat, quantities):
            ax.grid(alpha=.15)
            ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.4g}"))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.4g}"))
            if mode == "age":
                ax.set(xlabel="Reported age (trillion yr)", xlim=(0, 3.65))
                if key in ["L_Lsun", "R_Rsun"]:
                    ax.set_yscale("log")
            else:
                ax.set(xlabel="Central hydrogen mass fraction",
                       xlim=(.699, values["central_X"][-1]))
                # Autoscaling must use the displayed abundance interval.
                ys = [d[key][(d["central_X"] <= .699) &
                             (d["central_X"] >= values["central_X"][-1])]
                      for _, d, _, _ in tracks]
                lo, hi = min(y.min() for y in ys), max(y.max() for y in ys)
                scale = dict((k, s) for k, _, s in quantities)[key]
                pad = .08 * (hi - lo)
                ax.set_ylim(scale * (lo - pad), scale * (hi + pad))
        handles = [Line2D([], [], color=color, lw=1.7 if label == "Ember" else .85,
                         alpha=1 if label == "Ember" else .45)
                   for label, _, color, _ in tracks]
        labels = [label for label, _, _, _ in tracks]
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.53, .995),
                   ncol=3, frameon=False, fontsize=8)
        stem = "evolution_history" if mode == "age" else "f77_matched_hydrogen"
        fig.savefig(HERE / (stem + ".pdf"), metadata={"CreationDate": None})
        fig.savefig(HERE / (stem + ".png"), dpi=180)
        plt.close(fig)
    source = "local archives with verified joins" if args.from_archives else "published CSV"
    print(f"Ember from {source}: {len(rows)} states, {values['age_yr'][-1]/1e12:.4g} trillion yr; F77: {len(f77_rows)} accepted states")


if __name__ == "__main__":
    main()
