#!/usr/bin/env python3
"""Rebuild figures from the published CSV, or recheck local raw histories."""
import argparse
import csv
import gzip
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
FILES = ["evolution-cold-remnant-forward-512-3300gyr-v1.json.gz",
         "evolution-cold-remnant-forward-512-3400gyr-v1.json.gz"]


def archived_history():
    rows = []
    columns = None
    for name in FILES:
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
    assert rows[0][0] == 0 and rows[-1][0] == 3.4e12
    return columns, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from-archives", action="store_true",
                        help="recheck local raw-history joins and regenerate the CSV")
    args = parser.parse_args()
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
        assert rows[0][0] == 0 and rows[-1][0] == 3.4e12
    values = dict(zip(columns, np.asarray(rows).T))
    age = values["age_yr"] / 1e12
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42})
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.0), constrained_layout=True)
    ax = axes[0, 0]
    for key, label, color in [("central_X", "$^1$H", "#ad4d28"),
                              ("central_Y3", "$^3$He", "#2e7a89")]:
        ax.plot(age, values[key], label=label, color=color)
    ax.plot(age, 0.98-values["central_X"]-values["central_Y3"],
            label="$^4$He", color="#625a9c")
    ax.axhline(.2, color="0.55", ls=":", lw=.9)
    ax.set(ylabel="Baryonic mass fraction", ylim=(0, .8))
    ax.legend(frameon=False, ncol=3, fontsize=8)
    axes[0, 1].plot(age, values["L_Lsun"]*1e3, color="#ad4d28")
    axes[0, 1].set(ylabel="$L / (10^{-3} L_\\odot)$")
    axes[1, 0].plot(age, values["Teff_K"], color="#ad4d28")
    axes[1, 0].axhline(3400, color="0.55", ls=":", lw=.9)
    axes[1, 0].set(ylabel="$T_{\\rm eff}$ (K)")
    axes[1, 1].plot(age, values["R_Rsun"], color="#ad4d28")
    axes[1, 1].set(ylabel="$R/R_\\odot$")
    for ax in axes.flat:
        ax.set(xlabel="Age from static initial model (trillion yr)", xlim=(0, 3.4))
        ax.grid(alpha=.15)
    fig.savefig(HERE / "evolution_history.pdf")
    fig.savefig(HERE / "evolution_history.png", dpi=180)
    plt.close(fig)
    source = "local archives with verified joins" if args.from_archives else "published CSV"
    print(f"History from {source}: {len(rows)} states, {age[-1]:.2f} trillion yr")


if __name__ == "__main__":
    main()
