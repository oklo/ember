# Solar grid controller

Run `python3 scripts/solar_grid/coordinate.py` against the existing true-solar campaign. The default campaign directory is `out/solar-grid-z01648-nov23`. Set `EMBER_GRID_ROOT` to use another copy with the same sibling input archives. The controller requires the campaign manifest, recovery plan, executable and physics inputs used for that campaign. The recovery plan contains installation-specific commands and paths. Review those before moving a campaign to another machine.

The coordinator reconciles existing workers and prioritizes runnable masses. When the mass queue is empty it can dispatch bounded recovery jobs from `recovery-plan.json`. An atmosphere job opts into checkpoint recovery through its `resume_column` settings for the native column directory, Python executable, CPU seconds per attempt and maximum restarts. A retry uses a changed, finite profile as its starting guess. A completed native atmosphere still has to pass the independent scientific assessment before the pilot can rejoin the mass queue.

`status.json` reports controller activity and the state of each mass. A working controller can be running only a recovery calculation. Count endpoints from track records separately. When the configured recovery paths are exhausted the controller writes `needs-attention.json`. Creating `PAUSE_NATIVE` in the campaign directory holds new dispatch while existing workers finish.

Run `python3 tests/test_campaign_column_restart.py` to exercise checkpoint recovery in a temporary campaign.
