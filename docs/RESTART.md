# Continuing a stellar calculation

The evolution driver can save and restore the complete internal model:
mass mesh, log radius/density/temperature, luminosity, all eight abundances,
abundance basis, elemental inventory, age, accepted/rejected counters, and
the next proposed timestep. Decimal values use enough digits for an exact
double-precision round trip. The initial static relaxation is skipped on restart.

Append `--checkpoint FILE` to the usual nine positional arguments to write
the initial model and then atomically replace that file after each accepted
step. Start with a new output path. Append `--restart FILE` to continue;
the duration argument is the target stellar age, and the saved next timestep
supersedes the initial-step argument. Input and output checkpoint paths must
be distinct because an existing checkpoint output is never overwritten at
the start of a new invocation.

The driver limits each invocation to 10000 accepted steps, stops after more
than 100 consecutive rejected attempts, and requires a time step of at least
one year. A successful step resets the consecutive rejection count. The total
accepted and rejected counts remain in the saved state and output for checking
the full history. Earlier executables instead stopped after 100 rejections
accumulated over the star's entire evolution, which interrupted the trial
near 3.544 trillion years despite continued accepted steps. Changing this rule
does not loosen the accuracy checks or permit a new executable to read an old
checkpoint; the revised calculation starts from its initial model. A second
check in the checkpoint reader/writer still rejected lifetime counts above
the old limits and stopped the refined trial near 3.548 trillion years.
That check has also been corrected. The restart test now verifies identical
physical evolution and saved-state recovery after starting with 15000 accepted
and 102 rejected steps in a test fixture, and rejects negative or overflowing
counters. Research checkpoints are never edited to bypass executable identity.

For example, with the current gas boundary:

```sh
build/apps/ember-evolve 512 2.5e12 1e8 10 sfii-svh wd nongrey:data/atmosphere/nongrey_gs98_z020_extended_tau100.dat metal:data/eos/freeeos300_gs98_z020.dat ledoux-diffusive --checkpoint out/gas-2500gyr.restart > out/gas-2500gyr.json
build/apps/ember-evolve 512 2.75e12 1e8 10 sfii-svh wd nongrey:data/atmosphere/nongrey_gs98_z020_extended_tau100.dat metal:data/eos/freeeos300_gs98_z020.dat ledoux-diffusive --restart out/gas-2500gyr.restart --checkpoint out/gas-2750gyr.restart > out/gas-2750gyr.json
```

Source-domain guards still apply. This continuation has now completed under
the recorded source and executable hashes; see
[the 2.75-trillion-year result](results/evolution_metal_m010_2750gyr_gas.json).
Before continuation, the checkpoint-enabled executable reproduced the full
2.5-trillion-year output byte for byte. Old
JSON profiles are not native checkpoints: they omit internal logarithms,
five metal carrier abundances, and the timestep-controller state.

Restart requires identical physical selections, tolerances, mesh, mass,
executable bytes and original table bytes. The binary may be copied to a new
location. Invoke it through an explicit file path when using checkpoints.
The native FNV-1a fingerprints detect accidental changes; they are not
cryptographic. `scripts/run_evolution_snapshot.py` additionally records
SHA-256 hashes of restart inputs and checkpoint outputs and checks that the
restart input did not change during execution. General code or physics changes require a new calculation. The explicit,
limited hot-opacity extension described below has additional checks.

Restart JSON explicitly identifies its starting age and contains only the
continued history; the final profile is complete. The summary script labels
that interval and its helium-3 maximum accordingly. An endpoint clips the
last timestep, so extending a calculation stopped at an earlier target need
not reproduce the timestep sequence of a run with a later target from the
outset. This remains subject to the ordinary timestep accuracy checks.

For an exact trajectory test, `--checkpoint-after N` writes once after the
Nth accepted step of this invocation while the original calculation keeps
running. `tests/test_evolution_restart.py` continues a 70-Myr checkpoint to
1 Gyr and verifies exact equality of every subsequent accepted history row
and the final profile, including after relocating the executable. It also
checks a final checkpoint round trip and rejection of changed tolerances,
changed table bytes, and truncated input. This validates restart mechanics,
not a stellar lifetime or any physical approximation.

## Adding lower-hydrogen hot opacity

`--opacity-extension-restart FILE` is an explicit continuation mode for adding
only lower-hydrogen planes to the hot TOPS family. It requires
`--restart-source-executable FILE`, `--restart-source-opacity DIR`, and the new
`--opacity-directory DIR`. The original checkpoint is checked against its
original executable and data. Every old hot-table entry and temperature/density
coordinate must remain identical; cool TOPS and AESOPUS files must be unchanged.
EOS, atmosphere, mass, mesh, tolerances and physical selections must also match.
No checkpoint is rewritten to make its identity pass.

This mode does not prove compatibility between executables. That requires an
independent trajectory comparison. The actual 3.600–3.605-trillion-year test has
identical accepted histories and final profiles with the old executable/opacity
and with the new executable/extended opacity. Its
[record](results/opacity_extension_restart_overlap_v1.json) preserves hashes.
The full test suite also rejects changes to old hot or cool entries and verifies
that the extended calculation subsequently supports ordinary exact restart.

The output labels the extension explicitly. A new checkpoint is bound to the
new executable and data. This mode does not permit changing the atmosphere or
EOS; calculations using those new inputs begin from a consistent fresh model.
