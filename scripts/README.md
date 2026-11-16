# Tools

| Task | Entry point |
| --- | --- |
| Supervise a stellar study and recover failed segments | `ember_supervise.py`; [protocol](../docs/AUTONOMY.md) |
| Package a configuration and its material inputs | `package_runtime_data.py` |
| Verify, bundle or install external test fixtures | `package_test_data.py`; [testing](../tests/README.md) |
| Rebuild EOS, opacity, conduction and atmosphere inputs | [Source index](../docs/DATA_REPRODUCTION.md) |
| Generate and recover an individual atmosphere column | `solve_column.py`; [source calculations](../docs/NONGREY.md) |
| Replot retained stellar histories | `../research/record/reproduce.py` |

Source generators, importers and their independent audits belong here.
Use their documented specifications and `--help`; data from different source
versions are not interchangeable. C++ and Fortran probes support source checks
and regression tests. External packages retain their own licenses.

Dated campaign launchers, manuscript builders and unused experiment drivers
are archived in Git history. New study-specific orchestration, raw outputs
and figures belong under `out/`. Promote a tool into this directory when its
inputs, failure behavior and reproduction requirements are explicit.
