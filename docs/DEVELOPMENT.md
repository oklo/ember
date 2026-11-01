# Working on Ember

Start with the equations and the scientific question. Keep physical
approximations visible in the configuration and in the resulting record.
Numerical convergence, agreement with an independent calculation and physical
validity are separate claims.

| Responsibility | Entry point |
| --- | --- |
| Structural block elimination | `src/henyey.cpp`, [method](HENYEY.md) |
| Evolution and time-step selection | `src/evolution.cpp`, `src/controller.cpp` |
| Lifetime configuration and restart | `apps/lifetime_driver.cpp`, [configuration](LIFETIME_DRIVER.md) |
| Material table readers and joins | `src/eos_*.cpp`, `src/opacity_table.cpp`, `include/ember/opacity_radiative.hpp` |
| Retained study verification and figures | `research/record/reproduce.py` |
| Offline source generation | `scripts/`, indexed by [physics inputs](DATA_REPRODUCTION.md) |

Use the smallest function that expresses a repeated operation. Name quantities
for their physical or numerical role, with units where the convention is not
obvious. Comments should explain assumptions, equations, coordinate changes
and reasons for numerical choices. Avoid restating assignments or promising
future features. Preserve arithmetic order when a change is intended only to
improve readability.

Run the existing tests that exercise the changed behavior. For material joins,
check values, derivatives, abundance conventions and domain boundaries. For the
Henyey solver, retain the independent dense solve and original-equation residual
checks. Do not replace scientific checks with tests that merely repeat the
implementation.

The dated scripts under `out/` are experiment records. They sometimes depend on
specific source copies or continuation chains. Keep those identities intact;
promote a workflow into maintained code only when its inputs and behavior are
understood. The compact research-record renderer deliberately reads saved data
without importing campaign supervisors.

## Evidence used in the November 2026 cleanup

Cotroneo, Improta and Liguori's [ISSRE study](https://arxiv.org/abs/2508.21634)
found repetition and unused constructs in generated Python/Java samples.
[GitClear's observational report](https://www.gitclear.com/ai_assistant_code_quality_2025_research)
reports increasing duplication and declining code movement in its sample.
[CodeRabbit's PR study](https://www.coderabbit.ai/blog/state-of-ai-vs-human-code-generation-report)
reports readability and dependency problems, while acknowledging uncertain
authorship classification. These are useful review prompts, not a reliable
test of who wrote any particular function or evidence that a coding style
caused a defect.

Documentation evidence is mixed. Guelman and colleagues' [comment study](https://arxiv.org/abs/2408.14007)
found most generated Javadoc equivalent or superior to the originals under
expert review. Katzy and colleagues' [multilingual study](https://arxiv.org/abs/2505.15469)
identified failures of correctness and informativeness that automatic metrics
did not reliably distinguish. The practical criterion here is accurate,
useful explanation, rather than making prose appear less machine-written.

For documentation, the useful standard is whether a reader can find the method,
its limits and a working command. Google's [review guidance](https://google.github.io/eng-practices/review/reviewer/looking-for.html)
supports concise explanations of intent, proportionate abstraction and keeping
documentation consistent with behavior. The cleanup applies those criteria
while retaining the explicit Astra/Opus development lineage.
