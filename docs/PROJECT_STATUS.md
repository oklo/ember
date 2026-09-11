# Progress and next work — 11 September 2026

The goal is the complete life of an initially 0.1-solar-mass star, through
hydrogen exhaustion, helium-remnant cooling past 100 K and disappearance under
explicit nucleon-decay scenarios. Stable baryons and different environments
are alternatives, not uncertainties hidden inside a single predicted track.

## Current results

- **Evolution:** computed through **3.890 trillion years**, with surface
  temperature **5546 K**, central hydrogen **X = 0.0008157** and central density
  **9993 g/cm³**. All 9169 structures pass the complete-chain and endpoint checks. Central
  hydrogen crosses **X = 0.001** near **3.884 trillion years**.
- **Immediate constraint:** dense interior radiative opacity. The new source
  requires a consistent treatment of the plasma cutoff and spectral averaging;
  more density points alone do not fix its discontinuities. The dense EOS
  extension passes its component checks but awaits acceptance and installation.
- **Atmospheres:** the accepted grid through **5600 K** is already in use.
  Source models through **6400 K** exist; one depth comparison and the remaining
  grid checks still prevent selecting the full extension.
- **Paper:** the [16-page PDF](reports/2026-09-11/ember_status_and_future.pdf)
  now ends with a four-page lifetime timeline. Stellar plots show the
  checked **3.890-trillion-year** track. Later events and decay curves are
  explicitly conditional, with reproducible reference calculations.
- **Parallel research:** Fable has acknowledged the
  [SPH and literature brief](research/fable/ASSIGNMENT.md) and is running
  degenerate-fluid WD encounter pilots with the existing SWIFT code. No encounter result
  is claimed yet. The [coordination protocol](research/fable/COORDINATION.md)
  handles asynchronous reviews, resource allocation and token-limit pauses.

## Priorities

1. Resolve the dense-opacity prescription, accept the dense EOS and continue
   hydrogen exhaustion and the cooling turn.
2. Complete forward atmosphere checks; introduce diffusion and coupled grain
   physics before they control the cooling structure.
3. Follow residual burning, crystallization and thermal evolution through
   1000, 500 and 100 K, then use reduced models where comparisons justify them.
4. In parallel, measure or bound unbound encounter heating with Fable's adaptive
   survey and literature review; connect it to uncertain environmental histories.
5. Evaluate competing nuclear, accretion and dark-matter scenarios, then the
   loss of the remnant under explicit nucleon-decay assumptions.

Use idle compute for useful independent calculations while giving the stellar
trajectory priority. Choose resolution by measured cost and sensitivity, retain
compact results, and match model detail to the
[approximation remit](APPROXIMATION_REMIT.md). The [handoff](../HANDOFF.md)
and Fable's local STATUS.md contain changing job information.
