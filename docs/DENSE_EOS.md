# Dense-matter preparation for the helium remnant

No dense remnant EOS has been installed in the evolved 0.1 Msun calculation.
The present work tests components needed for that installation. It does not
establish a white-dwarf cooling age or cold-atmosphere coverage.

## Cold electron thermal response

The analytic electron component integrates relativistic Fermi-Dirac
occupations. Its original first thermal derivatives subtracted two large
terms to recover a small heat capacity. Pure-He tests at 100 K and
rho=1e3..1e6 g/cm3 showed relative heat-capacity discrepancies up to .233%
against the independently known degenerate limit:

```text
cv = pi^2 n_e k_B^2 T / (rho p_F v_F)
dP/dlnT = rho T cv (2+x_F^2) / [3(1+x_F^2)]
x_F = p_F/(m_e c)
```

See the electron heat-capacity expression in
[Baiko & Yakovlev (2019)](https://www.ioffe.ru/astro/Stars/Paper/baiko_yakovlev19mn.pdf).

`src/fermi.cpp` now pairs quadrature points at opposite offsets from the
Fermi surface when eta>100. Kernel differences and the density-conserving
chemical-potential response are evaluated without subtracting large nearly
equal numbers. This also fixes cancellation in the mixed temperature/density
derivative. These are equivalent evaluations of the same electron integrals,
not a replacement with the leading asymptotic formula. Warm-state arithmetic
is retained. `IdealEos` and `ElectronGas` use the same cold response.

`tests/test_eos.cpp` checks electrons separately from ions, including heat
capacity, thermal pressure, both heat-capacity derivatives and value/response
agreement at twelve states. `results/cold_electron_response_v1_audit.json`
preserves before/after calculations. The largest remaining difference from
the leading Sommerfeld heat capacity is 8.20e-7 at T=10000 K and rho=1000,
where finite-temperature corrections are resolved. At 100 K all four density
tests agree to about 1e-10. The current source also evaluates electron occupation entropy directly on
the cold Fermi shell, and supplies ionic Sackur–Tetrode entropy with species
mixing. First-law and Maxwell tests cover warm and cold regimes, including
the branch join. The earlier response report predates this entropy addition;
the September 10 full-suite validation tests the current implementation.
The production metal-bearing potential EOS remains a separate implementation.

## Ioffe source pressure-derivative audit

The pinned [EOS EIP source](https://www.ioffe.ru/astro/EIP/eipintr.html) is
based on Potekhin & Chabrier and the quantum-ion update of
[Baiko & Chugunov (2022)](https://doi.org/10.1093/mnras/stab3613).
`build_eip_probe.py --phase-probe` builds an isolated wrapper around
`EOSFI22`, at the exact requested density and explicitly selected liquid or
solid phase. It includes ideal ions and ion-ion, ion-electron and electron
exchange/correlation contributions. Ideal electrons and radiation are
excluded, avoiding the native mixture driver's approximate density inversion.

`audit_eip_derivatives.py` independently differences free energy, pressure
and internal energy at three step sizes. It checks entropy, heat capacity,
pressure, both pressure derivatives and F=U-TS. Some forced branches are
metastable or outside their intended physical regime; this tests the
formulae's internal derivatives, not their stable-phase applicability.

Two issues were found:

- `EOSFI22` adds the quantum **temperature** derivative `PDTQL` to its
  density derivative `PDRi`. Substitution of `PDRQL` removes the largest
  defect, but does not fully restore consistency.
- `LIQUBC`'s density derivative is not the full derivative of its own free
  energy when the coefficients C_i depend on the ion density parameter.
  The correction below is independently derived from that free energy;
  it differs from the printed equation (39), as well as the distributed
  implementation. The source's energy, free energy and pressure are retained.

For each quantum mode, put y=C_i(r_s) T_p/T, D_i=dln C_i/dln r_s,
a_i=1/2-D_i/3 and D_i'=dD_i/dln r_s. Since r_s is proportional to
rho^(-1/3), and u_i=y df_i/dy:

```text
p_i = a_i u_i
p_i + dp_i/dlnT   = a_i c_i
p_i + dp_i/dlnrho = (a_i+a_i^2+D_i'/9) u_i - a_i^2 c_i
D_1' = -D_1(1-D_1)
D_2' = 0
D_3' = D_3(3D_1-2D_3-1)
```

The optional `--repair-quantum-pressure-derivatives` build applies both
changes, retaining the upstream source, exact patch and separate executable.
Original and corrected v2 audit reports are in `results/`. Across sixteen
forced-phase states, the maximum density-derivative defect drops from
6.77 to 1.05e-6 in ideal-ion pressure units after Richardson removal of
finite-difference truncation. All other tested identities agree to 1.42e-8
or better in their corresponding ideal-ion units. The residual is recorded,
not treated as zero or as an error bound on total stellar pressure.

No source correction has been submitted upstream or installed into stellar
evolution. Quantum helium phase equilibrium, mixed compositions, caloric
reference matching to the partially ionized EOS, and consistent coexistence
remain outstanding. The native Gamma=175 switch is a documented source
choice; it is not a self-consistent helium melting calculation.
