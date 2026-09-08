# Central boundary and Henyey relaxation

Ember can now relax a complete model on a fixed Lagrangian mass mesh. The
first verified solution is an explicitly controlled radiative polytrope,
compared with an independent Lane–Emden integration. This is a solver
benchmark, not a physical 0.1 solar-mass stellar model or an evolutionary run.

## Closing the central sphere

The mesh variables include `ln r`, so a mesh point at `r=0` cannot represent
the center. The solver requires

```
0 < m[0] < m[1] < ... < m[N-1] = M.
```

The unresolved central sphere uses the leading regular expansion, with
density and heating evaluated at the first mesh point:

```
f_r = ln r[0] - (ln(3 m[0]/(4 pi)) - ln rho[0])/3
f_L = L[0] - m[0] * (eps_nuc[0] + eps_grav[0]).
```

Its two constraints leave pressure and temperature free for the global
boundary-value problem. They are first-point quantities, approximating their
central limits, rather than additional unknowns at zero radius. Accuracy
must be checked as the central sphere shrinks. For the n=3 benchmark,
the radius residual of an exact Lane–Emden profile is `-xi[0]^2/15` at
leading order; the test verifies this coefficient and its quadratic decay.

`central_residual` supplies an analytic `2x4` Jacobian. The gravitational
heating expression is shared with the zone equations, including the fixed
previous state and the existing backward energy difference. This does not
upgrade the order of that time discretization.

## Block elimination

There are `4N` unknowns, two inner constraints, `4(N-1)` zone equations,
and two surface constraints. `solve_henyey` solves `J dy = -f`:

1. Stack the two carried constraints with four equations from the next
   zone. These six rows contain four current and four next-point variables.
2. Eliminate the four current variables, pivoting across all six active
   equations. Store four triangular relations for back substitution.
3. Pass the remaining two constraints to the next point.
4. At the surface, combine them with the two atmospheric constraints,
   solve the final four equations, and substitute inward.

The block solve takes linear work and storage in mesh size. A boundary need
not have an invertible subblock in a prescribed pair of variables: pivoting
can use the zone equations. Rows are equilibrated before elimination. The
current small-block kernel is portable C++; it does not call Accelerate,
which remains linked and available for future linear algebra.

The correction is checked against the original matrix using the maximum
componentwise backward error

```
|f_i + (J dy)_i| / (|f_i| + sum_j |J_ij dy_j|).
```

If the error exceeds `1e-10`, up to three iterative-refinement solves use
the original coefficients to correct `J*dy+f`. The correction is then checked
again against the original matrix. An error still above `1e-10`, a singular
pivot, or non-finite arithmetic reports failure. No diagonal perturbation or
silent regularization changes the problem. Tests compare with an independent dense Gauss–Jordan solve and a
known solution, force row pivoting, vary equation scales by 160 decades,
and exercise a 2048-point system.

## Newton iteration and units

`relax` copies its input and returns the last accepted model with an explicit
`converged` flag and diagnostics. Mesh and composition are fixed. Non-positive
`dt` is static; positive `dt` requires a previous model on the identical mass
mesh. The routine does not update age or composition or choose a time step.

Luminosity units are fixed from the initial model for the entire solve:

```
L_ref = max_i |L_i|
L_unit[i] = max(|L_i|, L_ref * m[i]/M).
```

There is no solar luminosity floor or division by the evolving local
luminosity. Interior flux may pass through zero. Each luminosity column uses
its local unit; the other columns already represent dimensionless logarithms.
Zone mass, hydrostatic, and transport rows are multiplied by `dm`. Energy
rows use `dm/max(L_unit[i], L_unit[i+1])`, and the inner luminosity constraint
uses `1/L_unit[0]`. These are fixed units, used consistently for both the
Jacobian and line-search merit. The temperature row is never scaled by
convective efficiency.

Newton updates are capped in logarithmic variables and scaled luminosity,
then backtracked until the infinity norm of the residual satisfies an Armijo
decrease. Surface luminosity stays positive. A trial outside a physics
table, or a failed trial atmosphere integration, may trigger backtracking;
it never triggers extrapolation or substitution of another physics module.
Exhausted backtracking returns failure and includes the last rejection.
Invalid input or invalid initial physics throws directly.

Convergence requires both the residual tolerance (default `1e-9`) and the
**undamped** remaining Newton correction tolerance (default `1e-8`). A tiny
damped step alone cannot certify success. The iteration history records
residual, correction, accepted damping, and linear backward error; its final
converged entry has damping zero because no further update is needed.

## Reproducible continuum benchmark

```
cmake --build build
ctest --test-dir build --output-on-failure
mkdir -p out
build/apps/ember-polytrope 256 > out/polytrope-256.json
```

The executable emits JSON containing its assumptions, convergence history,
CGS profile, and radius error. A failed solve exits nonzero. Its fixtures
live in `examples/radiative_polytrope.hpp`, outside the production library.

For a monatomic ideal gas, constant heating `epsilon`, and constant opacity,
`L=epsilon*m` and `P=B*T^4` give `grad_rad=1/4` when

```
B = 4 pi a_rad c G / (3 kappa epsilon).
```

The artificial atmospheric boundary sets `T=Teff` and `P=B*Teff^4`.
Together with the ideal gas law this is an n=3 polytrope,
`P=K*rho^(4/3)`, with `K^3=R_gas^4/B`. An independent RK4 integration of
the Lane–Emden ODE supplies the continuum profile. Total mass is 0.1 solar
masses; the sphere is truncated at `xi=4`, and Stefan–Boltzmann matching
fixes its central density. These artificial inputs are not intended to
represent hydrogen-burning or photospheric physics.

The mesh samples `xi` quadratically from `1e-3` to 4. Starting from smooth
perturbations of the continuum profile, the measured results are:

| Points | Accepted Newton updates | Radius error | Final scaled residual |
|---:|---:|---:|---:|
| 64 | 6 | 3.140% | 7.13e-11 |
| 128 | 7 | 0.8147% | 6.28e-15 |
| 256 | 7 | 0.2049% | 7.06e-15 |

The radius error falls by approximately four per doubling. The tests also
check `L(m)=epsilon*m`, shrinking the central sphere, analytic central
derivatives under degeneracy and contraction, equilibrium preservation with
a fixed previous model, and failure behavior. A small Newton residual is
not a measurement of spatial accuracy; this benchmark keeps both visible.

## Stellar validation and next physics

The actual pp/FD/MLT/grey modules now converge a fixed-composition 0.5 Msun
equilibrium with AESOPUS/OPAL opacity at 128, 256 and 512 points. The finer
meshes required the linear iterative refinement described above. Integrated
nuclear heating balances surface luminosity, and an independently integrated
virial error decreases by about four per mesh doubling. This is a numerical
validation with approximate envelope physics, not yet a realistic M dwarf.
The newer CMS19/AESOPUS/TOPS combination also converges a 0.1 Msun model
through 4096 points, with fine-mesh R/L agreement and virial convergence.
It remains experimental: a local EOS pressure/entropy consistency defect
reaches 24%, and internal energy is explicitly disabled because of a source
join defect. Positive-dt structure calls reject that EOS. See
[`EQUILIBRIUM.md`](EQUILIBRIUM.md) for commands, results and exact limitations.
Physical atmosphere grids, consistent caloric/composition EOS physics,
conduction, adaptive meshes, and composition/time evolution remain pending.
