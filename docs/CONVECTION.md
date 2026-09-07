# Mixing-length transport

The interior uses local, optically thick Böhm–Vitense mixing-length theory
with the Schwarzschild criterion. `Physics::alpha_mlt` is the ratio of mixing
length to pressure scale height. Its default of 1.9 is a starting value, not
an ember solar calibration.

The coefficient convention is BV58: `(a, b, c) = (1/8, 1/2, 24)` in
[Salaris & Cassisi (2008), equations 1–3 and table 1](https://arxiv.org/pdf/0807.0863).
For pressure scale height `H_P = P/(rho*g)` and mixing length `l = alpha*H_P`,
the buoyancy velocity and enthalpy flux are

```
v^2    = l^2 g delta (grad - grad_element) / (8 H_P)
F_conv = rho cp T v l (grad - grad_element) / (2 H_P).
```

The dimensionless radiative-loss parameter is

```
U = 3 a_rad c_light T^3 / (cp rho^2 kappa l^2) sqrt(8 H_P / (g delta)).
```

Small `U` means efficient convection. Stable layers (`grad_rad <= grad_ad`),
including inward luminosity, retain exactly the radiative gradient.

## A cubic with a bounded root

For unstable layers put `W = grad_rad - grad_ad` and
`q = sqrt(grad - grad_element)`. Combining the element's radiative cooling
with conservation of total flux gives

```
grad_element - grad_ad = 2 U q
W = q^2 + 2 U q + 9 q^3/(8 U).
```

This is algebraically the usual MLT cubic, expressed in a variable that
avoids subtracting `sqrt(grad-grad_ad+U^2) - U`. That subtraction loses the
element contrast when convection is inefficient.

Scale `q` by the smallest of the three single-term roots:

```
q0 = min(sqrt(W), W/(2 U), cbrt(8 U W/9))
t = q/q0
A t^2 + B t + C t^3 = 1.
```

Each coefficient is at most one, with one exactly one, and the positive root
lies between 1/3 and 1. A bracketed Newton iteration solves this monotone
equation. Logarithms form the scale and coefficients without overflowing
`U^2` or `U*W`. Nonconvergence throws.

`ConvectionState` retains the element contrast, superadiabatic excess, and
convective contribution in gradient units separately. Those can remain useful
after the corresponding full gradients round to the same floating-point
number. It also supplies analytic partial derivatives with respect to
`grad_rad`, `grad_ad`, and `ln U`, obtained by implicit differentiation of
the cubic. The analytic zone Jacobian combines these with EOS transport
responses and the state dependence of the pressure scale height and opacity;
see `JACOBIAN.md`. Numerical differences remain a test reference.

## Structure equation and checks

The zone uses centered thermodynamic quantities and always evaluates

```
f_transport = Delta(ln T)/Delta m - grad * Delta(ln P)/Delta m.
```

There is no multiplication of this row by `grad/grad_rad`. A structure test
reaches that ratio below `1e-6` and checks that the response to an outer
temperature perturbation at fixed pressure remains approximately one in
`Delta m * f_transport`.

The convection suite checks a hand-solvable finite-efficiency case, both
asymptotic limits, dimensional flux conservation, element cooling, the
`alpha^-2` scaling of `U`, continuity at the Schwarzschild boundary, and
analytic derivatives against independent differences. A sweep covers
`U = 1e-250 ... 1e250` and radiative excesses from `1e-12` to `1e12`.
Structure tests check both endpoints' transport derivatives in radiative,
superadiabatic, and efficient regimes.

## Scope

This is an interior diffusion prescription. Optically thin radiative losses
and the atmosphere boundary are not implemented here. Neither are Ledoux
composition terms, overshooting, composition mixing, or a geometrical cap on
`l`. The default `alpha` cannot be transferred to another MLT convention or
atmosphere without checking its calibration. Physical inputs to `U` must be
positive and finite; `g=0` is outside the local formula, so the structure
evaluates it at zone midpoints away from the center.
