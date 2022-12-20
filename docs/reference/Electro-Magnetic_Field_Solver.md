---
layout: single
classes: wide
title: Electro-Magnetic Field Solver
permalink: /reference/emf_solver
usemathjax: true

sidebar:
  nav: "ref"
---

This section (`emf_solver`) configures the Electro-Magnetic field solver settings and is
optional. If not present the code will use the default values.

There are various types of solvers allowed, which are specified in the [el_mag_fld](Electro-Magnetic_Field.md) portion of the input deck. The currently available types are *"yee"*, *"custom"*, *"stencil"*, *"ndfx"* and *"lehe"*, each of which are described below.

## Yee solver

This solver type is selected by setting `solver` to *"yee"* in the [el_mag_fld](Electro-Magnetic_Field.md) portion of the input deck. It does not accept any additional parameters.

## Customized solver

This solver type is selected by setting `solver` to *"custom"* in the [el_mag_fld](Electro-Magnetic_Field.md) portion of the input deck.  It accepts the following data:

- **type**, character(\*), default = "standard"
- **solver_ord**, integer, default = 2

<!-- -->

- **n_coef**, integer, defualt = 1
- **coef_e**(n_coef), real, default = ( 1.0, 0.0, 0.0, ... )
- **coef_b**(n_coef), real, default = ( 1.0, 0.0, 0.0, ... )
- **kl**, real, default = 0.0
- **ku**, real, default = 0.0
- **dk**, real, default = 0.0
- **weight_n**, integer, default = 10
- **weight_w**, integer, default = Inf

<!-- -->

- **filter_current**, logical, default = .false.
- **filter_limit**, real, default = 1.0
- **filter_width**, real, default = 0.0
- **n_damp_cell**, integer, default = 0

<!-- -->

- **correct_current**, logical, default = .false.

<!-- -->

- **taper_bnd**, logical, default = .false.
- **taper_order**, integer, default = 2

**type** - specifies the type of customized solver to be used.  Valid options are:

- *"standard"* - This higher-order finite-difference solver extends the coefficient stencil in the `x1`-direction while retaining the same second-order precision in the `x2`- and `x3`-directions as the regular Yee solver. This higher-order solver has higher precision and better numerical dispersion than Yee.  The following additional parameters should be set: `solver_ord`.
- *"bump"* - This solver is based on the standard higher-order solver and is usually used to suppress the numerical Cherenkov radiation (NCR) emitted by relativistic charged particles moving in the `x1`-direction. The NCR suppression is achieved by slightly modifying the numerical dispersion relation locally \[1\].  The following additional parameters should be set: `solver_ord`, `n_coef`, `kl`, `ku` and `dk`. *Note*: A recommended setting is `solver_ord = 16`, `n_coef = solver_ord`, `kl = 0.1`, `ku = 0.3` and `dk = 0.001`, which should apply to most scenarios. Read Ref. \[1\] for more details of how to choose these parameters.
- *"xu"* - This solver is designed to eliminate the unphysical electromagnetic fields surrounding the relativistic charged particles moving in the `x1`-direction \[2\]. It is also nearly free of numerical dispersion as long as the coefficient stencil is wide enough, i.e., `n_coef` is large enough. The following additional parameters should be set: `solver_ord`, `n_coef`, `weight_n`, and `weight_w`.
- *"dual"* - This solver uses different finite-difference operators for electric and magnetic fields \[3\]. It is designed to eliminate the spurious force felt by the relativistic charged particles moving in the `x1`-direction. This unphysical force stems from the fact that the electric and magnetic fields are alternately staggered on the time axis. It is also nearly free of numerical dispersion as long as the coefficient stencil is wide enough, i.e., `n_coef` is large enough. Some numerical tests indicates that this solver may also efficiently suppress NCR. The following additional parameters should be set: `solver_ord`, `n_coef`, `weight_n`, and `weight_w`.
- *"customized-coef"* - This option allows for the user to manually specify the coefficients for the electric and magnetic field stencils.  The following additional parameters should be set: `n_coef`, `coef_e` and `coef_b`.

**solver_ord** - specifies the order of the solver.  Must be positive and even.

**n_coef** - number of coefficients to use.  This number must be larger than half of the `solver_ord`.

**coef_e**, **coef_b** - arrays containing specified coefficients for the E and B stencils.  This parameter is only considered for the *"customized-coef"* solver type.

**kl**, **ku**, **dk** - specify parameters for the *"bump"* solver.  The `kl` and `ku` parameters are lower and upper limits of the region in `k1` space whose numerical dispersion relation is perturbated. Its value ranges from 0.0 to 0.5, which corresponds to the lower and upper limits of the fundamental Brillouin zone in `k1` space.  The `dk` parameter sets the perturbation amplitude of the `k1` region defined by `kl` and `ku`.

**weight_n**, **weight_w** - specify parameters for the *"xu"* and *"dual"* solvers.  The `weight_n` parameter sets the order of the weight function that is used to fit the desired numerical dispersion relation. The default value is 10 and should not be modified except for experienced users.  The `weight_w` parameter sets the width of the weight function that is used to fit the desired numerical dispersion relation. The default value is 0.3 and should not be modified except for experienced users.

**filter_limit** - TODO: rest of parameters.

## Stencil solver

This solver type is selected by setting `solver` to *"stencil"* in the [el_mag_fld](Electro-Magnetic_Field.md) portion of the input deck.  It accepts the following data:

- **k1**, **k2**, real, default = 0.0

**k1**, **k2** - specify the parameters for the *"stencil"* field solver as defined in \[4\].
Interesting values for these parameters are `k1 \< 0`, and `k2 = 0` or `k2 =
2*k1`. Setting `k1 = k2 = 0` recovers the Yee solver, and setting `k1 =
-1/8`, `k2 = 0`, implements a 4th-order-accurate spatial derivative
approximation.

## Numerical dispersion free in x direction (NDFX) solver

This solver type is selected by setting `solver` to *"ndfx"* in the [el_mag_fld](Electro-Magnetic_Field.md) portion of the input deck.  It is implemented as described in Ref. \[5\] and aims to eliminate dispersion errors along the `x1`-direction.  Note: the time step condition is changed to `dt \< dx1`.  It does not accept any additional parameters.

## Lehe solver

This solver type is selected by setting `solver` to *"lehe"* in the [el_mag_fld](Electro-Magnetic_Field.md) portion of the input deck.  It is implemented as described in Ref. \[6\] and aims to reduce spurious Cherenkov radiation and corresponding emittance growth in accelerated beams.  Note: the time step condition is changed to `dt \< dx1`.  It does not accept any additional parameters.

## References

\[1\] [F. Li, et al., “Controlling the numerical Cerenkov instability in PIC simulations using a customized finite difference Maxwell solver and a local FFT based current correction,” Comput Phys Commun, vol. 214, pp. 6-17, 2017.](https://doi.org/10.1016/j.cpc.2017.01.001)

\[2\] [X. Xu, et al., “On numerical errors to the fields surrounding a relativistically moving particle in PIC codes,” J Comput Phys, vol. 413, p. 109451, 2020.](https://doi.org/10.1016/j.jcp.2020.109451)

\[3\] [F. Li, et al., “A new field solver for modeling of relativistic particle-laser interactions using the particle-in-cell algorithm,” Comput Phys Commun, vol. 258, p. 107580, 2021.](https://doi.org/10.1016/j.cpc.2020.107580)

\[4\] [A. D. Greenwood, et al., "On the elimination of numerical Cerenkov radiation in PIC simulations," J Comput Phys, vol. 201, no. 2, pp. 665-684, 2004.](https://doi.org/10.1016/j.jcp.2004.06.021)

\[5\] [A. Pukhov, “Three-dimensional electromagnetic relativistic particle-in-cell code VLPL (Virtual Laser Plasma Lab),” J Plasma Phys, vol. 61, no. 3, pp. 425-433, 1999.](https://doi.org/10.1017/S0022377899007515)

\[6\] [R. Lehe, et.al., "Numerical growth of emittance in simulations of laser-wakefield acceleration," Phys Rev ST Accel Beams, vol. 16, no. 2, p. 021301, 2013](https://doi.org/10.1103/PhysRevSTAB.16.021301)
