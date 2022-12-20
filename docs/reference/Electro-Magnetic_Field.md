---
layout: single
classes: wide
title: Electro-Magnetic Field
permalink: /reference/el_mag_fld
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures the Electro-Magnetic field settings and is
optional. If not present the code will use the default values. It
accepts the following data:

- **smooth_type**, character(\*), default = "none"
- **smooth_niter**, integer, default = 1
- **smooth_nmax**, integer, default = -1

<!-- -->

- **type_init_b**, character(\*), default = "uniform"
- **type_init_e**, character(\*), default = "uniform"
- **init_b0**(3), float, default = 0.0
- **init_e0**(3), float, default = 0.0
- **init_b_mfunc**(3), character(\*), default = "NO_FUNCTION_SUPPLIED!"
- **init_e_mfunc**(3), character(\*), default = "NO_FUNCTION_SUPPLIED!"
- **init_dipole_b_m**(3), real, default = 0.0
- **init_dipole_b_x0**(1:x_dim), real, default = 0.0
- **init_dipole_b_r0**, real, default = 0.001
- **init_dipole_e_p**(3), real, default = 0.0
- **init_dipole_e_x0**(1:p_x_dim), real, default = 0.0
- **init_dipole_e_r0**, real, default = 0.001

<!-- -->

- **ext_fld**, character(\*), default = "none"
- **type_ext_b**, character(\*), default = "uniform"
- **type_ext_e**, character(\*), default = "uniform"
- **ext_b0**(3), float, default = 0.0
- **ext_e0**(3), float, default = 0.0
- **ext_b_mfunc**(3), character(\*), default = "NO_FUNCTION_SUPPLIED!"
- **ext_e_mfunc**(3), character(\*), default = "NO_FUNCTION_SUPPLIED!"
- **ext_dipole_b_m**(3), real, default = 0.0
- **ext_dipole_b_x0**(1:x_dim), real, default = 0.0
- **ext_dipole_b_r0**, real, default = 0.001
- **ext_dipole_e_p**(3), real, default = 0.0
- **ext_dipole_e_x0**(1:p_x_dim), real, default = 0.0
- **ext_dipole_e_r0**, real, default = 0.001

<!-- -->

- **solver**, character(\*), default = "yee"

**smooth_type** controls the type of smoothing to be applied to the EM
fields. The actual smoothing parameters are set in a smooth section that
must follow this section. It accepts the following options:

- *"none"* - No smoothing is applied to the EM fields. This is the
  default.
- *"stand"* - The EM fields are smoothed only for particle interpolation.
  The smoothing process will not affect the actual EM fields used in the
  field solver. Smoothing is applied at every time step.
- *"local"* - The smoothing is applied to the simulation EM fields. This
  will affect all EM fields in the simulation and will have a cumulative
  effect, because the EM fields being smoothed in a given iteration will
  have already been smoothed previously. See also **smooth_niter** and
  **smooth_nmax**.
- *"nci"* - Numerical Cherenkonv Instability filter \[1\].
  This is a special filter that is applied in the x1 direction only, and only to fields
  seen by the particles.

**smooth_niter** specifies the number of iterations between smooth
passes for the "local" smooth. The default is 1, meaning the local
smooth will be applied at every time step. **smooth_nmax** specifies the
final iteration at which the "local" smooth will be applied. After this
iteration, the smooth type will be set to "none". The default is -1,
meaning the smoothing will not stop.

**type_init_b** specifies the type of initial field to be used for each
component of the magnetic field. Valid values are:

- *"uniform"* - uniform field (default)
- "math func" - field defined by a user supplied analytical expression.
  This expression can be a function of the spatial coordinates (x1, x2,
  and x3).
- *"dipole"* - dipolar field (not available for 1D simulations).

**type_init_e** specifies the type of initial field to be used for each
component of the electric field. This parameter has the same form as
type_init_b.

**init_b0** specifies the values for an initial magnetic field of type
"uniform".

**init_e0** specifies the values for an initial electric field of type
"uniform".

**init_b_mfunc** specifies the analytical expression to be used for
calculating the initial magnetic field of type "math func". This
expression can be a function of any of the following, which represent
the physical coordinates of the position being calculated: x1 (in 1D, 2D
and 3D), x2 (in 2D and 3D) and x3 (in 3D). See the documentation on the
analytical function parser for details on the mathematical expression.

**init_e_mfunc** specifies the analytical expression to be used for
calculating the initial electric field of type "math func". This
parameter has the same form as ext_b_mfunc.

**init_dipole_b_m** specifies the magnetic dipole moment of the dipole
contributing for an initial magnetic field of type "dipole". The dipole
definition follows the definition given in Jackson, Classical
Electrodynamics, 3rd edition, eq. (5.64), Wiley & Sons, 1998. The
parameters *init_dipole_b_x0* and *init_dipole_b_r0* define the
remaining parameters for the dipole.

**init_dipole_b_x0** specifies the position of the dipole contributing
for an external magnetic field of type "dipole". See the parameters
*init_dipole_b_m* and *init_dipole_b_r0* for further details.

**init_dipole_b_r0** specifies the minimum distance from the dipole
position below which the field contributing for an external magnetic
field of type "dipole" is set to the value at the dipole position. See
the parameters dipole_b_m and dipole_b_x0 for further details.

**init_dipole_e_p** specifies the electric dipole moment of the dipole
contributing for an initial electric field of type "dipole". The dipole
definition follows the definition given in Jackson, Classical
Electrodynamics, 3rd edition, eq. (4.20), Wiley & Sons, 1998. The
parameters *init_dipole_e_x0* and *init_dipole_e_r0* define the
remaining parameters for the dipole.

**init_dipole_e_x0** specifies the position of the dipole contributing
for an initial electric field of type "dipole". See the parameters
*init_dipole_e_p* and *init_dipole_e_r0* for further details.

**init_dipole_e_r0** specifies the minimum distance from the dipole
position below which the field contributing for an initial electric
field of type "dipole" is set to the value at the dipole position. See
the parameters *init_dipole_e_m* and *init_dipole_e_x0* for further
details.

**ext_fld** Controls the use of external electrical and magnetic fields
in the simulation. Valid options are:

- *"none"* - Use only the self generated electric and magnetic fields in
  the simulation. This is the default.
- *"static"* - Use also an externally supplied electric and/or magnetic
  field that is constant in time.
- *"dynamic"* - Same as the previous option, but allows for the external
  field to vary in time (e.g. used in boosted frame simulations). This
  will only affect external fields of type "math func"

**type_ext_b** specifies the type of external field to be used for each
component of the magnetic field. Valid values are:

- *"uniform"* - uniform field (default)
- "math func" - field defined by a user supplied analytical expression.
  This expression can be a function of the spatial coordinates (x1, x2,
  and x3) and, if the ext_fld parameter is set to "dynamic", also of
  time (t).
- *"dipole"* - dipolar field (not available for 1D simulations).

**type_ext_e** specifies the type of external field to be used for each
component of the electric field. This parameter has the same form as
type_ext_b.

**ext_b0** specifies the values for an external magnetic field of type
"uniform".

**ext_e0** specifies the values for an external electric field of type
"uniform".

**ext_b_mfunc** specifies the analytical expression to be used for
calculating the external magnetic field of type "math func". This
expression can be a function of any of the following, which represent
the physical coordinates of the position being calculated: x1 (in 1D, 2D
and 3D), x2 (in 2D and 3D) and x3 (in 3D). For 'dynamic' external fields
the variable t (time) is also allowed. See the documentation on the
analytical function parser for details on the mathematical expression.

**ext_e_mfunc** specifies the analytical expression to be used for
calculating the external electric field of type "math func". This
parameter has the same form as ext_b_mfunc.

**ext_dipole_b_m** specifies the magnetic dipole moment of the dipole
contributing for an external magnetic field of type "dipole". The dipole
definition follows the definition given in Jackson, Classical
Electrodynamics, 3rd edition, eq. (5.64), Wiley & Sons, 1998. The
parameters *ext_dipole_b_x0* and *ext_dipole_b_r0* define the
remaining parameters for the dipole.

**ext_dipole_b_x0** specifies the position of the dipole contributing
for an external magnetic field of type "dipole". See the parameters
*ext_dipole_b_m* and *ext_dipole_b_r0* for further details.

**ext_dipole_b_r0** specifies the minimum distance from the dipole
position below which the field contributing for an external magnetic
field of type "dipole" is set to the value at the dipole position. See
the parameters dipole_b_m and dipole_b_x0 for further details.

**ext_dipole_e_p** specifies the electric dipole moment of the dipole
contributing for an external electric field of type "dipole". The dipole
definition follows the definition given in Jackson, Classical
Electrodynamics, 3rd edition, eq. (4.20), Wiley & Sons, 1998. The
parameters *dipole_e_x0* and *dipole_e_r0* define the remaining
parameters for the dipole.

**ext_dipole_e_x0** specifies the position of the dipole contributing
for an external electric field of type "dipole". See the parameters
*ext_dipole_e_p* and *ext_dipole_e_r0* for further details.

**ext_dipole_e_r0** specifies the minimum distance from the dipole
position below which the field contributing for an external electric
field of type "dipole" is set to the value at the dipole position. See
the parameters *dipole_e_m* and *dipole_e_x0* for further details.

**solver** allows the user to choose from variations of the field
solver. This has an impact on the actual dispersion relation that is
solved. Valid options are:

- *"yee"* - Use the standard Yee solver.
- *"custom"* - Recommended option if not using Yee solver.  Allows for using solvers that eliminate numerical Cherenkov radiation \[2\], remove unphysical fields surrounding relativistic charged particles and reduce numerical dispersion errors in `x1` \[3\], remove dispersion errors and correct for time-stagger errors in the Lorentz force \[4\], or to use an explicitly defined set of coefficeints.  See the [emf_solver](Electro-Magnetic_Field_Solver.md) section for more information.
- *"stencil"* - Use a linear combination of 3 methods for calculating the
  spatial derivatives, using different grid points, with relative
  weights defined by the **k1** and **k2** parameters defined below \[5\].
- *"ndfx"* - Use a Numerical Dispersion Free (along X1 direction) solver.
  **IMPORTANT:** Time step condition is changed to dt \< dx1, not the
  usual Courant condition used for other solvers \[6\].
- *"lehe"* - Use the solver developed by R. Lehe. **IMPORTANT:** time step
  condition is changed to dt \< dx1 \[7\].

Here's an example of a el_mag_fld section that specifies an external
magnetic field of 1.0 normalized units along x3.

```text
el_mag_fld
{
  ext_fld = "static", 
  type_ext_b(3) = "uniform", ! not really necessary since it is the default
  b0(3) = 1.0d0,
}
```

## References

\[1\] [B. B. Godfrey and J.-L. Vay, “Suppressing the numerical Cherenkov instability in FDTD PIC codes,” Journal of Computational Physics, vol. 267, pp. 1–6, Jun. 2014.](https://arxiv.org/abs/1401.0838)

\[2\] [F. Li, et al., “Controlling the numerical Cerenkov instability in PIC simulations using a customized finite difference Maxwell solver and a local FFT based current correction,” Comput Phys Commun, vol. 214, pp. 6-17, 2017.](https://doi.org/10.1016/j.cpc.2017.01.001)

\[3\] [X. Xu, et al., “On numerical errors to the fields surrounding a relativistically moving particle in PIC codes,” J Comput Phys, vol. 413, p. 109451, 2020.](https://doi.org/10.1016/j.jcp.2020.109451)

\[4\] [F. Li, et al., “A new field solver for modeling of relativistic particle-laser interactions using the particle-in-cell algorithm,” Comput Phys Commun, vol. 258, p. 107580, 2021.](https://doi.org/10.1016/j.cpc.2020.107580)

\[5\] [A.D. Greenwood et al., "On the elimination of numerical Cerenkov radiation in PIC simulations", Journal of Computational Physics 201 (2004) 665–684](https://doi.org/10.1016/j.jcp.2004.06.021)

\[6\] [A. Pukhov, "Three-dimensional electromagnetic relativistic particle-in-cell code VLPL (Virtual Laser Plasma Lab)", J. Plasma Physics (1999), vol. 61, part 3, pp. 425–433](https://www.cambridge.org/core/journals/journal-of-plasma-physics/article/threedimensional-electromagnetic-relativistic-particleincell-code-vlpl-virtual-laser-plasma-lab/7FBA476D599E2F19DFEA3F0F2F84FFAB)

\[7\] [R. Lehe, et.al., "Numerical growth of emittance in simulations of laser-wakefield acceleration", Physical Review Special Topics-Accelerators and Beams, vol. 16, no. 2, p. 021301, Feb. 2013](https://doi.org/10.1103/PhysRevSTAB.16.021301)
