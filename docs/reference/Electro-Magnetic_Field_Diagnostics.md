---
layout: single
classes: wide
title: Electro-Magnetic Field Diagnostics
permalink: /reference/diag_emf
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures the electro-magnetic field diagnostic settings
and is optional. If not present the code will not do any
electro-magnetic field diagnostics. It accepts the
following data:

- **ndump_fac**, integer, default = 0
- **ndump_fac_ave**, integer, default = 0
- **ndump_fac_lineout**, integer, default = 0
- **ndump_fac_ene_int**, integer, default = 0
- **ndump_fac_charge_cons**, integer, default = 0
- **prec**, integer, default = 4
- **n_tavg**, integer, default = -1
- **n_ave**(x_dim), integer, default = -1
- **reports**(:), character(\*), default = "-"

**ndump_fac** - controls the frequency of full grid diagnostics. This
value is multiplied by the *ndump* value specified in the *time_step*
section to determine the number of iterations between each diagnostic
dump. If set to 0 the writing of this diagnostic information is
disabled.

**ndump_fac_ave** - controls the frequency of spatial average / envelope
grid diagnostics. This value is multiplied by the *ndump* value
specified in the *time_step* section to determine the number of
iterations between each diagnostic dump. If set to 0 the writing of this
diagnostic information is disabled.

**ndump_fac_lineout** - controls the frequency of lineout / slice
diagnostics. This value is multiplied by the *ndump* value specified in
the *time_step* section to determine the number of iterations between
each diagnostic dump. If set to 0 the writing of this diagnostic
information is disabled.

**ndump_fac_ene_int** specifies the frequency at which to write
spatially integrated electro-magnetic field energy diagnostics. This
value is multiplied by the *ndump* value specified in the *time_step*
section to determine the number of iterations between each diagnostic
dump. If set to 0 the writing of this diagnostic information is
disabled.

**ndump_fac_charge_cons** specifies the frequency at which to calculate
the error in charge conservation i.e. $F = \nabla \cdot \bf{E} - \rho$.
This value is multiplied by the *ndump* value specified in the time_step
section to determine the number of iterations between each calculation.
When calculated the maximum absolute error is reported to stdout. To
save the spatially resolved F the user must choose the "chargecons"
report below.

**n_tavg** specifies the number of time steps to be used when
calculating the time averaged diagnostics. The frequency of these
diagnostics is controlled by the *ndump_fac* parameter described above.

**n_ave** number of gridpoints on each direction to average over for
spatially averaged dumps. The frequency of these diagnostics is
controlled by the *ndump_fac_ave* parameter described above.

**prec** controls the numerical precision used for grid diagnostics. The
default is to save data in single precision (prec = 4) . If the user
wants data to be saved in double precision this parameter must be set to
8. This option is ignored if OSIRIS is compiled in single precision.

**reports** specifies the grid quantities to report, including
spatial/time averaging, lineouts, etc., as described in the [grid diagnostics section](../other/grid_diagnostics). The
available quantities are:

- *"e1", "e2", "e3"* - Electric Field components
- *"b1" "b2", "b3"* - Magnetic Field components
- *"ext_e1" "ext_e2", "ext_e3"* - External Electric Field components
  (when applicable)
- *"ext_b1" "ext_b2", "ext_b3"* - External Magnetic Field components
  (when applicable)
- *"part_e1" "part_e2", "part_e3"* - Electric Field components seen by
  the particles i.e. smooth(self-generated) + external (when applicable)
- *"part_b1" "part_b2", "part_b3"* - Magnetic Field components seen by
  the particles i.e. smooth(self-generated) + external (when applicable)
- *"ene_e1" "ene_e2", "ene_e3"* - Energy in Electric Field components
  (field component squared)
- *"ene_b1" "ene_b2", "ene_b3"* - Energy in Magnetic Field components
  (field component squared)
- *"ene_e"* - Total energy in Electric Field ( $E_1^2+E_2^2+E_3^2$ )
- *"ene_b"* - Total energy in Magnetic Field ( $B_1^2+B_2^2+B_3^2$ )
- *"ene_emf"* - Total energy in Electro-Magnetic Field ( $E^2+B^2$ )
- *"div_e"* - Electric Field divergence
- *"div_b"* - Magnetic Field divergence
- *"chargecons"* - Charge conservation diagnostic
- *"psi"* - $\Psi_x$ diagnostic (pseudopotential)
- *"s1" "s2", "s3"* - Poynting flux components

Here's an example of a diag_emf section that will write diagnostics
information every 2\*ndump iterations for all the components of the
electric field:

```text
  diag_emf
  {
    ndump_fac = 2,
    reports = "e1", "e2", "e3", "b3, tavg",      
  }
```
