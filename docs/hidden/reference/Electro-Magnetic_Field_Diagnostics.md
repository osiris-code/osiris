# Electro-Magnetic Field Diagnostics

This section configures the electro-magnetic field diagnostic settings
and is optional. If not present the code will not do any
electro-magnetic field diagnostics. Starting with r357 the input file
format has changed and the code will not work with previous input files.
If you are looking for the documentation for previous releases see the
[*old version*](#old_version) section below. It accepts the
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
spatial/time averaging, lineouts, etc., as described in the [grid
diagnostics section](:Reference_Guide:_Grid_Diagnostics "wikilink"). The
available quantities are:

- "e1", "e2", "e3" - Electric Field components
- "b1", "b2", "b3" - Magnetic Field components
- "ext_e1", "ext_e2", "ext_e3" - External Electric Field components
  (when applicable)
- "ext_b1", "ext_b2", "ext_b3" - External Magnetic Field components
  (when applicable)
- "part_e1", "part_e2", "part_e3" - Electric Field components seen by
  the particles i.e. smooth(self-generated) + external (when applicable)
- "part_b1", "part_b2", "part_b3" - Magnetic Field components seen by
  the particles i.e. smooth(self-generated) + external (when applicable)
- "ene_e1", "ene_e2", "ene_e3" - Energy in Electric Field components
  (field component squared)
- "ene_b1", "ene_b2", "ene_b3" - Energy in Magnetic Field components
  (field component squared)
- "ene_e" - Total energy in Electric Field ( $E_1^2+E_2^2+E_3^2$ )
- "ene_b" - Total energy in Magnetic Field ( $B_1^2+B_2^2+B_3^2$ )
- "ene_emf" - Total energy in Electro-Magnetic Field ( $E^2+B^2$ )
- "div_e" - Electric Field divergence
- "div_b" - Magnetic Field divergence
- "chargecons" - Charge conservation diagnostic
- "psi" - $\Psi_x$ diagnostic (pseudopotential)
- "s1", "s2", "s3" - Poynting flux components

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

## Old Version

This is the file format used in releases up to r356. All users are urged
to upgrade their input files to the new version as soon as possible.
This section of the documentation will be removed in the near future.

- **ndump_fac_all**, integer, default = 0
- **ndump_fac_ave**, integer, default = 0
- **ndump_fac_ave_ene**, integer, default = 0
- **ndump_fac_charge_cons**, integer, default = 0
- **ndump_fac_ene**, integer, default = 0
- **ndump_fac_ene_int**, integer, default = 0
- **ndump_fac_ext**, integer, default = 0
- **ndump_fac_lineout**, integer, default = 0
- **ndump_fac_part**, integer, default = 0
- **ndump_fac_tavg**, integer, default = 0
- **n_tavg**, integer, default = -1
- **merge_data**, integer, default = 1
- **n_ave**(x_dim), integer, default = 1
- **ifdmp_efl**(3), bool, default = .false.
- **ifdmp_bfl**(3), bool, default = .false.
- **ifenv_efl**(3), bool, default = .false.
- **ifenv_bfl**(3), bool, default = .false.
- **ifene_efl**(3), bool, default = .false.
- **ifene_bfl**(3), bool, default = .false.
- **ifene_tot_efl**, bool, default = .false.
- **ifene_tot_bfl**, bool, default = .false.
- **ifene_tot_emf**, bool, default = .false.
- **ifene_sum_efl**(3), bool, default = .false.
- **ifene_sum_bfl**(3), bool, default = .false.
- **if_dmp_dive**, bool, default = .false.
- **if_dmp_divb**, bool, default = . false.
- **lineouts_b**(:), character(\*), default = "-"
- **lineouts_e**(:), character(\*), default = "-"

**ndump_fac_all** specifies the frequency at which to write
self-generated electro-magnetic field diagnostics. This value is
multiplied by the *ndump* value specified in the *time_step* section to
determine the number of iterations between each diagnostic dump. If set
to 0 the writing of this diagnostic information is disabled. See also
*ifdmp_efl*, *ifdmp_bfl*, *if_dmp_dive*, and *if_dmp_divb*.

**ndump_fac_ave** specifies the frequency at which to write spatially
averaged/envelope self-generated electro-magnetic field diagnostics.
This value is multiplied by the *ndump* value specified in the
*time_step* section to determine the number of iterations between each
diagnostic dump. If set to 0 the writing of this diagnostic information
is disabled. See also *ifenv_efl* and *ifenv_bfl*.

**ndump_fac_ave_ene** specifies the frequency at which to write
spatially averaged/envelope electro-magnetic field energy diagnostics.
This value is multiplied by the *ndump* value specified in the
*time_step* section to determine the number of iterations between each
diagnostic dump. If set to 0 the writing of this diagnostic information
is disabled. See also *ifenv_efl* and *ifenv_bfl*.

**ndump_fac_charge_cons** specifies the frequency at which to write the
charge conservation diagnostic i.e. div E - charge. This value is
multiplied by the *ndump* value specified in the time_step section to
determine the number of iterations between each diagnostic dump.

**ndump_fac_ene** specifies the frequency at which to write
electro-magnetic field energy diagnostics. Note that these diagnostic
only report the energy in the self-generated fields (i.e. the external
fields are not taken into account). This value is multiplied by the
ndump value specified in the time_step section to determine the number
of iterations between each diagnostic dump. If set to 0 the writing of
this diagnostic information is disabled. See also ifene_efl, ifene_bfl,
ifene_tot_efl, ifene_tot_bfl, ifene_sum_efl, ifene_sum_bfl, and
ifene_tot_emf.

**ndump_fac_ene_int** specifies the frequency at which to write
spatially integrated electro-magnetic field energy diagnostics. This
value is multiplied by the *ndump* value specified in the *time_step*
section to determine the number of iterations between each diagnostic
dump. If set to 0 the writing of this diagnostic information is
disabled.

**ndump_fac_ext** specifies the frequency at which to write external
electro-magnetic field diagnostics. This value is multiplied by the
*ndump* value specified in the time_step section to determine the number
of iterations between each diagnostic dump. If set to 0 the writing of
this diagnostic information is disabled. See also *ifdmp_efl* and
*ifdmp_bfl*.

**ndump_fac_part** specifies the frequency at which to write the fields
as seen by simulation particles i.e. smooth(self-generated) + external
(when applicable). This value is multiplied by the *ndump* value
specified in the time_step section to determine the number of iterations
between each diagnostic dump. If set to 0 the writing of this diagnostic
information is disabled. See also *ifdmp_efl* and *ifdmp_bfl*.

**ndump_fac_lineout** specifies the frequency at which to write lineouts
of the electro-magnetic field energy. This value is multiplied by the
*ndump value* specified in the *time_step* section to determine the
number of iterations between each diagnostic dump. If set to 0 the
writing of this diagnostic information is disabled. See also
*lineouts_b* and *lineouts_e*.

**ndump_fac_tavg** specifies the frequency at which to write time
averaged self-generated electro-magnetic field diagnostics. The chosen
electric and magnetic field components are averaged over *n_tavg*
timesteps before being saved to disk. Note that osiris must hold the
temporary information for each field component being reported in memory
so this increases the memory requirements for an osiris simulation. This
value is multiplied by the ndump value specified in the *time_step*
section to determine the number of iterations between each diagnostic
dump. If set to 0 the writing of this diagnostic information is
disabled. See also *ifdmp_efl* and *ifdmp_bfl*.

**n_ave** number of gridpoints on each direction to average over for
spatially averaged dumps.

**n_tavg** specifies the number of time steps to be used when
calculating the time averaged electro-magnetic field diagnostics
specified by the *ndump_fac_tavg* parameter.

**ifdmp_efl** specifies whether to do a raw diagnostic dump for the
electric field for each direction

**ifdmp_bfl** specifies whether to do a raw diagnostic dump for the
magnetic field for each direction

**ifenv_efl** specifies wheter to take the envelope (.true.) or the
spatial average (.false.) in the specified direction for ave diagnostics
for the electric field.

**ifenv_bfl** specifies wheter to take the envelope (.true.) or the
spatial average (.false.) in the specified direction for ave diagnostics
for the magnetic field.

**ifene_efl** specifies whether to do a field energy diagnostic dump for
the electric field for each direction - E(i)^2

**ifene_bfl** specifies whether to do a field energy diagnostic dump for
the electric field for each direction - B(i)^2

**ifene_tot_efl** specifies whether to do a field energy diagnostic dump
for the electric field summed over all directions - E^2 = E1^2 + E2^2 +
E3^2

**ifene_tot_bfl** specifies whether to do a field energy diagnostic dump
for the magnetic field summed over all directions - B^2 = B1^2 + B2^2 +
B3^2

**ifene_tot_emf** specifies whether to do a field energy diagnostic dump
for the full electro-magnetic field - E^2 + B^2

**ifene_sum_efl** specifies whether to do a field energy diagnostic dump
for the electric field summed over some of the directions. Each element
of **ifene_sum_efl** corresponds to E1^2 + E2^2, E1^2 + E3^2, and E2^2 +
E3^2, respectively

**ifene_sum_bfl** specifies whether to do a field energy diagnostic dump
for the magnetic field summed over some of the directions. Each element
of **ifene_sum_efl** corresponds to B1^2 + B2^2, B1^2 + B3^2, and B2^2 +
B3^2, respectively

**ifdmp_dive** specifies whether to do an electric field divergence
diagnostic. The frequency of this diagnostic is controled by the
*ndump_fac_all* parameter.

**ifdmp_divb** specifies whether to do a magnetic field divergence
diagnostic. The frequency of this diagnostic is controled by the
*ndump_fac_all* parameter.

**lineouts_b** specifies which lineouts to extract for the magnetic
field. Each required lineout is specified as a string in the form "fi,
xi, position" where "fi" can be one of "f1", "f2" or "f3" and specifies
the field component to select, "xi" can be one of "x1", "x2" and "x3"
specifying the direction along which to take the lineout, and "position"
specifies the coordinates of the lineout in terms of perpendicular cell
indexes:

- (2D) "f1, x2, 120" - Extract a lineout of B1, along x2, for ix1 = 120.
- (3D) "f2, x1, 100, 70" - Extract a lineout of B2, along x1, for ix2 =
  100, and ix3 = 70.

Lineouts are not available for 1D runs. The file names include the
lineout component, direction and index (the order in which the lineout
was specified) but not the position.

**lineouts_e** specifies which lineouts to extract for the electric
field. See *lineouts_b* for details.

Here's an example of a diag_emf section that will write diagnostics
information every 2\*ndump iterations for all the components of the
electric field:

```
  diag_emf
  {
    ndump_fac_all = 2,
    ifdmp_efl(1:3) = .true. , .true. , .true. ,
  }
```
