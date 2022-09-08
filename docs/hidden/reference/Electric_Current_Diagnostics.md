# Electric Current Diagnostics

This section configures the electrical current diagnostic settings and is optional. If not present the code will not do any electrical current diagnostics. Starting with r357 the input file format has changed and the code will not work with previous input files. If you are looking for the documentation for previous releases see the [*old version*](#old_version) section below. It accepts the following data:

- **ndump_fac**, integer, default = 0
- **ndump_fac_ave**, integer, default = 0
- **ndump_fac_lineout**, integer, default = 0
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

- "j1", "j2", "j3" - Electric current components
- "div_j" - Divergence of current

Here's an example of a diag_current section that will write diagnostics
information every 20\*ndump iterations for all the components of the
electrical current:

```text
diag_current
{
  ndump_fac = 20,
  reports = "j1", "j2", "j3" , 
}
```

## Old Version

This is the file format used in releases up to r356. All users are urged
to upgrade their input files to the new version as soon as possible.
This section of the documentation will be removed in the near future.

- **ndump_fac_all**, integer, default = 0
- **n_tavg**, integer, default = -1
- **ndump_fac_tavg**, integer, default = 0
- **ndump_fac_ave**, integer, default = 0
- **n_ave**(x_dim), integer, default = 1
- **ifdmp_current**(x_dim), bool, default = .false.
- **ifenv_current**(x_dim), bool, default = .false.
- **ifdmp_div**, bool, default = .false.
- **ndump_fac_lineout**, integer, default = 0
- **lineouts**(:), character(\*), default = "-"

**ndump_fac_all** specifies the frequency at which to write electrical
current diagnostics. This value is multiplied by the *ndump* value
specified in the time_step section to determine the number of iterations
between each diagnostic dump. If set to 0 the writing of this diagnostic
information is disabled. See also *ifdmp_phy_field*.

**n_tavg** specifies the number of time steps to be used when
calculating the time averaged electric current diagnostics specified by
the *ndump_fac_tavg* parameter.

**ndump_fac_tavg** specifies the frequency at which to write time
averaged electric current diagnostics. The electric current is averaged
over *n_tavg* timestep's before being saved to disk. Note that osiris
must hold the temporary information for the average electric current in
memory so this increases the memory requirements for an osiris
simulation. This value is multiplied by the *ndump* value specified in
the *time_step* section to determine the number of iterations between
each diagnostic dump. If set to 0 the writing of this diagnostic
information is disabled. See also *ifdmp_phy_field*.

**ndump_fac_ave** specifies the frequency at which to write spatially
averaged/envelope charge diagnostics. This value is multiplied by the
*ndump* value specified in the *time_step* section to determine the
number of iterations between each diagnostic dump. If set to 0 the
writing of this diagnostic information is disabled. See also
*ifenv_current*.

**n_ave** number of gridpoints on each direction to average over for
spatially averaged dumps.

**ifdmp_current** specifies whether to do an electric current diagnostic
dump for for each direction

**ifenv_current** specifies wheter to take the envelope (.true.) or the
spatial average (.false.) in the specified direction for ave diagnostics
for the electric current.

**ifdmp_div** specifies whether to do an electric current divergence
diagnostic dump

**ndump_fac_lineout** specifies the frequency at which to write lineouts
of the electric current. This value is multiplied by the *ndump value*
specified in the *time_step* section to determine the number of
iterations between each diagnostic dump. If set to 0 the writing of this
diagnostic information is disabled. See also *lineouts*.

**lineouts** specifies which lineouts to extract for the electric
current. Each required lineout is specified as a string in the form "fi,
xi, position" where "fi" can be one of "f1", "f2" or "f3" and specifies
the electric current component to select, "xi" can be one of "x1", "x2"
and "x3" specifying the direction along which to take the lineout, and
"position" specifies the coordinates of the lineout in terms of
perpendicular cell indexes:

- (2D) "f1, x2, 120" - Extract a lineout of J1, along x2, for ix1 = 120.
- (3D) "f2, x1, 100, 70" - Extract a lineout of J2, along x1, for ix2 =
  100, and ix3 = 70.

Lineouts are not available for 1D runs. The file names include the
lineout component, direction and index (the order in which the lineout
was specified) but not the position.

Here's an example of a diag_phy_field section that will write
diagnostics information every `20 * ndump` iterations for all the
components of the electrical current:

```
diag_current
{
  ndump_fac_all = 20,
  ifdmp_current = .true. , 
}
```
