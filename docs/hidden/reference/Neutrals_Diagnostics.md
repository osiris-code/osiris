# Neutrals Diagnostics

This section configures the neutral diagnostic settings and is optional.
If not present the code will not do any neutral diagnostics. Starting
with r357 the input file format has changed and the code will not work
with previous input files. If you are looking for the documentation for
previous releases see the [*old version*](#old-version)
section below. It accepts the following data:

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

- "ion_charge" - Background Ion charge density ( particle density times
  ionization level )
- "neut_den" - Density of the initial background neutral gas

Here's an example of a diag_neutral section that will write diagnostics
information every 20\*ndump iterations for ion charge density.

```text
diag_current
{
  ndump_fac = 20,
  reports = "ion_charge" , 
}
```

## Old Version

This section configures neutral diagnostics settings. It accepts the
following data:

- **ionlev_dump_fac**, integer, default = 1
- **neutral_dump_fac**, integer, default = 0
- **merge_data**, integer, default = 1

**ionlev_dump_fac** specifies the dump factor for the ionization level.

**neutral_dump_fac** specifies the dump factor for neutral density.

**merge_data** specifies wheter the code should merge all diagnostics
data on the first node of the simulation. If set to 0 spatially resolved
diagnostics will save data regarding each local simulation volume on
each node.

Here's an example of a diag_neutral section:

```text
diag_neutral
{
  ionlev_dump_fac = 1,
  neutral_dump_fac = 1,
}
```
