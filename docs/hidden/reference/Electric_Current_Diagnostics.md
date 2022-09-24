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
