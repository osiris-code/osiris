# Species Diagnostics

This section configures the particle species diagnostic settings and is
optional. If not present the code will not do any electro-magnetic field
diagnostics. Starting with r357 the input file format has changed and
the code will not work with previous input files. If you are looking for
the documentation for previous releases see the [*old version*](#old-version) section below. It accepts the
following data:

- **ndump_fac**, integer, default = 0
- **ndump_fac_ave**, integer, default = 0
- **ndump_fac_lineout**, integer, default = 0
- **prec**, integer, default = 4
- **n_tavg** , integer, default = 0
- **n_ave**(1:p_x_dim), integer, default =0
- **reports**(:), character(\*), default = "-"
- **rep_cell_avg**(:), character(\*), default = "-"
- **rep_udist**(:), character(\*), default = "-"
- **ndump_fac_ene**, integer, default = 0
- **ndump_fac_temp**, integer, default = 0
- **ndump_fac_raw**, integer, default = 0
- **raw_gamma_limit**, real, default = 0.0
- **raw_fraction**, real, default = 1.0
- **raw_math_expr**, character(\*), default = "NO_FUNCTION_SUPPLIED!"
- **ndump_fac_tracks**, integer, default = 0
- **niter_tracks**, integer, default = 1
- **file_tags**, character(\*), default = ""
- **ifdmp_tracks_efl**(3), bool, default = .false.
- **ifdmp_tracks_bfl**(3), bool, default = .false.
- **ifdmp_tracks_psi**, bool, default = .false.
- **ndump_fac_pha**, integer, default = 0
- **ndump_fac_pha_tavg**, integer, default = 0
- **ps_xmin**(1:p_x_dim), real, default = 0.0
- **ps_xmax**(1:p_x_dim), real, default = 0.0
- **ps_pmin**(1:p_p_dim), real, default = 0.0
- **ps_pmax**(1:p_p_dim), real, default = 0.0
- **if_ps_p_auto**(1:p_p_dim), bool, default = .false.
- **ps_lmin**(1:3), real, default = 0.0
- **ps_lmax**(1:3), real, default = 0.0
- **if_ps_l_auto**(1:3), bool, default = .false.
- **ps_gammamin**, real, default = 1.0
- **ps_gammamax**, real, default = 0.0
- **if_ps_gamma_log**, bool, default = .false
- **if_ps_gamma_auto**, bool, default = .false.
- **ps_nx**(1:p_x_dim), integer, default = 64
- **ps_nx_3D**(1:p_x_dim), integer, default = ps_nx
- **ps_np**(1:p_x_dim), integer, default = 64
- **ps_np_3D**(1:p_x_dim), integer, default = ps_np
- **ps_nl**(1:p_x_dim), integer, default = 64
- **ps_nl_3D**(1:p_x_dim), integer, default = ps_np
- **ps_ngamma**, integer, default = 64
- **n_ene_bins**, integer, default = 0
- **ene_bins**(:), real, default = 0
- **phasespaces**(:), character(\*), default = "-"
- **pha_ene_bin**(:), character(\*), default = "-"
- **pha_cell_avg**(:), character(\*), default = "-"
- **pha_time_avg**(:), character(\*), default = "-"

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

**reports** specifies the spatially resolved quantities to report,
including spatial/time averaging, lineouts, etc., as described in the
[grid diagnostics
section](:Reference_Guide:_Grid_Diagnostics "wikilink"). The available
quantities are:

- "charge" - Species charge
- "m" - Species mass ( charge\*rqm )
- "ene" - Local kinetic energy
- "q1", "q2", "q3" - Heat flux component
- "j1", "j2", "j3" - Electric current component for this species. Please
  note that this current is not exactly the same (but very close) as the
  current used in the OSIRIS algorithm. This current is calculated
  depositing q\*v on a grid and the OSIRIS current is calculated using a
  charge conserving scheme. The differences between the two are however
  minimal.

**rep_cell_avg** specifies the cell averaged quantities to report,
including spatial/time averaging, lineouts, etc., as described in the
[grid diagnostics
section](:Reference_Guide:_Grid_Diagnostics "wikilink"). The available
quantities are the same as for the *reports* item. After depositing the
selected quantity, it will be divided by the absolute cell density.
These diagnostics are controlled by the same parameters as the ones
specified by the *report* item.

**rep_udist** specifies grid resolved fluid / thermal momenta
diagnostics, including spatial/time averaging, lineouts, etc., as
described in the [grid diagnostics
section](:Reference_Guide:_Grid_Diagnostics "wikilink"). The available
quantities are:

- "ufl1", "ufl2", "ufl3" - Average fluid momenta component.
- "uth1", "uth2", "uth3" - Momentum distribution width along the
  specified direction.

**ndump_fac_ene** specifies the frequency at which to write total
particles species energy and momumtum flux diagnostics. This value is
multiplied by the ndump value specified in the time_step section to
determine the number of iterations between each diagnostic dump. If set
to 0 the writing of this diagnostic information is disabled.

**ndump_fac_temp** specifies the frequency at which to write total
particles species temperature diagnostics. This value is multiplied by
the ndump value specified in the time_step section to determine the
number of iterations between each diagnostic dump. If set to 0 the
writing of this diagnostic information is disabled.

**ndump_fac_raw** specifies the frequency at which to write particle
species raw diagnostics. This diagnostic dumps all particle information
(position, momenta and charge) to file. This value is multiplied by the
ndump value specified in the time_step section to determine the number
of iterations between each diagnostic dump. If set to 0 the writing of
this diagnostic information is disabled. See also gamma_limit and
particle_fraction.

**raw_gamma_limit** ( \>= 1.0 ) minimal relativistic gamma for raw
diagnostic. Only particle data from particles with gamma \>= gamma_limit
will be saved. Note that this parameter has nothing to do with the gamma
distribution diagnostic.

**raw_fraction** ( \[0.0, 1.0\] ) fraction of particles to dump for raw
diagnostic. Particles are selected randomly testing for random \<=
particle_fraction so that approximately only particle data from
particle_fraction of the total particles are saved.

**raw_math_expr** specifies a mathematical function for particle
selection in the RAW diagnostic. This expression can be a function of
any of the following:

- x{1\|2\|3} which represent the physical coordinates of the particle.
- p{1\|2\|3} which represent the generalized momenta of the particle.
- g, which represents the particle Lorentz gamma factor.
- t, which represents the current simulation time.

See the documentation on the analytical function parser for details on
the mathematical expression.

**ndump_fac_tracks** specifies the frequency at which to write particle
tracking information to file. This value is multiplied by the ndump
value specified in the time_step section to determine the number of
iterations between each diagnostic dump. All the values for particles
whose tags are in the **file_tags** file will be saved to memory every
**niter_tracks** timesteps. This information will be flushed to disk at
a frequency defined by this parameter. Note that for this diagnostic you
also need to set add_tags in the parent species.

**niter_tracks** specifies the number of iterations between saving a
track point. For example, setting this parameter to 10 would mean that
track points would be added at every 10 timesteps.

**file_tags** specifies the name of the file that holds the list of tags
of the particles we want to follow.

**ifdmp_tracks_efl** specifies whether to save the electric field in the
tracks file for each direction

**ifdmp_tracks_bfl** specifies whether to save the magnetic field in the
tracks file for each direction

**ifdmp_tracks_psi** specifies whether to save the Psi quantity in the
tracks file for each direction

**ndump_fac_pha** specifies the frequency at which to write phasespace
particle species diagnostics. This value is multiplied by the ndump
value specified in the time_step section to determine the number of
iterations between each diagnostic dump. If set to 0 the writing of this
diagnostic information is disabled. The individual phasespace dumps are
turned on and off through the phasespaces, pha_ene_bin and pha_cell_avg
parameters. See these parameters for details.

**ndump_fac_pha_tavg** specifies the frequency at which to write time
averaged phasespace particle species diagnostics. This value is
multiplied by the ndump value specified in the time_step section to
determine the number of iterations between each diagnostic dump. If set
to 0 the writing of this diagnostic information is disabled. The
individual phasespace dumps are turned on and off through the
pha_time_avg parameters. See these parameters for details.

**ps_xmin**, **ps_xmax** specify the lower and upper limit, for every
direction, to be considered for phasespace diagnostics using the x1, x2
or x3 (in 3D) coordinate.

**ps_nx** specifies the number of points, for every direction, to be
considered for 1D or 2D phasespace diagnostics using the x1, x2 or x3
(in 3D) coordinate.

**ps_nx_3D** specifies the number of points, for every direction, to be
considered for 3D phasespace diagnostics using the x1, x2 or x3 (in 3D).
If not specified defaults to the same value as ps_nx. This parameter
allows for these 3D phasespaces to be done at lower resolutions while
maintaining a high resolution for 1D and 2D phasespaces.

**ps_pmin**, **ps_pmax** specify the lower and upper limit, for every
direction, to be considered for 1D and 2D phasespace diagnostics using
the p1, p2 or p3 quantities (linear momentum).

**ps_np** specifies the number of points, for every direction, to be
considered for 1D or 2D phasespace diagnostics using the p1, p2 or p3
quantities.

**ps_np_3D** specifies the number of points, for every direction, to be
considered for 3D phasespace diagnostics using the the p1, p2 or p3
quantities. If not specified defaults to the same value as ps_np. This
parameter allows for these 3D phasespaces to be done at lower
resolutions while maintaining a high resolution for 1D and 2D
phasespaces.

**if_ps_p_auto** specifies whether the code should try to autoscale the
momenta phasespace limits set by ps_pmin and ps_pmax so that all
particles in the species are present in the phasespace diagnostic. The
code will autoscale by increasing the number of cells in the phasespace
diagnostic. The limits set by ps_pmin and ps_pmax will always be present
in the phasespace (the code will never shrik below these parameters).
Note that you can select autoscaling in each momentum direction
independently, each component of if_ps_p_auto corresponds to
autoscaling the momentum in that direction.

**ps_lmin**, **ps_lmax** specify the lower and upper limit, for every
direction, to be considered for 1D and 2D phasespace diagnostics using
the l1, l2 or l3 quantities (angular momentum).

**ps_nl** specifies the number of points, for every direction, to be
considered for 1D or 2D phasespace diagnostics using the l1, l2 or l3
quantities.

**ps_nl_3D** specifies the number of points, for every direction, to be
considered for 3D phasespace diagnostics using the the l1, l2 or l3
quantities. If not specified defaults to the same value as ps_nl. This
parameter allows for these 3D phasespaces to be done at lower
resolutions while maintaining a high resolution for 1D and 2D
phasespaces.

**if_ps_l_auto** specifies whether the code should try to autoscale the
momenta phasespace limits set by ps_lmin and ps_lmax so that all
particles in the species are present in the phasespace diagnostic. The
code will autoscale by increasing the number of cells in the phasespace
diagnostic. The limits set by ps_lmin and ps_lmax will always be present
in the phasespace (the code will never shrik below these parameters).
Note that you can select autoscaling in each momentum direction
independently, each component of if_ps_l_auto corresponds to
autoscaling the angular momentum in that direction.

**ps_gammamin**,**ps_gammamax** specify the lower and upper limit to be
considered for the gamma distribution diagnostics.

**if_ps_gamma_auto** specifies whether the code should try to autoscale
the gamma distribution limits set by ps_gammamin and ps_gammamax so that
all particles in the species are present in this diagnostic. The code
will autoscale by increasing the number of cells in the phasespace
diagnostic. The limits set by ps_gammamin and ps_gammamax will always be
present in the diagnostic (the code will never shrik below these
parameters).

**phasespaces** specifies the phasespace diangnostics to perform. This
parameter is an array of strings, where each string specifies the
parameters for each of the diagnostics to perform. The strings specify
which quantities to use in the phasespace simply by using the quantity
name i.e. "x1" will do a 1D phasespace with the x1 position, "p1x1" will
do a standard 2D momentum/position (p1/x1) phasespace and "p3x2x1" will
do a 3D phasespace with p3, x2 and x1. Valid quantities for use as
phasespace axis are as follows:

- x{1\|2\|3}– particle position along directions 1, 2 and 3
  respectively. Specifying x3 in a 2D run or specifying x3 or x2 in a 1D
  run results in a run-time error.
- p{1\|2\|3}– particle linear momentum along directions 1, 2 and 3
  respectively.
- l{1\|2\|3}– particle angular momentum along directions 1, 2 and 3
  respectively, calculated referenced to the center of the simulation
  box. If using a moving window the angular momentum is calculated
  referencing the center of the current simulation window.
- g, gl - lorentz factor of the particle and logarithm of the lorentz
  factor of the particle.

By default the charge of the particle is used in the phasespace, but
other quantities can be used by appending “_QUANT” to the phase-space
string, where QUANT stands for the quantity being deposited. Quantities
that can be deposited are the folowing:

- charge (default) – deposits the charge of the particle, ex: “x1” or
  “x1_charge”.
- m – deposits the mass of the particle, ex: “x2x1_m”.
- ene – deposits the kinetic energy of the particle, ex: “x1_ene”.
- \|charge\| - deposits the absolute value of the charge of the
  particle, ex: “x3x2x1_\|charge\|”.
- q{1\|2\|3} – deposits the heat flux of the particle in the given
  direction, ex: “x2x1_q2” would deposit the heat flux along x2 in a
  x2x1 phasespace.
- j{1\|2\|3} – deposits the electrical current of the particle in the
  given direction, ex: “x3x1_j1” would deposit the current along x1 in a
  x3x1 phasespace.

The maximum number of phasespaces that a user can specify is defined by
the constant p_max_phasespaces, set at compile time in os-dspec.f90. The
default is 64.

**n_ene_bins** specifies the number of values that the ene_bins
parameter holds. The actual number of energy bins used for the energy
binned phasespace diagnostics defined by the pha_ene_bin parameter are
actually n_ene_bins + 1. See ene_bins and pha_ene_bins for details.

**ene_bins** specifies the energies defining the energy bins to be used
by the energy binned phasespace diagnostics defined by the pha_ene_bin
parameter. The energy for each particle is calculated as gamma - 1
(particle kinetic energy normalized to rest mass) where gamma the
relativistic Lorentz factor. The first bin will hold the particles with
energies \<= than ene_bins(1), the last bin will hold particles with
energies \> ene_bins(n_ene_bins) and the i-th bin will hold particles
with ene_bin(i-1) \< energy \<= ene_bin(i). The results for all bins are
saved in the same file.

**pha_ene_bin** specifies the energy binned phasespace diagnostics to
perform. This parameter has exactly the same structure as the
phasespaces parameter. The difference is that each of these phasespaces
is separated into a set of phasespaces that only take into account the
particles with energies falling into the corresponding energy bin. The
number of energy bins and the ranges for each bin are specified through
the n_ene_bins and ene_bins parameters.

**pha_cell_avg** specifies the cell averaged phasespace diagnostics to
perform. This parameter has exactly the same structure as the
phasespaces parameter. The difference is that each of these phasespaces
stores the average value of the quantity being deposited for every
particle inside each cell (rather than the total sum of this quantity).
If no particles are present inside a given cell the value in that cell
is set to zero.

**pha_time_avg** specifies the time averaged phasespace diagnostics to
perform. This parameter has exactly the same structure as the
phasespaces parameter. The difference is that each of these phasespaces
stores the average value of the phasespace over a number of time steps
specified by the n_tavg parameter. Please note that i) it is not
possible to make time averaged phasespaces if any of the axis used are
being autoranged, as the ranges may change while data is being averaged
and ii) for each time averaged phasespace that the user requests, the
code will need to hold its grid in memory, which means that using these
phasespaces increases the memory requirements for osiris significantly.

Here's an example of a diag_species section:

```text
diag_species 
{
 ndump_fac = 1, 
 reports = "charge",
 
 ndump_fac_pha = 1, 
 ndump_fac_raw = 1,
 ps_xmin(1:2) =  0.0, 0.0, 
 ps_xmax(1:2) = 15.0, 3.5,  
 ps_nx(1:2)   = 4096, 512, 

 ps_pmin(1:3) = -2.0, -2.0, -2.0,
 ps_pmax(1:3) = 40.0,  2.0,  2.0,
 ps_np(1:3)   = 1024,  256,  256,
 if_ps_p_auto(1:3) = .true., .true., .true., 
  
 ps_gammamin = 1.0, 
 ps_gammamax = 50.0,
 ps_ngamma = 8192,
 if_ps_gamma_auto = .true.,
     
 phasespaces = "x1", "gl_m", "p1x1", "x2x1_ene", "p2x1x2",

 raw_gamma_limit = 10.0,
 raw_fraction = 1.0, 
}
```

## Old Version

This is the file format used in releases up to r356. All users are urged
to upgrade their input files to the new version as soon as possible.
This section of the documentation will be removed in the near future.

- **ndump_fac_charge**, integer, default = 0
- **ndump_fac_charge_ave**, integer, default = 0
- **n_ave**(1:p_x_dim), integer, default =0
- **ndump_fac_pha**, integer, default = 0
- **ndump_fac_pha_tavg**, integer, default = 0
- **ndump_fac_ene**, integer, default = 0
- **ndump_fac_raw**, integer, default = 0
- **ndump_fac_charge_lineout**, integer, default = 0
- **ps_xmin**(1:p_x_dim), real, default = 0.0
- **ps_xmax**(1:p_x_dim), real, default = 0.0
- **ps_pmin**(1:p_p_dim), real, default = 0.0
- **ps_pmax**(1:p_p_dim), real, default = 0.0
- **if_ps_p_auto**(1:p_p_dim), bool, default = .false.
- **ps_lmin**(1:3), real, default = 0.0
- **ps_lmax**(1:3), real, default = 0.0
- **if_ps_l_auto**(1:3), bool, default = .false.
- **ps_gammamin**, real, default = 1.0
- **ps_gammamax**, real, default = 0.0
- **if_ps_gamma_log**, bool, default = .false
- **if_ps_gamma_auto**, bool, default = .false.
- **ps_nx**(1:p_x_dim), integer, default = 64
- **ps_nx_3D**(1:p_x_dim), integer, default = ps_nx
- **ps_np**(1:p_x_dim), integer, default = 64
- **ps_np_3D**(1:p_x_dim), integer, default = ps_np
- **ps_nl**(1:p_x_dim), integer, default = 64
- **ps_nl_3D**(1:p_x_dim), integer, default = ps_np
- **ps_ngamma**, integer, default = 64
- **dep_sch**, integer, default = 1
- **raw_gamma_limit**, real, default = 0.0
- **raw_fraction**, real, default = 1.0
- **raw_math_expr**, character(\*), default = "NO_FUNCTION_SUPPLIED!"
- **n_ene_bins**, integer, default = 0
- **ene_bins**(:), real, default = 0
- **phasespaces**(:), character(\*), default = "-"
- **pha_ene_bin**(:), character(\*), default = "-"
- **pha_cell_avg**(:), character(\*), default = "-"
- **pha_time_avg**(:), character(\*), default = "-"
- **n_tavg** , integer, default = 0
- **ndump_fac_tracks**, integer, default = 0
- **niter_tracks**, integer, default = 1
- **file_tags**, character(\*), default = ""
- **charge_lineouts**(:), character(\*), default = "-"
- **ifdmp_tracks_efl**(3), bool, default = .false.
- **ifdmp_tracks_bfl**(3), bool, default = .false.

**ndump_fac_charge** specifies the frequency at which to write charge
density particle species diagnostics. This value is multiplied by the
ndump value specified in the time_step section to determine the number
of iterations between each diagnostic dump. If set to 0 the writing of
this diagnostic information is disabled. Whenever possible this should
be used instead of 'x1' (1D), 'x2x1' (2D), and 'x3x2x1' (3D)
phasespaces, as it uses significantly less memory and respects the code
algorithm, including interpolation type and periodic boundaries.

**ndump_fac_charge_ave** specifies the frequency at which to write
spatially averaged charge density particle species diagnostics. This
value is multiplied by the ndump value specified in the time_step
section to determine the number of iterations between each diagnostic
dump. If set to 0 the writing of this diagnostic information is
disabled. The actual averaging is controlled by the 'n_ave' parameter,
see this parameter for details.

**n_ave** specifies the number of gridpoints on each direction to
average over for spatially averaged charge density dumps. See also the
'ndump_fac_charge_ave' parameter.

**ndump_fac_pha** specifies the frequency at which to write phasespace
particle species diagnostics. This value is multiplied by the ndump
value specified in the time_step section to determine the number of
iterations between each diagnostic dump. If set to 0 the writing of this
diagnostic information is disabled. The individual phasespace dumps are
turned on and off through the phasespaces, pha_ene_bin and pha_cell_avg
parameters. See these parameters for details.

**ndump_fac_pha_tavg** specifies the frequency at which to write time
averaged phasespace particle species diagnostics. This value is
multiplied by the ndump value specified in the time_step section to
determine the number of iterations between each diagnostic dump. If set
to 0 the writing of this diagnostic information is disabled. The
individual phasespace dumps are turned on and off through the
pha_time_avg parameters. See these parameters for details.

**ndump_fac_ene** specifies the frequency at which to write total
particles species energy and momumtum flux diagnostics. This value is
multiplied by the ndump value specified in the time_step section to
determine the number of iterations between each diagnostic dump. If set
to 0 the writing of this diagnostic information is disabled.

**ndump_fac_ene** specifies the frequency at which to write total
particles species energy and momumtum flux diagnostics. This value is
multiplied by the ndump value specified in the time_step section to
determine the number of iterations between each diagnostic dump. If set
to 0 the writing of this diagnostic information is disabled. See also
gamma_limit and particle_fraction.

**ndump_fac_raw** specifies the frequency at which to write particle
species raw diagnostics. This diagnostic dumps all particle information
(position, momenta and charge) to file. This value is multiplied by the
ndump value specified in the time_step section to determine the number
of iterations between each diagnostic dump. If set to 0 the writing of
this diagnostic information is disabled. See also gamma_limit and
particle_fraction.

**ndump_fac_charge_lineout** specifies the frequency at which to write
lineouts of the charge. This value is multiplied by the *ndump value*
specified in the *time_step* section to determine the number of
iterations between each diagnostic dump. If set to 0 the writing of this
diagnostic information is disabled. See also *charge_lineouts*.

**ps_xmin**, **ps_xmax** specify the lower and upper limit, for every
direction, to be considered for phasespace diagnostics using the x1, x2
or x3 (in 3D) coordinate.

**ps_nx** specifies the number of points, for every direction, to be
considered for 1D or 2D phasespace diagnostics using the x1, x2 or x3
(in 3D) coordinate.

**ps_nx_3D** specifies the number of points, for every direction, to be
considered for 3D phasespace diagnostics using the x1, x2 or x3 (in 3D).
If not specified defaults to the same value as ps_nx. This parameter
allows for these 3D phasespaces to be done at lower resolutions while
maintaining a high resolution for 1D and 2D phasespaces.

**ps_pmin**, **ps_pmax** specify the lower and upper limit, for every
direction, to be considered for 1D and 2D phasespace diagnostics using
the p1, p2 or p3 quantities (linear momentum).

**ps_np** specifies the number of points, for every direction, to be
considered for 1D or 2D phasespace diagnostics using the p1, p2 or p3
quantities.

**ps_np_3D** specifies the number of points, for every direction, to be
considered for 3D phasespace diagnostics using the the p1, p2 or p3
quantities. If not specified defaults to the same value as ps_np. This
parameter allows for these 3D phasespaces to be done at lower
resolutions while maintaining a high resolution for 1D and 2D
phasespaces.

**if_ps_p_auto** specifies whether the code should try to autoscale the
momenta phasespace limits set by ps_pmin and ps_pmax so that all
particles in the species are present in the phasespace diagnostic. The
code will autoscale by increasing the number of cells in the phasespace
diagnostic. The limits set by ps_pmin and ps_pmax will always be present
in the phasespace (the code will never shrik below these parameters).
Note that you can select autoscaling in each momentum direction
independently, each component of if_ps_p_auto corresponds to
autoscaling the momentum in that direction.

**ps_lmin**, **ps_lmax** specify the lower and upper limit, for every
direction, to be considered for 1D and 2D phasespace diagnostics using
the l1, l2 or l3 quantities (angular momentum).

**ps_nl** specifies the number of points, for every direction, to be
considered for 1D or 2D phasespace diagnostics using the l1, l2 or l3
quantities.

**ps_nl_3D** specifies the number of points, for every direction, to be
considered for 3D phasespace diagnostics using the the l1, l2 or l3
quantities. If not specified defaults to the same value as ps_nl. This
parameter allows for these 3D phasespaces to be done at lower
resolutions while maintaining a high resolution for 1D and 2D
phasespaces.

**if_ps_l_auto** specifies whether the code should try to autoscale the
momenta phasespace limits set by ps_lmin and ps_lmax so that all
particles in the species are present in the phasespace diagnostic. The
code will autoscale by increasing the number of cells in the phasespace
diagnostic. The limits set by ps_lmin and ps_lmax will always be present
in the phasespace (the code will never shrik below these parameters).
Note that you can select autoscaling in each momentum direction
independently, each component of if_ps_l_auto corresponds to
autoscaling the angular momentum in that direction.

**ps_gammamin**,**ps_gammamax** specify the lower and upper limit to be
considered for the gamma distribution diagnostics set by if_gamma.

**if_ps_gamma_auto** specifies whether the code should try to autoscale
the gamma distribution limits set by ps_gammamin and ps_gammamax so that
all particles in the species are present in this diagnostic. The code
will autoscale by increasing the number of cells in the phasespace
diagnostic. The limits set by ps_gammamin and ps_gammamax will always be
present in the diagnostic (the code will never shrik below these
parameters).

**dep_sch** specifies the deposition scheme used for all the phasespace
and gamma distribution diagnostics. This parameter sets the deposition
scheme to nearest grid point (0) or linear deposition (1). NGP is only
maintained for regression test purposes and should not be used.

**if_gamma** specifies whether to do the gamma distribution diagnostic.
Particles will be deposited onto a 1D grid with ps_ngamma points and
ranging from ps_gammamin and ps_gammamax. Particles with a gamma value
greater than ps_gammamax or smaller than ps_gammamin will not be
deposited. See also the if_ps_gamma_auto parameter.

**raw_gamma_limit** ( \>= 1.0 ) minimal relativistic gamma for raw
diagnostic. Only particle data from particles with gamma \>= gamma_limit
will be saved. Note that this parameter has nothing to do with the gamma
distribution diagnostic.

**raw_fraction** ( \[0.0, 1.0\] ) fraction of particles to dump for raw
diagnostic. Particles are selected randomly testing for random \<=
particle_fraction so that approximately only particle data from
particle_fraction of the total particles are saved.

**raw_math_expr** specifies a mathematical function for particle
selection in the RAW diagnostic. This expression can be a function of
any of the following:

- x{1\|2\|3} which represent the physical coordinates of the particle.
- p{1\|2\|3} which represent the generalized momenta of the particle.
- g, which represents the particle Lorentz gamma factor.
- t, which represents the current simulation time.

See the documentation on the analytical function parser for details on
the mathematical expression.

**phasespaces** specifies the phasespace diangnostics to perform. This
parameter is an array of strings, where each string specifies the
parameters for each of the diagnostics to perform. The strings specify
which quantities to use in the phasespace simply by using the quantity
name i.e. "x1" will do a 1D phasespace with the x1 position, "p1x1" will
do a standard 2D momentum/position (p1/x1) phasespace and "p3x2x1" will
do a 3D phasespace with p3, x2 and x1. Valid quantities for use as
phasespace axis are as follows:

- x{1\|2\|3}– particle position along directions 1, 2 and 3
  respectively. Specifying x3 in a 2D run or specifying x3 or x2 in a 1D
  run results in a run-time error.
- p{1\|2\|3}– particle linear momentum along directions 1, 2 and 3
  respectively.
- l{1\|2\|3}– particle angular momentum along directions 1, 2 and 3
  respectively, calculated referenced to the center of the simulation
  box. If using a moving window the angular momentum is calculated
  referencing the center of the current simulation window.
- g, gl - lorentz factor of the particle and logarithm of the lorentz
  factor of the particle.

By default the charge of the particle is used in the phasespace, but
other quantities can be used by appending “_QUANT” to the phase-space
string, where QUANT stands for the quantity being deposited. Quantities
that can be deposited are the folowing:

- charge (default) – deposits the charge of the particle, ex: “x1” or
  “x1_charge”.
- m – deposits the mass of the particle, ex: “x2x1_m”.
- ene – deposits the kinetic energy of the particle, ex: “x1_ene”.
- \|charge\| - deposits the absolute value of the charge of the
  particle, ex: “x3x2x1_\|charge\|”.
- q{1\|2\|3} – deposits the heat flux of the particle in the given
  direction, ex: “x2x1_q2” would deposit the heat flux along x2 in a
  x2x1 phasespace.
- j{1\|2\|3} – deposits the electrical current of the particle in the
  given direction, ex: “x3x1_j1” would deposit the current along x1 in a
  x3x1 phasespace.

The maximum number of phasespaces that a user can specify is defined by
the constant p_max_phasespaces, set at compile time in os-dspec.f90. The
default is 64.

**n_ene_bins** specifies the number of values that the ene_bins
parameter holds. The actual number of energy bins used for the energy
binned phasespace diagnostics defined by the pha_ene_bin parameter are
actually n_ene_bins + 1. See ene_bins and pha_ene_bins for details.

**ene_bins** specifies the energies defining the energy bins to be used
by the energy binned phasespace diagnostics defined by the pha_ene_bin
parameter. The energy for each particle is calculated as m \* gamma
where m is the particle mass and gamma the relativistic lorentz factor,
and is specified in simulation units. The first bin will hold the
particles with energies \<= than ene_bins(1), the last bin will hold
particles with energies \> ene_bins(n_ene_bins) and the i-th bin will
hold particles with ene_bin(i-1) \< energy \<= ene_bin(i). The results
for all bins are saved in the same file.

**pha_ene_bin** specifies the energy binned phasespace diagnostics to
perform. This parameter has exactly the same structure as the
phasespaces parameter. The difference is that each of these phasespaces
is separated into a set of phasespaces that only take into account the
particles with energies falling into the corresponding energy bin. The
number of energy bins and the ranges for each bin are specified through
the n_ene_bins and ene_bins parameters.

**pha_cell_avg** specifies the cell averaged phasespace diagnostics to
perform. This parameter has exactly the same structure as the
phasespaces parameter. The difference is that each of these phasespaces
stores the average value of the quantity being deposited for every
particle inside each cell (rather than the total sum of this quantity).
If no particles are present inside a given cell the value in that cell
is set to zero.

**n_tavg** specifies the number of time steps to be used when
calculating the phasespaces specified by the pha_time_avg parameter.

**pha_time_avg** specifies the time averaged phasespace diagnostics to
perform. This parameter has exactly the same structure as the
phasespaces parameter. The difference is that each of these phasespaces
stores the average value of the phasespace over a number of time steps
specified by the n_tavg parameter. Please note that i) it is not
possible to make time averaged phasespaces if any of the axis used are
being autoranged, as the ranges may change while data is being averaged
and ii) for each time averaged phasespace that the user requests, the
code will need to hold its grid in memory, which means that using these
phasespaces increases the memory requirements for osiris significantly.

**ndump_fac_tracks** specifies the frequency at which to write particle
tracking information to file. This value is multiplied by the ndump
value specified in the time_step section to determine the number of
iterations between each diagnostic dump. All the values for particles
whose tags are in the **file_tags** file will be saved to memory every
**niter_tracks** timesteps. This information will be flushed to disk at
a frequency defined by this parameter. Note that for this diagnostic you
also need to set add_tags in the parent species.

**niter_tracks** specifies the number of iterations between saving a
track point. For example, setting this parameter to 10 would mean that
track points would be added at every 10 timesteps.

**file_tags** specifies the name of the file that holds the list of tags
of the particles we want to follow.

**charge lineouts** specifies which lineouts to extract for the charge.
Each required lineout is specified as a string in the form "xi,
position" where "xi" can be one of "x1", "x2" and "x3" specifying the
direction along which to take the lineout, and "position" specifies the
coordinates of the lineout in terms of perpendicular cell indexes:

- (2D) "x2, 120" - Extract a lineout of charge along x2, for ix1 = 120.
- (3D) "x1, 100, 70" - Extract a lineout of charge along x1, for ix2 =
  100, and ix3 = 70.

Lineouts are not available for 1D runs. The file names include the
direction and index (the order in which the lineout was specified) but
not the position.

**ifdmp_tracks_efl** specifies whether to save the electric field in the
tracks file for each direction

**ifdmp_tracks_bfl** specifies whether to save the magnetic field in the
tracks file for each direction

Here's an example of a diag_species section:

```text
diag_species 
{
 ndump_fac_charge = 1, 
 ndump_fac_pha = 1, 
 ndump_fac_raw = 1,
 ps_xmin(1:2) =  0.0, 0.0, 
 ps_xmax(1:2) = 15.0, 3.5,  
 ps_nx(1:2)   = 4096, 512, 

 ps_pmin(1:3) = -2.0, -2.0, -2.0,
 ps_pmax(1:3) = 40.0,  2.0,  2.0,
 ps_np(1:3)   = 1024,  256,  256,
 if_ps_p_auto(1:3) = .true., .true., .true., 
  
 ps_gammamin = 1.0, 
 ps_gammamax = 50.0,
 ps_ngamma = 8192,
 if_ps_gamma_auto = .true.,
     
 phasespace = "x1", "gl_m", "p1x1", "x2x1_ene", "p2x1x2",

 raw_gamma_limit = 10.0,
 raw_fraction = 1.0, 
}
```
