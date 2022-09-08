# Input file changes

This page lists changes to the input file structure that may cause your
previous input files to stop working. They are listed according to the
code version where the change occurred:

## Revision 599

Collisions now use the global density / omega_p0 value defined in the
simulation section instead of the locally defined norm_charge_density
parameter which has been removed. Here's a quick example of the new
syntax:

```text
 simulation
 {
   n0 = 2.0d23, ! reference density in cm^-3 
 }
 
 collisions
 {
    n_collide = 1,                   ! collide every time step
    nx_collision_cells(1:2) = 10,10, ! each collision cells has 10x10 PIC cells
 
    ! This is now defined in the simulation section
    !norm_charge_density = 2.0d23, ! real density in cm^-3
 
    coulomb_logarithm_automatic = .true.,
  }
```

## Revision 595

Boundary condition types are now defined using text labels in the input
file e.g. "thermal", "conducting". This will most likely break ALL input
files. Here's a quick comparison between old and new syntax for the
particle species boundaries:

```text
spe_bound  {
   ! old version
   ! type(1:2,1) = 50, 50

   ! new version
   type(1:2,1) =   "thermal",  "thermal",

   ! other parameters stay the same
   uth_bnd(1:3,1,1) = 0.0701, 0.0701, 0.0701,
   uth_bnd(1:3,2,1) = 0.0701, 0.0701, 0.0701,
 }
```

Please see the full documentation for the following sections:

- [Electro-Magnetic Field Boundaries](reference/Electro-Magnetic_Field_Boundaries.md)

<!-- -->

- [Species Boundaries](reference/Species_Boundary.md)

Also, please note that there is a [new section](user/Boundary_Conditions.md) on the user guide
discussing the boundary conditions in the code.

## Version 486

The behavior of some zpulse parameters have changed. The focal plane
(per_focus) and laser pulse (lon_start) positions are now defined in
absolute simulation positions (before they were relative to the boundary
the laser is propagating to).
