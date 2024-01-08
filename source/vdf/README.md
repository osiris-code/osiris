# VDF Documentation


## Creating a VDF object

VDF objects must be initialized through a call to the `new` method:

```fortran
type(t_vdf) :: f
call f % new( x_dim, f_dim, nx, gc_num, dx, zero )
```

The parameters are the following:

```fortran
  integer, intent(in) :: x_dim, f_dim
  integer, dimension(:), intent(in) :: nx
  integer, dimension(2,:), intent(in) :: gc_num
  real(p_double), dimension(:), intent(in) :: dx
  logical, intent(in) :: zero
```

* `x_dim` - The number of spatial dimensions of the grid. Must be one of 1, 2 or 3.
* `f_dim` - The number of field components on each point of the grid. Must be >= 1. Set to 1 for a scalar field.
* `gc_num(pos,dir)` - The number of guard cells in each direction, **dir**, for each boundary, **pos** (1 for lower boundary, 2 for upper boundary). Must be >= 0.
* `dx` - cell dimensions for the grid. This is meant to be used by code that requires knowing this quantity, and can be safely ignored (set to 0).
* `zero` - Set this parameter to `.true.` to zero the grid after initialization.

## Creating a VDF object from another VDF object

VDF objects can also be initialized through a call to the `new` method but giving another VDF object as a source parameter:

```fortran
type(t_vdf) :: f, source
call f % new( source )
```

This will create a new VDF object with the exact same structure as the source object. By default the data in the new object is not initialized. Available options are:

```fortran
  class(t_vdf), intent(in) :: source
  integer, intent(in), optional :: f_dim
  logical, intent(in), optional :: zero
  logical, intent(in), optional :: copy
  integer, intent(in), optional :: fc
```

* `source` - Souce VDF to copy structure (and possibly data) from.
* `f_dim` - Number of field components of the new vdf. Defaults to the same number of components as the source.
* `zero` - Controls whether to zero the new vdf. Defaults to `.false.`.
* `copy` - Controls whether to copy the vdf data from source. Defaults to `.false.`.
* `fc` - Controls which field component to copy from source. Defaults copying all components. When set the routine will copy the selected field component into the first component of the newly created VDF object.


## Accessing the data

To access the VDF data you must use the f1, f2 or f3 member variables, depending on whether you are using a 1D, 2D or 3D grid, respectively. These grids will always have 1 more dimension than the number of spatial dimensions: the first component of the grid is used for the field component (use 1 for scalar grids):

```fortran
  ! 1D vdf
  a = vdf_1d % f1( fc, ix )

  ! 2D vdf
  b = vdf_2d % f2( fc, ix, iy )

  ! 3D vdf
  c = vdf_1d % f1( fc, ix, iy, iz )
```

Grid indexes will go from `1`to `f_dim` for the first component, and `1 - gc_num(1,dim)` to `nx(dim) + gc_num(2,dim)` for the spatial dimensions. The data in the range `1 .. nx(dim)` is considerer to belong to the local parallel node, while index values below and above this range are considered as guard cells.

## Writing data to disk

The VDF class provides the `write` method to write VDF data to disk (currently using HDF5). When running in parallel, this routine will gather data from all nodes, and write a single file to disk. The file includes relevant metadata described by a `t_vdf_report` structure:

```fortran
type( t_vdf ) :: data
type( t_vdf_report ) :: info

! ...

! Prepare the metadata
info%name = 'cell_load'
info%ndump = ndump

info%xname  = (/'x1', 'x2', 'x3'/)

info%xlabel = (/'x_1', 'x_2', 'x_3'/)
info%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)

info%time_units = '1 / \omega_p'
info%dt         = 1.0

info%fileLabel = ''
info%path  = trim(path_mass) // 'LOAD' // p_dir_sep // 'CELL'

info%label = 'Particles per cell'
info%units = 'particles'

info%n = n
info%t = t

info%path = trim(path_mass) // 'LOAD' // p_dir_sep // 'GRID'  // p_dir_sep
info%filename = 'cell_load-'// idx_string( n/ndump, p_time_length )

! write the file
call data % write( info, 1, g_space, grid, no_co )
```

Additional diagnostics, such as line / sliceouts and spatial / temporal averages are available in the `m_vdf_report` module.

## Complex data VDF objects

A VDF object that holds complex data is of type `t_cmplx_vdf`. It is initialized exactly the same way as a VDF object:

```fortran
type(t_cmplx_vdf) :: z
call z % new( x_dim, f_dim, nx, gc_num, dx, zero )
```

This is implemented as a subclass of `t_vdf` so instead of the f1, f2, f3 member variables (that point to real values), you must access the complex data using the z1, z2 and z3 member variables:

```fortran
  ! sets the grid point to (1.0 - j)
  cmplx_vdf_2d % f2( 1, 20, 30 ) = cmplx( 1.0, -1.0, kind = p_double )
```

### Parallel Communications

The VDF communication routines defined in module `m_vdf_comm` work with Complex VDF objects, just use them the same way they are used for normal VDF data.

### Digital Filtering (smoothing)

The VDF smoothing routines defined in module `m_vdf_smooth` work with Complex VDF objects, just use them the same way they are used for normal VDF data. However, it should be noted that these have not yet been optimized for complex data, so if your code depends heavily on these you may want to write a new version of the smoothing module.

### Writing complex VDF data to disk

There are two routines available to write complex VDF data to disk. The first is the `write` method described above. When called with a complex VDF object, this method will save the real part of the data.

Additionally, you can use the `write_cmplx` method, that gives the developer control over which part of the complex data to save:

```fortran
type(t_cmplx_vdf) :: z
call z % write_cmplx( report, fc, part, g_space, grid, no_co  )
```

The parameters are the same as the ones used by the `write` method, with an additional `part` parameter:

```fortran
  type( t_vdf_report ), intent(in) :: report   ! File name, path and metadata
  integer,              intent(in) :: fc       ! field component to save
  integer,              intent(in) :: part     ! Part of the data to save
  type( t_space ),      intent(in) :: g_space  ! Spatial information
  class( t_grid ),       intent(in) :: grid     ! Global grid information
  class( t_node_conf ),  intent(in) :: no_co    ! Parallel configuration
```

The **part** parameter allows developer to choose which part of the complex data to save, and may take one of three values:

* _p\_cmplx\_real_ - save the real part of the data
* _p\_cmplx\_imag_ - save the imaginary part of the data
* _p\_cmplx\_abs_ - save the complex magnitude of the data
