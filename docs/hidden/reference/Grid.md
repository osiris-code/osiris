# Grid

This section configures the numerical grid and coordinate settings and
must be present in the input file. It accepts the following data:

- **nx_p**(x_dim), integer, default = 1
- **coordinates**, character(\*), default = "cartesian"
- **load_balance**(x_dim), logical, default = .false.
- **lb_type**, character(\*), default = "none"
- **lb_gather**, character(\*), default = "sum"
- **n_dynamic**, integer, default = 0
- **start_load_balance**, integer, default = -1
- **max_imbalance**, float, default = 0.0
- **cell_weight**, float, default = 0.0
- **ndump_global_load**, integer, default = 0

<!-- -->

- **ndump_node_load**, integer, default = 0
- **ndump_grid_load**, integer, default = 0

**nx_p** specifies the number of grid cells to use in each direction for
the simulation.

**coordinates** specifies the coordinate system to use and must be one
of the following:

- "cartesian" - use cartesian coordinates.
- "cylindrical" - use cylindrical coordinates with B1 defined on the
  symmetry axis.

**load_balance** - specifies in which directions the code will attempt
to improve load balance by shifting node boundaries. When choosing any
load balance type (other than none) at least 1 direction must be set to
.true.

**lb_type** specifies the type of parallel load balancing to use for
distributing the computational load among the multiple nodes specified.
This parameter may take one of the following values:

- "none" - Divide the grid cells as evenly as possible among all nodes.
- "static" - Divide the computational load as evenly as possible along
  the nodes in the directions specified by *load_balance* in the
  beginning of the simulation.
- "dynamic" - Divide the computational load as evenly as possible along
  the nodes in the directions specified by *load_balance* in the
  beginning of the simulation, and then dynamically adjust the load
  distribution during the simulation at a frequency defined by the
  *n_dynamic* parameter.

**lb_gather** specifies how the load along the load balance direction is
gathered from all nodes. Possible values are:

- "sum" - Sum data from nodes in transverse directions. This is the
  default.
- "max" - Use the maximum value from nodes in transverse direction. In
  situations where the load varies significantly transversely this may
  yield a better partition.

**n_dynamic** When using the *dynamic* load balance type specifies at
which frequency the code should try to adjust node boundaries in order
to improve load balance.

**start_load_balance** When using the *dynamic* load balance, this
parameter specifies the iteration at which dynamic load balance will
start. The default is -1 meaning that dynamic load balancing will begin
at simulation start.

**max_imbalance** specifies a threshold below which the code will not
reshape node boundaries when doing dynamic load balancing. For example,
if the user sets this parameter to 0.1, the code will not change node
boundaries if the difference between maximum load and minimum load on
all the nodes is below 10%. This parameter must be in the range
\[0..1\[.

**cell_weight** specifies the relative weight of grid cells in relation
to particles to be considered for load calculations. Setting this
parameter to 0. means that cells are not considered in load calculations
and only particles are considered. Setting it to 1.5 means that each
cell will be counted as 1.5 particles.

**ndump_global_load** specifies the frequency at which to save global
load information: total, max, min and average particles per node. Output
is saved on a text file in MS/LOAD/GLOBAL/global_particle_load

**ndump_node_load** specifies the frequency at which to save node load
information i.e. the number of particles per node for all nodes. Output
is saved in hdf5 files using a grid with the same dimensions as the
parallel partition. Files are stored in MS/LOAD/NODE.

**ndump_grid_load** specifies the frequency at which to save the number
of particles per cell for all simulation cells. Output is saved in hdf5
files using a grid with the same dimensions as simulation grid. Files
are stored in MS/LOAD/CELL.

Here is an example of a grid section for a 2D run, using a 4096 × 512
cell grid and cartesian coordinates.

```text
grid 
{
  nx_p(1:2) = 4096, 512,
  coordinates = "cartesian",
}
```

Here is another example of a run using dynamic load balancing along x1:

```text
grid 
{
  nx_p(1:2) =  1024, 128,

  ! Load balance parameters
  load_balance(1:2) = .true., .false.,
  lb_type = "dynamic",
  n_dynamic = 32,
  cell_weight = 1.2, 
}
```
