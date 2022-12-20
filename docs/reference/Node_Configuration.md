---
layout: single
classes: wide
title: Node configuration
permalink: /reference/node_conf
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures the parallel nodes and periodic boundary
settings and must be present in the input file. It accepts the following
data:

- **node_number**(x_dim), integer, default = 1
- **n_threads**, integer, default = 1
- **if_periodic**(x_dim), boolean, default = .false.
- **topology**, character(\*), default = "mpi"

**node_number** specifies the number of nodes to use in each direction
for the simulation. The total number of nodes will be the product of the
number of nodes for each direction. A single node run will be specified
by setting all the items to 1. It is not required that the number of
grid points on a given direction is evenly dividable by the number of
nodes specified for that direction, but this will guaranty better load
balancing. Note that only longitudinal partitions are allowed when
algorithm = 'pgc' in the general simulation parameters section.

**n_threads** specifies the number of threads per node to use. This
option uses a shared memory algorithm inside each MPI node, that can
significantly improve load imbalance issues. To use this option the code
must be compiled with OpenMP support. Setting this value to 1 will use
the standard distributed memory algorithm only.

**if_periodic** specifies if the boundary conditions for each direction
will be periodic boundary conditions. If set, they will override any
boundary conditions you specify for Electro-Magnetic fields or particle
species.

**topology** specifies how osiris should map simulation nodes to
parallel nodes. Currently available options are:

- *"mpi"* - Use MPI functions (namely MPI_CART_CREATE) to generate the
  topology. MPI can use knowledge about network details to choose the
  best topology. This is the default.
- *"old"* - Use the old OSIRIS algorithm to generate the topology. This
  should only be used for testing purposes and comparing to old versions
  of the code.
- *"bgq"* - (available on BlueGene/Q systems only). This generates a
  topology based on the 5D torus topology of the BlueGene/Q system. This
  should provide the best performance on these systems.

Here's an example of a node_conf section for a 2D run, using 8 nodes in
the $x_1$ direction and 4 nodes in the $x_2$ direction,
and launching 2 threads per node. The run will use periodic boundaries
in the $x_2$ direction.

```text
node_conf 
{
  node_number(1:2) = 8, 4,
  n_threads = 2,
  if_periodic(1:2) = .false., .true.,
}
```
