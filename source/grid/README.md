#Grid object documentation

The grid object will be moved from a Fortran 90 module structure to a full Fortran 2003 class in the near future, this documents the currently available methods.

**module `m_grid` functions**


* `setup` - Initialize data structures following reading input parameters
* `restart_write` - write checkpoint information
* `restart_read` - read checkpoint information. Gets called from setup if restarting from a previous simulation.
* `cleanup` - Cleanup data structure

* `g_nx` - Gets the global simulation grid size
* `num_partitions` - Gets the number of partitions in each direction. You can also 

* `my_nx_p` - Gets position of lower/upper boundary of local grid node on the global simulation grid. For convenience it also returns the size of the local grid node.
* `my_nx_p_min` - Gets the position of the lower boundary of the local grid node on the global simulation grid.
* `nx_p` - Gets position of lower/upper boundary of the specified node on the global simulation grid. For convenience it also returns the size of the specified grid node.
* `get_node_limits` - Gets position of lower/upper boundary of all nodes on the global simulation grid. For convenience it also returns the size of all the grid nodes.

* `nx_p_part` - Gets position of lower/upper boundary of a given partition / direction on the global simulation grid. For convenience it also returns the size of the specified partition.

* `my_nx` - Gets the size (number of cells in each direction) of the local grid node
get this information from the node_conf object.
* `get_part_width` - Gets all grid node sizes along the specified direction
* `get_max_nx` - Get maximum number of cells in a given direction on any grid node
* `local_vol` - Gets the total size (product of the number of cells in each direction) of the local grid node


* `copy` - Creates a copy of the object
* `equal_my_nx` - Checks if the position of the local grid node on the global simulation grid is the same for 2 grid objects
* `operator(==)` - Checks if two grid objects are identical
* `operator(/=)` - Checks if two grid objects are not identical
* `x_dim` - Returns the number of dimensions of the the grid
* `coordinates` - Returns the type of coordinates being used (cartesian or r-z cylindrical)
* `test_partition` - Tests if the parallel partition is compatible with the grid parameters.


Functions defined in `os-grid.f90` outside of the `m_grid` module:

* `read_input_grid` - Reads input parameters for grid object


#Node conf object documentation


* `read_nml` - read configuration from input file
* `setup` - setup data structures following reading input parameters
* `restart_write` - write checkpoint information
* `restart_read` - read checkpoint information. Gets called from setup if restarting from a previous simulation.
* `cleanup` - Cleanup data structure

* `no_num` - Gets total number of nodes
* `nx` - Gets number of nodes in all directions / specific direction
* `neighbor` - Gets id of the neighboring parallel nodes
* `my_aid` - Get id of local node
* `ngp` - Returns the position of the requested node on the parallel partition
* `my_ngp` - Returns the position of the local node on the parallel partition
* `ifpr` - Returns true if single node periodic along the specified direction
* `ngp_id` - returns a unique id integer based on the node grid position
* `send_ping` - Sends a short message (ping) to another node (blocking). Used for synchronization between 2 processes.
* `recv_ping` - Receives a short message (ping) to another node (blocking). Used for synchronization between 2 processes.
* `reduce_array` - Array parallel reduce operation
* `reduce_array_size` - Array parallel reduce operation, using only the specified size
* `reduce` - Scalar parallel reduce operation
* `gather_array` - Array gather operation
* `allgather` - Scalar all gather operation
* `broadcast` - Scalar broadcast operation
* `wait_for_all` - Global barrier operation
* `periodic` - Returns true if parallel partition is periodic along specified direction
* `new_slice` - Creates a new node_conf object representing a slice of the old one.
* `comm` - Returns the MPI communicator used

* `on_edge` - Returns true if the local node is on the specified edge of the parallel partition
* `phys_bound` - Returns true if the specified / local node is a physical boundary of the simulation
* `n_threads` - Returns selected number of threads per node
* `root` - Returns true if this is the root node of the simulation



