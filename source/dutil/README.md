# Diagnostic I/O subsytem documentation

The new I/O subsystem for diagnostics replaces the previous HDF5 based module. The two main features are:

- Parallel I/O algorithm selection at runtime. 4 options are available, MPI, independent, MPI/IO independent and MPI/IO collective.
- Support for multiple file types: HDF5 and ZDF. HDF5 is no longer a pre-requisite for compiling OSIRIS

Additionally, a two-level I/O algorithm for parallel I/O that merges groups of grid cells prior to file output, greatly improving performance, has been implemented. The details are given in the grid and vdf sections.

### Parallel I/O algorithm

The parallel I/O algorithm can now be selected at runtime using the `parallel_io` parameter in the `simulation` section, e.g.:

```
simulation{
    parallel_io = "mpiio-coll",
}
```

Currently implemented options are:

- "_mpi_" - **[default]** Have root node write all the data. Data is sent by all nodes to root node sequentially, overlapping communication and I/O. This is the same algorithm currently used in OSIRIS. 
- "_indep_" - Have each node independently write its data to the file. This requires that all nodes have access to a shared file system. (not supported by HDF5 files)
- "_mpiio-indep_" - Use MPI/IO to have each node independently write its data to the file.
- "_mpiio-coll_" - Use MPI/IO to collectively write the data to the file. Depending on the file size, number of nodes and the (parallel) filesystem, MPI/IO will automatically choose the best algorithm for output. This is optimized for large files.

### Support for multiple-file types

OSIRIS can now do diagnostic output in multiple file types selectable at runtime. The default is to use the ZDF format, that was developed for the [ZPIC project](https://github.com/zambzamb/zpic). **visXD** has been updated to support these files, and there are also Python routines available for processing ZDF files.

HDF5, as always, is fully supported.. but no longer a required external dependency so that OSIRIS can be easily used in situations where HDF5 is not available. The Makefile will check the configuration file for the `H5_FCOMPILEFLAGS` variable; if it is defined then it assumes that HDF5 support is required, and will compile the necessary files, also making HDF5 the default option. Additionally, there is a new configuration file option, `H5_HAVE_PARALLEL`, that should be set to an arbitrary value if the HDF5 library was compiled with parallel I/O support. The `PARALLEL_IO` configuration file option is now obsolete and will be ignored.

The choice of file format can be done at runtime using the `file_format` parameter in the `simulation` section, e.g.:

```
simulation{
    file_format = "zdf",
}
```

Currently implemented options are:

- "_zdf_" - Use the ZDF format. This option is always available and is the default if HDF5 support was not selected at compile time. 
- "_hdf5_" - Use the HDF5 format. This option is only available if HDF5 support was selected at compile time, and if so is the default.

Additional file types (e.g. PnetCDF) will be added in the future.

## Notes for developers

The new new I/O subsystem for diagnostics is based on a `t_diag_file` class that implements the following interface methods:

- `init` - Initializes the `diag_file` object with data from `g_space`, `grid` and `no_co` objects. These are just used to get information regarding the global simulation data. 
- `open` - Open the diagnostic file for writing
- `close` - Close the diagnostic file
- `cleanup` - Free any temporary data allocated by the `diag_file` object
- `add_dataset` - Adds a dataset to the diagnostics file. Data may be from a single node or be distributed along multiple parallel nodes.

The `diag_file` class is organized using objects coming from the fortran interfaces defined in `dutil/os-diagfile.f03` and the use of a `grid_io` object to describe parallel configuration options. This allows us to greatly simplify the high-level interface.

### Adding a single node grid dataset

You begin by specifying the grid description in the `diag_file` object. In this example we will be saving a 2D [1000,512] dataset:

```fortran
class( t_diag_file ), allocatable :: diagFile

call create_diag_file( diagFile )

! File will be of type GRID
diagFile%ftype = p_diag_grid

! Iteration info
diagFile%iter%n         = 1000
diagFile%iter%t         = 100.0
diagFile%iter%time_units = "1 / \omega_p"

! Initialize grid file metadata
diagFile%grid%ndims = 2
diagFile%grid%name = 'charge'	! grid name - should be a proper computer variable name
diagFile%grid%label = '\rho'	! grid label - may include formatting commands (e.g. LaTeX)
diagFile%grid%units = 'e \omega_p^2 / c^2' ! grid units - may include formatting commands (e.g. LaTeX)

diagFile%grid%count(1) = 1024
diagFile%grid%count(2) = 512

diagFile%grid%axis(1)%min = 0.0
diagFile%grid%axis(1)%max = 1.0

diagFile%grid%axis(1)%name  = 'x1'           ! axis name - should be a proper computer variable name
diagFile%grid%axis(1)%label = 'x_1'          ! axis label - may include formatting commands (e.g. LaTeX)
diagFile%grid%axis(1)%units = 'c / \omega_p' ! axis units - may include formatting commands (e.g. LaTeX)

! (...) similar code for axis(2)
```

Once this has been done the dataset can written simply by calling `add_dataset`. The code will detect the data type automatically.

```fortran
! Notice that there is no file extension. This will be set from the file format.
diagFile % filename = 'test'
diagFile % filepath = '.'

! The code will default to the file format selected in the input file but you can
! also set it in the code:
! diagFile % fformat = p_hdf5_format

call diagFile % open( p_diag_create )
call diagFile % add_dataset( 'charge', buffer )
call diagFile % close( )
```

### Adding a chunked dataset

A chunked dataset is a dataset whose data is not saved all at once, but instead is split over multiple chunks. To add a chunked dataset you must follow three steps:

- Creating a chunked dataset (`start_cdset`). In this step you must define the name, global dimensions and datatype of the complete dataset.
- Writing individual chunks (`write_cdset`). A `t_diag_chunk` object is used to describe each individual chunk.
- Closing the dataset (`end_dataset`). This step is optional; if the dataset is not closed the reader will continue to search for chunks until end of file is reached.

Chunked datasets can be interleaved, i.e., you can start multiple datasets, write the chunks in any order, and close the datasets in any order.

The above example becomes:

```fortran
! diagFile is initialized previously
type( t_diag_datset ) :: dset
type( t_diag_chunk ) :: chunk

! Create chunked dataset
call diagFile%start_cdset( "data", 2, [1024,512], p_diag_float32, dset )

! Add first chunk
chunk % start(1:2) = [0,0]      ! start positions are 0 indexed
chunk % count(1:2) = [512,512]
chunk % stride(1:2) = [1,1]
chunk % data = c_loc( buffer1 )
call diagFile % write_cdset( dset, chunk )

! Add second chunk
chunk % start(1:2) = [512,0]
chunk % data = c_loc( buffer2 )
call diagFile % write_cdset( dset, chunk )

! Close dataset
call diagFile % end_cdset( dset )
```

### Extensible datasets

The module also allows for "extensible" datasets, i.e., datasets that are created with a given initial size, but that are later extended to a larger size. 

- `start_ext_cdset( name, dims, count, data_type, chunk_size, dset )`. The `chunk_size` parameter is only used by HDF5 format, and it defines how the dataset is extended.

#### Step 1: Create the extensible dataset

```fortran
! Open the file
call diagFile % open()

! Create extensible dataset
call diagFile % start_cdset( "ext_data", 2, [128,128], p_diag_float64, [64,64], dset )

! (optional) add some data before closing the file
chunk % count(1:2) = [64,64]
chunk % start(1:2) = [0,0]
chunk % stride(1:2) = [1,1]
chunk % data = c_loc( buffer )
call diagFile % write_cdset( dset, chunk )

! Close the file
call diagFile % close()
```

#### Step 2: Extend the dataset & write the data

```fortran
call diagFile % open( p_diag_update )

call diagFile % open_cdset( dset )

call diagFile % extend_cdset( dset, [256, 256] )

chunk % count(1:2) = [64,64]
chunk % start(1:2) = [128,128]
chunk % stride(1:2) = [1,1]
chunk % data = c_loc( buffer )
call diagFile % write_cdset( dset, chunk )

! If you open the dataset you should close it as some file formats require this
call diagFile % close_cdset( dset )

call diagFile % close()
```

# Parallel output

To save data distributed over a parallel MPI universe we begin by creating the file object and populating it with the usual data. We then create the file passing an additional parameter `comm` that specifies the MPI Communicator to use.

```fortran
class( t_diag_file ), allocatable :: diagFile

call create_diag_file( diagFile )

!(...)

diagFile % filename = 'test-parallel'
diagFile % filepath = '.'
diagFile % fformat = p_diag_format

! Notice the extra parameter 'comm' specifying the MPI communicator to use.
call diagFile % open( p_diag_create, comm )
```

When opening the file you may optionally set the Parallel I/O mode to use, by adding an additional parameter to the `diagFile % open()` call:

```fortran
! Notice the extra parameter 'DIAG_MPI' specifying the parallel I/O algorithm to use.
call diagFile % open( p_diag_create, comm, DIAG_MPI )
```

 Valid values for this parameter are:

- `DIAG_MPI` - Gather data on root node and write it. Communication and IO are overlapped i.e. while receiving the data from one node, the code is writing the data from the previous node.
- `DIAG_INDEPENDENT`- Have each node write its data independently. Requires that all nodes have access to the filesystem.
- `DIAG_MPIIO_INDEPENDENT` - Write data using MPI/IO independent file access.
- `DIAG_MPIIO_COLLECTIVE` - Write data using MPI/IO collective file access.

If this parameter is not set, the code will default to the value specified in the input file.

The next step is to create a chunked dataset that describes the global data. This is done with exactly the same interface:

```fortran
type( t_diag_dataset ) :: parDset

! Create chunked dataset describing global data
call diagFile%start_cdset( "data", 2, [1024,512], p_diag_float32, dset )
```

To save the parallel data we use a `t_diag_chunk`object, as we did earlier, to describe the data on the local node. 

```fortran
type( t_diag_chunk ) :: parChunk

! chunk contains local data
parChunk % local % count(1:2) = [16,16]
parChunk % local % start(1:2) = [128,0]
parChunk % local % stride(1:2) = 1
parChunk % local % data = c_loc( data )
```

Finally, write the data using the `write_par_cdset()` method. The module will merge the data from the multiple nodes using the parallel I/O algorithm defined earlier.

```fortran
! Write the data
call diagFile % write_par_cdset( parDset, parChunk )
```

The dataset and the file can then be closed in the usual way:

```fortran
! Close dataset
call diagFile % end_cdset( dset )

! Close the file
call diagFile % close()
```



