# IOPerf Tool

The IOPerf tool is used for measuring write performance of the OSIRIS I/O subsystem for parallel grids (VDFs). It uses exactly the same routines as osiris to create a parallel grid and save it to file, measuring the bandwidth of the process.

### Configuration and compilation

Since the tool uses the same routines from OSIRIS, the dimensionality (1D, 2D or 3D) is defined at compile time, using the same configuration file used for OSIRIS. The precision for the diagnostics file is set to double. To compile the tool, you follow the same process as if you were compiling OSIRIS, including the configuration step, but you specify `ioperf`as the target for `make`:

```shell
$ make ioperf
```

This will compile the IOPerf binary and place it on the usual directory.

### Running IOPerf

The IOPerf tool takes no command line parameters, and reads test run-time parameters from a file named `ioperf.params` that must be present in the directory where the tool is running. This file uses a similar format to the OSIRIS input file, but using only 2 sections:

* `node_conf` - Parallel node configuration
* `test_set` - IOPerf test parameters, described in detail below

The first section is exactly the same as in any OSIRIS simulation. The `periodic`and `n_threads` parameters will be ignored.

To run the tool, you launch it just as you do for OSIRIS, but using the `ioperf-*` binary, e.g.:

```shell
$ mpirun -np 16 ../bin/ioperf-2D.e
```

The IOPerf tool will output all the results to stdout and summary information to a text file name `ioperf.out`

##IOPerf Options

The current version of IOPerf is designed to do parameter scans. To do this, the user may specify up to 16 different values for each parameter, and IOPerf will test each possible combinations.

### Parameter scan options

* `str(*) :: grid` - Grid sizes to test specified as comma separated values e.g. `grid = "128,128,128", "256,256,256"`. The user must supply at least one value for this parameter, there is no default value.
* `str(*) :: format` - File format to use for the tests, must be one of `"zdf"` or `"hdf5"`. Default is `"zdf"` (tests `zdf` file format only)
* `str(*) :: parallel_io ` - Parallel I/O algorithm to use, must be one of "mpi" (use merge data on node 0, overlapping communication and I/O), "independent" (use MPI/IO independent writes), or "collective" (Use MPI/IO collective writes). Default is `"mpi"`
* `str(*) :: n_merge`- Number of parallel nodes in each direction that should merge data before doing I/O specified as comma separated values e.g. `n_merge = "2,2,1","4,1,1",` 
* `str(*) :: merge_type` - Algorithm used for merging data before doing I/O, must be one of `"none"` (no merging), `"point2point"` (use point to point communications) or `"gather"` (use gather operations). Default is `"none"`.

### Global test options

* `str :: path` - Path in which to save the test files. Defaults to `""`(empty string)
* `integer :: repeat` - Number of times that each individual test is repeated to improve measurement. Defaults to 3. 
* `logical :: keepFiles` - Keeps files at the end of each test. Defaults to `.false.`



## Example

The file below saves a 3D grid from a 4096 ($16 ^3$) node partition. The code will test all possible combinations of the parameters specified leading to a total of 48 different tests. Each test will be repeated 11 times to improve statistics, and the files will be removed after testing.

```
! This file should be named ioperf.params

node_conf {  
	node_number(1:3) =  16, 16, 16
}


test_set {

  ! Grid dimensions 128^3 and 256^3
  grid = "128,128,128", "256,256,256",
  
  ! Test both hdf5 and zdf formats
  format = "hdf5","zdf",
  
  ! Test all parallel io algorithms
  parallel_io = "mpi", "independent", "collective",
  
  ! Test all merging algorithms
  merge_type = "point2point", "gather",
  
  ! no merging and merge every 2 nodes along x1
  n_merge = "1,1,1","2,1,1",

	! Repeats each test 11 times
	repeat = 11,

}
```

