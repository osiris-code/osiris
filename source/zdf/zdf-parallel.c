
#define _POSIX_C_SOURCE 200809L
#define _FILE_OFFSET_BITS 64

#if defined(_WIN32) || defined(_WIN64) 
  
  #include <io.h>
  #include <fcntl.h>
  #include <sys/stat.h>
  #include <windows.h>
  #include <stdio.h>
  
/**
 * Windows does not have a direct 'fmemopen' replacement
     (i.e. and easy to use mem-mapper fit for large files..
     so we simulateit with a temp-file.
 */  
  FILE *fmemopen(void *buf, size_t len, const char *type)
  {
    int fd;
    FILE *fp;
    char tp[MAX_PATH - 13];
    char fn[MAX_PATH + 1];

    if (!GetTempPathA(sizeof(tp), tp))
      return NULL;

    if (!GetTempFileNameA(tp, "confuse", 0, fn))
      return NULL;

    fd = _open(fn,
      _O_CREAT | _O_RDWR | _O_SHORT_LIVED | _O_TEMPORARY | _O_BINARY,
      _S_IREAD | _S_IWRITE);
    if (fd == -1)
      return NULL;

    fp = _fdopen(fd, "w+");
    if (!fp) {
      _close(fd);
      return NULL;
    }

    fwrite(buf, len, 1, fp);
    rewind(fp);

    return fp;
  }

  // replace 'fseeko' amd 'ftello' with their Windows cognates
  #define fseeko _fseeki64
  #define ftello _ftelli64
#endif

#include "zdf-parallel.h"
#include <mpi.h>
#include <sys/types.h>
#include <stdlib.h>

/**
 * Size of temporary buffer for writing metadata
 */
#define TMP_BUFFER_SIZE 1024

/**
 * Returns MPI data type from ZDF data type
 * @param  data_type ZDF data type
 * @return           Corresponding MPI data type. Returns MPI_DATATYPE_NULL for unsupported
 *                   datatypes.
 */
MPI_Datatype zdf_mpi_type( int data_type ) {

	switch (data_type) {
		case( zdf_null ):    return( MPI_DATATYPE_NULL );
		case( zdf_float32 ): return( MPI_FLOAT    );
		case( zdf_float64 ): return( MPI_DOUBLE   );
		case( zdf_int8 ):    return( MPI_INT8_T   );
		case( zdf_int16 ):   return( MPI_INT16_T  );
		case( zdf_int32 ):   return( MPI_INT32_T  );
		case( zdf_int64 ):   return( MPI_INT64_T  );
		case( zdf_uint8 ):   return( MPI_UINT8_T  );
		case( zdf_uint16 ):  return( MPI_UINT16_T );
		case( zdf_uint32 ):  return( MPI_UINT32_T );
		case( zdf_uint64 ):  return( MPI_UINT64_T );
	}
	return( MPI_DATATYPE_NULL );
}

/**
 * Opens file for MPI parallel I/O
 * @param  zdf		Parallel ZDF file handle
 * @param  filename Filename of the ZDF file to open, including path
 * @return          Returns 1 on success, 0 on error
 */
int zdf_par_open_file_mpi( t_zdf_par_file* zdf, const char* filename ) {

	// Only root node opens file
	int rank;
	MPI_Comm_rank( zdf->comm, &rank );
	if ( rank == 0 )
		if( !zdf_open_file( & zdf->zdf_file, filename, zdf -> zdf_file.mode ) ) return(0);

	return(1);
}

/**
 * Opens file for Independent parallel I/O
 * @param  zdf		Parallel ZDF file handle
 * @param  filename Filename of the ZDF file to open, including path
 * @return          Returns 1 on success, 0 on error
 */
int zdf_par_open_file_indep( t_zdf_par_file* zdf, const char* filename ) {

	int rank; MPI_Comm_rank( zdf->comm, &rank );

	if ( !rank ) {
		// Root opens file using serial routine
		if ( zdf_open_file( &zdf->zdf_file, filename, zdf->zdf_file.mode ) == 0 ) {
			return(0);
		}
		if ( zdf->zdf_file.mode == ZDF_CREATE )	MPI_Barrier( zdf -> comm );
	} else {
		// Other nodes open file for updating after it has been created by root
		if ( zdf->zdf_file.mode == ZDF_CREATE )	MPI_Barrier( zdf -> comm );

		if ( zdf->zdf_file.mode == ZDF_READ ) {
			if (!(zdf->zdf_file.fp = fopen( filename, "r"))) {
				perror("(*error*) Unable to open ZDF file for reading.\n");
				return(0);
			}
		} else {
			if (!(zdf->zdf_file.fp = fopen( filename, "r+"))) {
				perror("(*error*) Unable to open ZDF file for reading / writing.\n");
				return(0);
			}
		}
	}

	return(1);
}

/**
 * Opens file for MPI/IO parallel I/O
 * @param  zdf		Parallel ZDF file handle
 * @param  filename Filename of the ZDF file to open, including path
 * @return          Returns 1 on success, 0 on error
 */
int zdf_par_open_file_mpiio( t_zdf_par_file* zdf, const char* filename ) {

	int rank;
	MPI_Comm_rank( zdf->comm, &rank );

	// MPI_MODE_UNIQUE_OPEN mode allows an implementation to optimize access by 
	// eliminating the overhead of file locking.
	int mpiio_amode = MPI_MODE_UNIQUE_OPEN;
	switch ( zdf -> zdf_file.mode ) {
		case( ZDF_CREATE ):
			mpiio_amode |= MPI_MODE_CREATE | MPI_MODE_WRONLY;
			break;
		case( ZDF_READ ):
			mpiio_amode |= MPI_MODE_RDONLY;
			break;
		case( ZDF_UPDATE ):
			mpiio_amode |= MPI_MODE_RDWR;
	}

	MPI_Info hints = MPI_INFO_NULL;
	MPI_Info_create(&hints);

	// Enable collective buffering optimization
	MPI_Info_set(hints , "romio_cb_write" , "enable" );	

	if ( MPI_File_open( zdf->comm, filename, mpiio_amode,  hints,
										&(zdf -> fh) ) != MPI_SUCCESS ) {
		fprintf(stderr, "(*error*) Unable to open ZDF file %s\n", filename );
		return(0);
	}

	// Free info object
	MPI_Info_free( &hints );

	// Truncate existing file upon creation
	if ( zdf -> zdf_file.mode == ZDF_CREATE ) {
		if ( MPI_File_set_size( zdf -> fh, 0 ) != MPI_SUCCESS ) return(0);
		// According to the ROMIO documentation we should perform a barrier following a truncation
		MPI_Barrier( zdf->comm );
	}

	if ( rank == 0 ) {
		if ( zdf -> zdf_file.mode == ZDF_CREATE ) {
			// Write magic number
			if ( MPI_File_write( zdf->fh, (void*) zdf_magic, ZDF_MAGIC_LENGTH, MPI_BYTE,
				MPI_STATUS_IGNORE) != MPI_SUCCESS ) {
				fprintf(stderr, "(*error*) Unable to write magic number.\n" );
				zdf_par_close_file( zdf );
				return(0);
			}
		} else {
			char test_magic[ZDF_MAGIC_LENGTH];
			// Read magic number
			if ( MPI_File_read( zdf -> fh, (void *) test_magic, ZDF_MAGIC_LENGTH, MPI_BYTE,
				MPI_STATUS_IGNORE ) != MPI_SUCCESS ) {
				fprintf(stderr, "(*error*) Unable to read magic number from ZDF file.\n");
				zdf_par_close_file( zdf );
				return(0);
			}

			// Check magic number
			for( int i = 0; i < ZDF_MAGIC_LENGTH; i++) {
				if ( test_magic[i] != zdf_magic[i] ) {
					fprintf(stderr, "(*error*) Invalid magic number, file is not a proper ZDF file.\n");
					zdf_par_close_file( zdf );
					return(0);
				}
			}
		}
	}

	return(1);
}

/**
 * Opens Parallel ZDF file
 * @param  zdf      Parallel ZDF file handle
 * @param  filename Filename of the ZDF file to open, including path
 * @param  mode     Can be one of ZDF_WRITE (for writing) or ZDF_READ (for reading)
 * @param  comm     MPI Communicator to use
 * @return          Returns 1 on success, 0 on error
 */
int zdf_par_open_file( t_zdf_par_file* zdf, const char* filename, enum zdf_file_access_mode amode,
	MPI_Comm comm, enum zdf_parallel_io_mode iomode ) {

	zdf -> zdf_file.mode = amode;
	zdf -> zdf_file.ndatasets = 0;
	zdf -> comm = comm;
	zdf -> iomode = iomode;

	zdf -> zdf_file.fp = NULL;
	zdf -> fh = MPI_FILE_NULL;

	switch ( zdf -> iomode ) {
		case ZDF_MPI:
			// Open file for MPI access
			if( !zdf_par_open_file_mpi( zdf, filename ) ) return(0);
			break;
		case ZDF_INDEPENDENT:
			// Open file for independent access
			if( !zdf_par_open_file_indep( zdf, filename ) ) return(0);
			break;
		case ZDF_MPIIO_INDEPENDENT:
		case ZDF_MPIIO_COLLECTIVE:
			// Open file for MPI/IO access
			if( !zdf_par_open_file_mpiio( zdf, filename ) ) return(0);
			break;
	}

	// Put file pointer after magic number
	zdf -> fpos = ZDF_MAGIC_LENGTH;

	return(1);
}


/**
 * Closes parallel ZDF file
 * @param  zdf Parallel ZDF file handle
 * @return     Returns 1 on success, 0 on error
 */
int zdf_par_close_file( t_zdf_par_file* zdf ) {

	switch ( zdf -> iomode ) {
		case ZDF_MPI:{
			int rank; MPI_Comm_rank( zdf->comm, &rank );
				if ( rank == 0 ) {
					if ( !zdf_close_file( &zdf -> zdf_file) ) return(0);
				}
			}
			break;
		case ZDF_INDEPENDENT:
			if ( !zdf_close_file( &zdf -> zdf_file) ) return(0);
			break;
		case ZDF_MPIIO_INDEPENDENT:
		case ZDF_MPIIO_COLLECTIVE: 
			if ( MPI_File_close( & (zdf -> fh) ) != MPI_SUCCESS ) return(0);
	}

	return(1);
}

/**
 * Adds string element to parallel ZDF file. Only root node writes data.
 * @param  zdf  Parallel ZDF file handle
 * @param  name Element name
 * @param  str  String value
 * @return      Returns 1 on success, 0 on error
 */
int zdf_par_add_string( t_zdf_par_file* zdf, const char* name, const char* str ) {

	size_t count;
	int rank; MPI_Comm_rank( zdf->comm, &rank );

	switch ( zdf -> iomode ) {
		case ZDF_MPI:
		case ZDF_INDEPENDENT:
			if ( !rank ) {
				if ( !(count =	zdf_add_string( &zdf->zdf_file, name, str ))) return(0);
				zdf -> fpos += count;
			}
			break;
		case ZDF_MPIIO_INDEPENDENT:
		case ZDF_MPIIO_COLLECTIVE: {
			uint8_t buf[TMP_BUFFER_SIZE];
			t_zdf_file tmp = { .fp = fmemopen( (void *) buf, TMP_BUFFER_SIZE, "w") };
			if ( !(count =	zdf_add_string( &tmp, name, str ))) return(0);
			fclose(tmp.fp);
			if ( !rank ) MPI_File_write( zdf->fh, (void*) buf, count, MPI_BYTE, MPI_STATUS_IGNORE);
			zdf -> fpos += count;
		}
	}

	return(1);
}

/**
 * Adds int32 element to parallel ZDF file. Only root node writes data.
 * @param  zdf    Parallel ZDF file handle
 * @param  name   Element name
 * @param  value  int32 value
 * @return      Returns 0 on error, other value on success.
 */
int zdf_par_add_int32( t_zdf_par_file* zdf, const char* name, const int32_t value ) {

	size_t count;
	int rank; MPI_Comm_rank( zdf->comm, &rank );

	switch ( zdf -> iomode ) {
		case ZDF_MPI:
		case ZDF_INDEPENDENT:
			if ( !rank ) {
				if ( !(count =	zdf_add_int32( &zdf->zdf_file, name, value ))) return(0);
				zdf -> fpos += count;
			}
			break;
		case ZDF_MPIIO_INDEPENDENT:
		case ZDF_MPIIO_COLLECTIVE: {
			uint8_t buf[TMP_BUFFER_SIZE];
			t_zdf_file tmp = { .fp = fmemopen( (void *) buf, TMP_BUFFER_SIZE, "w") };
			if (!(count =	zdf_add_int32( &tmp, name, value ))) return(0);
			fclose(tmp.fp);
			if ( !rank ) MPI_File_write( zdf->fh, (void*) buf, count, MPI_BYTE, MPI_STATUS_IGNORE);
			zdf -> fpos += count;
		}
	}

	return(1);
}

/**
 * Adds float64 element to parallel ZDF file. Only root node writes data.
 * @param  zdf    Parallel ZDF file handle
 * @param  name   Element name
 * @param  value  float64 value
  * @return      Returns 0 on error, other value on success.
*/
int zdf_par_add_double( t_zdf_par_file* zdf, const char* name, const double value ) {

	size_t count;
	int rank; MPI_Comm_rank( zdf->comm, &rank );

	switch ( zdf -> iomode ) {
		case ZDF_MPI:
		case ZDF_INDEPENDENT:
			if ( !rank ) {
				if ( !(count =	zdf_add_double( &zdf->zdf_file, name, value ))) return(0);
				zdf -> fpos += count;
			}
			break;
		case ZDF_MPIIO_INDEPENDENT:
		case ZDF_MPIIO_COLLECTIVE: {
			uint8_t buf[TMP_BUFFER_SIZE];
			t_zdf_file tmp = { .fp = fmemopen( (void *) buf, TMP_BUFFER_SIZE, "w") };
			if (!(count = zdf_add_double( &tmp, name, value ))) return(0);
			fclose(tmp.fp);
			if ( !rank ) MPI_File_write( zdf->fh, (void*) buf, count, MPI_BYTE, MPI_STATUS_IGNORE);
			zdf -> fpos += count;
		}
	}

	return(1);
}

/**
 * Adds iteration metadata group to parallel ZDF file. Only root node writes data.
 * @param  zdf  Parallel ZDF file handle
 * @param  iter Iteration info
 * @return      Returns 0 on error, other value on success.
 */
int zdf_par_add_iteration( t_zdf_par_file* zdf, const t_zdf_iteration* iter ) {

	size_t count;
	int rank; MPI_Comm_rank( zdf->comm, &rank );

	switch ( zdf -> iomode ) {
		case ZDF_MPI:
		case ZDF_INDEPENDENT:
			if ( !rank ) {
				if ( !(count = zdf_add_iteration( &zdf -> zdf_file, iter ))) return(0);
				zdf -> fpos += count;
			}
			break;
		case ZDF_MPIIO_INDEPENDENT:
		case ZDF_MPIIO_COLLECTIVE: {
			uint8_t buf[TMP_BUFFER_SIZE];
			t_zdf_file tmp = { .fp = fmemopen( (void *) buf, TMP_BUFFER_SIZE, "w") };
			if (!(count =	zdf_add_iteration( &tmp, iter ))) return(0);
			fclose(tmp.fp);
			if ( !rank ) MPI_File_write( zdf->fh, (void*) buf, count, MPI_BYTE, MPI_STATUS_IGNORE);
			zdf -> fpos += count;
		}
	}

	return(1);
}

/**
 * Adds parallel grid information metadata group to file. Only root node writes data.
 * @param  zdf  Parallel ZDF file handle
 * @param  grid Grid information
 * @return      Returns 0 on error, other value on success.
 */
int zdf_par_add_grid_info( t_zdf_par_file* zdf, const t_zdf_grid_info* grid ) {

	size_t count;
	int rank; MPI_Comm_rank( zdf->comm, &rank );

	switch ( zdf -> iomode ) {
		case ZDF_MPI:
		case ZDF_INDEPENDENT:
			if ( !rank ) {
				if ( !(count = zdf_add_grid_info( &zdf->zdf_file, grid ))) return(0);
				zdf -> fpos += count;
			}
			break;
		case ZDF_MPIIO_INDEPENDENT:
		case ZDF_MPIIO_COLLECTIVE: {
			uint8_t buf[TMP_BUFFER_SIZE];
			t_zdf_file tmp = { .fp = fmemopen( (void *) buf, TMP_BUFFER_SIZE, "w") };
			if (!(count = zdf_add_grid_info( &tmp, grid ))) return(0);
			fclose(tmp.fp);
			if ( !rank ) MPI_File_write( zdf->fh, (void*) buf, count, MPI_BYTE, MPI_STATUS_IGNORE);
			zdf -> fpos += count;
		}
	}

	return(1);
}

/**
 * Adds particle information metadata group to parallel ZDF file. Only root node writes data.
 * @param  zdf  Parallel ZDF file handle
 * @param  part Particle information
 * @return      Returns 0 on error, other value on success.
*/
int zdf_par_add_part_info( t_zdf_par_file* zdf, t_zdf_part_info* part ) {

	size_t count;
	int rank; MPI_Comm_rank( zdf->comm, &rank );

	switch ( zdf -> iomode ) {
		case ZDF_MPI:
		case ZDF_INDEPENDENT:
			if ( !rank ) {
				if ( !(count = zdf_add_part_info( &(zdf->zdf_file), part ))) return(0);
				zdf -> fpos += count;
			}
			break;
		case ZDF_MPIIO_INDEPENDENT:
		case ZDF_MPIIO_COLLECTIVE: {
			uint8_t buf[TMP_BUFFER_SIZE];
			t_zdf_file tmp = { .fp = fmemopen( (void *) buf, TMP_BUFFER_SIZE, "w") };
			if (!(count =	zdf_add_part_info( &tmp, part ))) return(0);
			fclose(tmp.fp);
			if ( !rank ) MPI_File_write( zdf->fh, (void*) buf, count, MPI_BYTE, MPI_STATUS_IGNORE);
			zdf -> fpos += count;
		}
	}

	return(1);
}

/**
 * Adds a chunked dataset header to parallel ZDF file. Only root node writes data.
 * @param  zdf     Parallel ZDF file handle
 * @param  dataset Dataset object
 * @return         0 on error, other value on success. The dataset object (on all nodes)
 *                 is also modified to include a unique id.
 */
int zdf_par_start_cdset( t_zdf_par_file* zdf, t_zdf_dataset* dataset ) {

	int rank; MPI_Comm_rank( zdf->comm, &rank );
	size_t count;

	switch ( zdf -> iomode ) {
		case ZDF_MPI:
		case ZDF_INDEPENDENT:
			if ( !rank ) {
				if (!(count = zdf_start_cdset( &(zdf->zdf_file), dataset ))) return(0);
				zdf -> fpos += count;
			} else zdf -> zdf_file.ndatasets++;
			dataset -> id = zdf -> zdf_file.ndatasets;
			break;
		case ZDF_MPIIO_INDEPENDENT:
		case ZDF_MPIIO_COLLECTIVE: {
			uint8_t buf[TMP_BUFFER_SIZE];
			t_zdf_file tmp = {
				.fp = fmemopen( (void *) buf, TMP_BUFFER_SIZE, "w"),
				.ndatasets = zdf -> zdf_file.ndatasets
			};
			if (!(count = zdf_start_cdset( &tmp, dataset ))) return(0);
			fclose(tmp.fp);
			if ( !rank ) MPI_File_write( zdf->fh, (void*) buf, count, MPI_BYTE, MPI_STATUS_IGNORE);
			zdf -> fpos += count;
			dataset -> id = ++zdf -> zdf_file.ndatasets;
		}
	}
	
	return(1);
}

/**
 * Extends dataset dimensions. After this call the file pointer is positioned at the
 * end of the file.
 *
 * @param  zdf       Parallel ZDF file handle
 * @param  dataset   Datasete object
 * @param  new_count New dimensions for the dataset, must be >= than the previous
 *                   dimensions.
 * @return           Returns 1 if successful, 0 otherwise.
 */
int zdf_par_extend_dataset( t_zdf_par_file* zdf, t_zdf_dataset* dataset, uint64_t* new_count ) {

    int rank, ierr; 

    MPI_Comm_rank( zdf->comm, &rank );
    if ( !rank ) {
        ierr = zdf_extend_dataset( &(zdf -> zdf_file), dataset, new_count );
    } else {
        for( unsigned i = 0; i < dataset -> ndims; i++ ) {
            if ( new_count[i] < dataset -> count[i] ) {
                fprintf(stderr,"(*error*) Invalid value for zdf_extend_dataset.\n");
                fprintf(stderr,"(*error*) New size is smaller than original size.\n");
                return(-1);
            }
            dataset -> count[i] = new_count[i];
        }
        ierr = 0;
    }
	return(ierr);
}

/**
 * Add an end marker for a chunked dataset to a parallel ZDF file. Only root node writes data.
 * @param  zdf     Parallel ZDF file handle
 * @param  dataset Dataset object
 * @return         1 on success, 0 on error
 */
int zdf_par_end_cdset( t_zdf_par_file* zdf, t_zdf_dataset* dataset ) {

	size_t count;
	int rank; MPI_Comm_rank( zdf->comm, &rank );

	switch ( zdf -> iomode ) {
		case ZDF_MPI:
		case ZDF_INDEPENDENT:
			if ( !rank ) {
				if (!(count = zdf_end_cdset( &(zdf->zdf_file), dataset ))) return(0);
				zdf -> fpos += count;
			}
			break;
		case ZDF_MPIIO_INDEPENDENT:
		case ZDF_MPIIO_COLLECTIVE: {
			uint8_t buf[TMP_BUFFER_SIZE];
			t_zdf_file tmp = { .fp = fmemopen( (void *) buf, TMP_BUFFER_SIZE, "w") };
			if (!(count = zdf_end_cdset( &tmp, dataset ))) return(0);
			fclose(tmp.fp);
			if ( !rank ) MPI_File_write( zdf->fh, (void*) buf, count, MPI_BYTE, MPI_STATUS_IGNORE);
			zdf -> fpos += count;
		}
	}

	return(1);
}


/**
 * Write data that is distributed over a parallel universe using MPI parallel I/O. Only root node
 * writes to disk, overlapping communication and file I/O.
 * @param  zdf     Parallel ZDF file handle
 * @param  dataset Dataset to write. Dataset parameters describe global data
 * @param  chunk   Local dataset chunk (includes pointer to data)
 * @return         Returns 1 if successful, 0 otherwise.
 */
size_t zdf_par_write_cdset_mpi( t_zdf_par_file* zdf, t_zdf_dataset* dataset, t_zdf_chunk* chunk ) {

	enum commm_tags { ping_tag = 1001, tile_tag, data_tag };
	int ping = 1234;

	const int ndims = dataset -> ndims;
	const size_t type_size = zdf_sizeof( dataset -> data_type );
	const MPI_Datatype mpi_type = zdf_mpi_type( dataset -> data_type );

	const MPI_Comm comm = zdf -> comm;
	int comm_rank;
	MPI_Comm_rank( comm, &comm_rank );

	if ( comm_rank == 0 ) {

		uint64_t size0 = 0, size1 = 0;

		t_zdf_chunk buffer0 = {.data = chunk -> data};
		for(int i = 0; i < ndims; i ++ ){
			buffer0.count[i]  = chunk -> count[i];
			buffer0.start[i]  = chunk -> start[i];
			buffer0.stride[i] = chunk -> stride[i];
		}
		t_zdf_chunk buffer1;

		int comm_size; MPI_Comm_size( comm, &comm_size );

		for( int src = 0, step = 0; src < comm_size; src++, step = !step ) {
			MPI_Request ping_req, data_req;

			// Start receiving data from next node
			if ( src + 1 < comm_size ) {

				MPI_Isend( (void *) &ping, 1, MPI_INT, src + 1, ping_tag, comm, &ping_req );

				t_zdf_chunk* buffer = (step)? &buffer0 : &buffer1;
				uint64_t*    size   = (step)? &size0 : &size1;

				MPI_Recv( (void *) buffer, 9, MPI_UINT64_T, src + 1, tile_tag, comm, MPI_STATUS_IGNORE);
				uint64_t recv_size = 1;
				for(int i = 0; i < ndims; i++)
					recv_size *= buffer -> count[i];

				if ( recv_size > *size ) {
					if ( *size > 0 ) free( buffer -> data );
					buffer -> data = malloc( recv_size * type_size );
					*size = recv_size;
				}

				MPI_Irecv( buffer->data, recv_size, mpi_type, src + 1, data_tag, comm, &data_req );
			}

			// Write data from current node
			t_zdf_chunk* buffer = (step)? &buffer1 : &buffer0;
			zdf_write_cdset( &zdf->zdf_file, dataset, buffer );

			// Wait for messages to complete
			if ( src + 1 < comm_size ) {
				MPI_Wait( &ping_req, MPI_STATUS_IGNORE );
				MPI_Wait( &data_req, MPI_STATUS_IGNORE );
			}

		}

		// Free communication buffer
		if ( size0 > 0 ) free( buffer0.data );
		if ( size1 > 0 ) free( buffer1.data );

	} else {
		// Wait for ping from root node
		MPI_Recv( (void *) &ping, 1, MPI_INT, 0, ping_tag, comm, MPI_STATUS_IGNORE );

		// Send local chunk information to root node
		int buffer_size = 1;
		t_zdf_chunk buffer = { .data = dataset -> data };
		for(int i = 0; i < ndims; i ++ ){
			buffer.count[i]  = chunk -> count[i];
			buffer.start[i]  = chunk -> start[i];
			buffer.stride[i] = chunk -> stride[i];
			buffer_size     *= chunk -> count[i];
		}
		MPI_Send( (void *) &buffer, 9, MPI_UINT64_T, 0, tile_tag, comm );

		// Send data to root node
		MPI_Send( chunk -> data, buffer_size, mpi_type, 0, data_tag, comm );

	}

	return(1);

}


/**
 * Write data that is distributed over a parallel universe, with each parallel node independently writing its
 * part to the file.
 * @param  zdf      Parallel ZDF file handle
 * @param  dataset  Dataset to write. Dataset parameters describe global data
 * @param  chunk    Local dataset chunk (includes pointer to data)
 * @param  offset_  Data offset in file layout (in units of data elements). If set to a negative value the code
 *                  will calculate this from individual chunk sizes.
 * @return          Returns 1 if successful, 0 otherwise.
 */
size_t zdf_par_write_cdset_indep( t_zdf_par_file* zdf, t_zdf_dataset* dataset, t_zdf_chunk* chunk, int64_t const offset_ ){

	// Get local rank
	int rank; MPI_Comm_rank( zdf->comm, &rank );

	// Only root node kept track of file position, broadcast it to every node
	int64_t fpos = zdf->fpos;
	MPI_Bcast( &fpos, 1, MPI_INT64_T, 0, zdf -> comm );

	// Get data offsets if needed
	off_t offset = (offset_<0) ? zdf_par_getoffset( chunk, dataset->ndims, zdf->comm ) : offset_;

	// Get the file offset for this rank
	size_t headerlen = size_zdf_chunk_header( dataset );
	size_t data_size = zdf_sizeof(dataset -> data_type);
	off_t file_offset = fpos  + offset * data_size + rank * headerlen;

	// Position file pointer
	if ( fseeko( zdf -> zdf_file.fp, file_offset, SEEK_SET ) == -1 ) {
		perror("(*error*) Unable to position file pointer for parallel I/O.");
		return(0);
	};

	// Write data
	if ( zdf_write_cdset( &zdf -> zdf_file, dataset, chunk ) == 0 ) {
		fprintf(stderr,"(*error*) Unable to write data using parallel I/O.");
		return(0);
	};

	// Set the file position variable to the end of the file on all nodes
	int size; MPI_Comm_size( zdf -> comm, &size );
	size_t count = 1; for( unsigned i = 0; i < dataset -> ndims; i++ ) count *= dataset -> count[i];
	zdf -> fpos = fpos + count * data_size + size * headerlen;

	// Set file position to the end of the file on root node
	if (!rank) {
		file_offset = zdf -> fpos;
		if ( fseeko( zdf -> zdf_file.fp, file_offset, SEEK_SET ) == -1 ) {
			perror("(*error*) Unable to position file pointer at the end of the file on root node.");
			return(0);
		};
	}

	return( 1 );
}

#define _USE_FILE_VIEW 1

#ifdef _USE_FILE_VIEW

/**
 * Write parallel dataset using MPI/IO and MPI_File_set_view.
 * @param  zdf      Parallel ZDF file handle
 * @param  dataset  Dataset to write. Dataset parameters describe global data
 * @param  chunk    Local dataset chunk (includes pointer to data)
 * @param  offset_  Data offset in file layout (in units of data elements). If set to a negative value the code
 *                  will calculate this from individual chunk sizes.
 * @param  coll     Use collective (_all) write operations if set to true
 * @return         	Returns 1 on sucesse, 0 on error.
 */


size_t zdf_par_write_cdset_mpiio( t_zdf_par_file* zdf, t_zdf_dataset* dataset, t_zdf_chunk* chunk, int64_t const offset_, int const coll ){

	// Get local rank
	int rank; MPI_Comm_rank( zdf->comm, &rank );

	// Write chunk header to memory buffer
	uint8_t buf[TMP_BUFFER_SIZE];
	t_zdf_file tmp = { .fp = fmemopen( (void *) buf, TMP_BUFFER_SIZE, "w") };
	size_t headerlen;
	if ( !(headerlen = zdf_write_chunk_header( &tmp, dataset, chunk ) ) ) return(0);
	fclose(tmp.fp);

	size_t count = 1;
	for(unsigned i = 0; i < dataset -> ndims; i++)
		count *= chunk -> count[i];

	MPI_Datatype type[2] = { MPI_BYTE, zdf_mpi_type( dataset -> data_type ) };
	int blocklen[2] = { headerlen, count };

	MPI_Aint disp[2];
	MPI_Get_address( buf, &disp[0]);
	MPI_Get_address( chunk -> data, &disp[1]);
	disp[1] -= disp[0];
	disp[0] = 0;
	MPI_Datatype mem_type;
	MPI_Type_create_struct (2, blocklen, disp, type, &mem_type);
	MPI_Type_commit( &mem_type );

	disp[1] = headerlen;
	MPI_Datatype file_type;
	MPI_Type_create_struct (2, blocklen, disp, type, &file_type);
	MPI_Type_commit( &file_type );

	// Get data and file offsets
	off_t offset = (offset_<0) ? zdf_par_getoffset( chunk, dataset->ndims, zdf->comm ) : offset_;
	size_t data_size = zdf_sizeof(dataset -> data_type);
	MPI_Offset file_offset = zdf -> fpos + offset * data_size + rank * headerlen;

	MPI_File_set_view( zdf -> fh, file_offset, mem_type, file_type, "native", MPI_INFO_NULL);
	if ( coll ) {
		MPI_File_write_all( zdf -> fh, (void *) buf, 1, mem_type, MPI_STATUS_IGNORE );
	} else {
		MPI_File_write( zdf -> fh, (void *) buf, 1, mem_type, MPI_STATUS_IGNORE );
	}

	// Restore default view
	count = 1;
	for(unsigned i = 0; i < dataset -> ndims; i++) count *= dataset -> count[i];
	int size; MPI_Comm_size( zdf->comm, &size );
	zdf -> fpos += count * data_size + size * headerlen;

	MPI_File_set_view( zdf -> fh, zdf -> fpos, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);

	return( 1 );
}


#else 

/**
 * Write parallel dataset using MPI/IO and MPI_File_write_at.
 * @param  zdf     Parallel ZDF file handle
 * @param  dataset  Dataset to write. Dataset parameters describe global data
 * @param  chunk    Local dataset chunk (includes pointer to data)
 * @param  offset_  Data offset in file layout (in units of data elements). If set to a negative value the code
 *                  will calculate this from individual chunk sizes.
 * @param  coll     Use collective (_all) write operations if set to true
 */

size_t zdf_par_write_cdset_mpiio( t_zdf_par_file* zdf, t_zdf_dataset* dataset, t_zdf_chunk* chunk, int64_t const offset_, int const all ){

	// Get local rank
	int rank; MPI_Comm_rank( zdf->comm, &rank );

	// Write chunk header to memory buffer
	uint8_t buf[TMP_BUFFER_SIZE];
	t_zdf_file tmp = { .fp = fmemopen( (void *) buf, TMP_BUFFER_SIZE, "w") };
	size_t headerlen;
	if ( !(headerlen = zdf_write_chunk_header( &tmp, dataset, chunk ) ) ) return(0);
	fclose(tmp.fp);

	// Get size of local chunk
	size_t count = 1;
	for(unsigned i = 0; i < dataset -> ndims; i++) count *= chunk -> count[i];

	// Create MPI datatype representing the chunk header and the chunk data
	MPI_Aint disp[2];
	MPI_Get_address( buf, &disp[0]);
	MPI_Get_address( chunk -> data, &disp[1]);
	disp[1] -= disp[0];
	disp[0] = 0;
	MPI_Datatype type[2] = { MPI_BYTE,  zdf_mpi_type( dataset -> data_type ) };
	int blocklen[2] = { headerlen, count };
	MPI_Datatype mem_type;
	MPI_Type_create_struct (2, blocklen, disp, type, &mem_type);
	MPI_Type_commit( &mem_type );

	// Get data and file offsets
	off_t offset = (offset_<0)?zdf_par_getoffset( chunk, dataset->ndims, zdf->comm ):offset_;
	size_t data_size = zdf_sizeof(dataset -> data_type);
	MPI_Offset file_offset = zdf -> fpos + offset * data_size + rank * headerlen;
	
	// Write header and data in a single call
	if (all) {
		MPI_File_write_at_all( zdf -> fh, offset, (void *) buf, 1, mem_type, MPI_STATUS_IGNORE );
	} else {
		MPI_File_write_at( zdf -> fh, offset, (void *) buf, 1, mem_type, MPI_STATUS_IGNORE );
	}
	// Free MPI datatype
	MPI_Type_free( &mem_type );

			MPI_Barrier( zdf->comm );

	// Set the file position variable to the end of the file on all nodes
	count = 1;
	for(unsigned i = 0; i < dataset -> ndims; i++) count *= dataset -> count[i];
	int size; MPI_Comm_size( zdf->comm, &size );
	zdf -> fpos += count * data_size + size * headerlen;

	return( 1 );
}

#endif

/**
 * Returns the local offset in the file layout, in number of data elements, of the parallel chunk. 
 * Note that this does not include the record headers.
 * @param  chunk (parallel) chunk data
 * @param  ndims Number of dimensions of chunk data
 * @param  comm  Parallel communicator
 * @return       Offset in number of data elements
 */
int64_t zdf_par_getoffset( t_zdf_chunk* chunk, uint32_t ndims, MPI_Comm comm ) {

	int64_t len = chunk -> count[0];

	for( unsigned i = 1; i < ndims; i ++ ) {
		len *= chunk -> count[i];
	}

	int rank; MPI_Comm_rank( comm, &rank );

	int64_t offset;
	MPI_Exscan( (void *) & len, (void *) & offset, 1, MPI_INT64_T, MPI_SUM, comm );

	// The value on rank 0 from MPI_Exscan is undefined so set it explicitly
	if ( rank == 0 ) offset = 0;

	return(offset);
}

/**
 * Adds dataset chunks from a parallel grid to file
 * @param  zdf      Parallel ZDF file handle
 * @param  dataset  Dataset object describing global data
 * @param  chunk    Parallel chunk object with local data
 * @param  offset   Data offset in file layout (in units of data elements). If set to a negative value the code
 *                  will calculate this from individual chunk sizes.
 * @return          Number of bytes written on success, 0 on error. The dataset object is also
 *                  modified to include a unique id.
 */
int zdf_par_write_cdset( t_zdf_par_file* zdf, t_zdf_dataset* dataset, t_zdf_chunk* chunk, const int64_t offset ){

    switch( zdf -> iomode ) {
    	case ZDF_MPI:
    		if( !zdf_par_write_cdset_mpi( zdf, dataset, chunk ) ) return(0);
    		break;
		case ZDF_INDEPENDENT:
    		if( !zdf_par_write_cdset_indep( zdf, dataset, chunk, offset ) ) return(0);
			break;
   		case ZDF_MPIIO_INDEPENDENT:
			if( !zdf_par_write_cdset_mpiio( zdf, dataset, chunk, offset, 0 ) ) return(0);
   			break;
   		case ZDF_MPIIO_COLLECTIVE:
   			if( !zdf_par_write_cdset_mpiio( zdf, dataset, chunk, offset, 1 ) ) return(0);
   			break;
		default :
			fprintf(stderr, "(*error*) Invalid parallel IO mode, unable to write dataset.\n");
			return(0);
    }
 	return( 1 );
}

/*****************************************************************************************************************
 * Fortran wrappers
 */

/**
 * Fortran wrapper for zdf_par_open_file
 * @param  zdf      Parallel ZDF file to open
 * @param  filename Filename of the ZDF file to open, including path
 * @param  mode     Can be one of ZDF_CREATE, ZDF_READ, or ZDF_UPDATE
 * @param  comm_f   Fortran MPI Communicator
 * @return          Returns 0 on success, -1 otherwise
 */
int zdf_par_open_file_f( t_zdf_par_file_f* zdf_f, char* filename, enum zdf_file_access_mode amode, 
	MPI_Fint comm_f, enum zdf_parallel_io_mode iomode ) {
	// Translate Fortran communicator into C communicator
	MPI_Comm comm = MPI_Comm_f2c( comm_f );

	t_zdf_par_file zdf;
	if ( !zdf_par_open_file( &zdf, filename, amode, comm, iomode ) ) return(0);

	zdf_f -> zdf_file = zdf.zdf_file;
	zdf_f -> iomode = iomode;
	zdf_f -> comm = comm_f;

	if (( iomode == ZDF_MPIIO_INDEPENDENT ) || ( iomode == ZDF_MPIIO_COLLECTIVE ) )
		zdf_f -> fh = MPI_File_c2f( zdf.fh );
	else
		zdf_f -> fh = -1;

	zdf_f -> fpos = zdf.fpos;

	return(1);
}


#define PAR_FILE_F2C( zdf_f ) { \
	.iomode = zdf_f -> iomode, \
	.comm = MPI_Comm_f2c( zdf_f -> comm ), \
	.fh = ((zdf_f->iomode == ZDF_MPIIO_INDEPENDENT)||(zdf_f->iomode == ZDF_MPIIO_COLLECTIVE)) ? MPI_File_f2c( zdf_f -> fh) : MPI_FILE_NULL, \
	.fpos = zdf_f -> fpos,\
	.zdf_file = zdf_f -> zdf_file\
}

int zdf_par_close_file_f( t_zdf_par_file_f* zdf_f ) {

	t_zdf_par_file zdf = PAR_FILE_F2C( zdf_f );

	if( ! zdf_par_close_file( &zdf ) ) return(0);
	zdf_f -> fpos = zdf.fpos;
	return(1);
}

int zdf_par_add_string_f( t_zdf_par_file_f* zdf_f, const char* name, const char* str ) {

	t_zdf_par_file zdf = PAR_FILE_F2C( zdf_f );

	if( ! zdf_par_add_string(&zdf, name, str ) ) return(0);
	zdf_f -> fpos = zdf.fpos;
	return(1);
}

int zdf_par_add_int32_f( t_zdf_par_file_f* zdf_f, const char* name, const int32_t value ) {

	t_zdf_par_file zdf = PAR_FILE_F2C( zdf_f );

	if( ! zdf_par_add_int32(&zdf, name, value ) ) return(0);
	zdf_f -> fpos = zdf.fpos;
	return(1);
}

int zdf_par_add_double_f( t_zdf_par_file_f* zdf_f, const char* name, const double value ) {

	t_zdf_par_file zdf = PAR_FILE_F2C( zdf_f );

	if( ! zdf_par_add_double(&zdf, name, value ) ) return(0);
	zdf_f -> fpos = zdf.fpos;
  return(1);
}

int zdf_par_add_iteration_f( t_zdf_par_file_f* zdf_f, const t_zdf_iteration* iter ) {

	t_zdf_par_file zdf = PAR_FILE_F2C( zdf_f );

	if( ! zdf_par_add_iteration(&zdf, iter ) ) return(0);
	zdf_f -> fpos = zdf.fpos;
	return(1);
}

int zdf_par_add_grid_info_f( t_zdf_par_file_f* zdf_f, const t_zdf_grid_info* grid ) {

	t_zdf_par_file zdf = PAR_FILE_F2C( zdf_f );

	if( ! zdf_par_add_grid_info(&zdf, grid ) ) return(0);
	zdf_f -> fpos = zdf.fpos;
	return(1);
}

int zdf_par_add_part_info_f( t_zdf_par_file_f* zdf_f, t_zdf_part_info* part ) {

	t_zdf_par_file zdf = PAR_FILE_F2C( zdf_f );

	if( ! zdf_par_add_part_info(&zdf, part ) ) return(0);
	zdf_f -> fpos = zdf.fpos;
	return(1);
}

int zdf_par_start_cdset_f( t_zdf_par_file_f* zdf_f, t_zdf_dataset* dataset ) {

	t_zdf_par_file zdf = PAR_FILE_F2C( zdf_f );

	if( ! zdf_par_start_cdset(&zdf, dataset ) ) return(0);
	zdf_f -> fpos = zdf.fpos;
	zdf_f -> zdf_file = zdf.zdf_file;
	return(1);

}

int zdf_par_write_cdset_f( t_zdf_par_file_f* zdf_f, t_zdf_dataset* dataset, t_zdf_chunk* chunk, const int32_t offset ){

	t_zdf_par_file zdf = PAR_FILE_F2C( zdf_f );

	int64_t off64 = offset;

	if( ! zdf_par_write_cdset(&zdf, dataset, chunk, off64 ) ) return(0);
	zdf_f -> fpos = zdf.fpos;
	return(1);

}

int zdf_par_end_cdset_f( t_zdf_par_file_f* zdf_f, t_zdf_dataset* dataset ) {

	t_zdf_par_file zdf = PAR_FILE_F2C( zdf_f );

	if( ! zdf_par_end_cdset(&zdf, dataset ) ) return(0);
	zdf_f -> fpos = zdf.fpos;
	return(1);

}
