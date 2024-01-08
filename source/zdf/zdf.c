/*
Copyright (C) 2017 Instituto Superior Tecnico

This file is part of the ZPIC Educational code suite

The ZPIC Educational code suite is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The ZPIC Educational code suite is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with the ZPIC Educational code suite. If not, see <http://www.gnu.org/licenses/>.
*/


/**
 * ZDF version 1
 *
 * This ZDF version is totally self contained. It does not depend on XDR or any other
 * external libraries.
 *
 * Current implementation should also work on big endian systems (untested)
 * This version is not compatible is version 0 (files have the opposite endianess)
 *
 */


#define _POSIX_C_SOURCE 200809L
#define _FILE_OFFSET_BITS 64

#include "zdf.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include <errno.h>

/**
 * On Windows we cannot use the POSIX.1 mkdir command, so use _mkdir instead
 */
#if defined(_MSC_VER) || defined(_WIN32) || defined(_WIN64)
#include <direct.h>
#define mkdir(path,mode) _mkdir(path)

// Windows also lacks the 'fseeko', 'ftello' functions.. so replace with Windows
//   equivalent.
#define fseeko _fseeki64
#define ftello _ftelli64

#endif

/**
 * Number of bytes required for writing data to ZDF file
 * (round up to multiple of BYTES_PER_ZDF_UNIT)
 */
#define RNDUP(x)  ((((x) + BYTES_PER_ZDF_UNIT - 1) / BYTES_PER_ZDF_UNIT) \
		    * BYTES_PER_ZDF_UNIT)


/**
 * Buffer size for converting big endian to little endian data
 */
#define ENDIAN_CONV_BUF_SIZE 1024

/**
 * Sizes of datatypes
 */
const unsigned size_zdf_int32   = 4;
const unsigned size_zdf_uint32  = 4;
const unsigned size_zdf_int64   = 8;
const unsigned size_zdf_uint64  = 8;
const unsigned size_zdf_double  = 8;
const unsigned size_zdf_float   = 4;


/**
 * IDs of ZDF records
 */
#define ZDF_INT32_ID     	0x00010000
#define ZDF_DOUBLE_ID    	0x00020000
#define ZDF_STRING_ID    	0x00030000

#define ZDF_DATASET_ID      0x00100002
#define ZDF_CDSET_START_ID  0x00110000
#define ZDF_CDSET_CHUNK_ID  0x00120000
#define ZDF_CDSET_END_ID	0x00130000

#define ZDF_ITERATION_ID 	0x00200001
#define ZDF_GRID_INFO_ID 	0x00210001
#define ZDF_PART_INFO_ID 	0x00220002
#define ZDF_TRACK_INFO_ID 	0x00230001

/* -----------------------------------------------------------------------------------------------
  recursively create path if required
-------------------------------------------------------------------------------------------------- */

/**
 * Recursively creates a path if required
 * @param  path path to create
 * @return      returns 0 on success, otherwise returns the value returned by mkdir
 */
int create_path( const char path[] )
{
	char uppath[256], *p;
	int ierr = 0;

	if (mkdir(path,S_IRWXU | (S_IRGRP | S_IXGRP ) | (S_IROTH | S_IXOTH) )) {
		switch (errno) {
			case ENOENT : // A component of the path does not exist 

				// get upper path
				strncpy(uppath, path, 256);
				p = uppath + strlen(uppath);
				while(*p!='/') p--;
				*p=0;

				//recursively build the path
				if ( !create_path( uppath ) ) ierr = create_path( path );
				break;

			case EEXIST : /* if directory already exists ignore the error */
				ierr = 0;
				break;
			default: ierr = errno;
		}
	}
	return ierr;
}


/**
 * Returns size of ZDF datatype
 * @param  data_type Data type id
 * @return           Returns type size in bytes or 0 for an invalid data type
 */
size_t zdf_sizeof( enum zdf_data_type data_type ) {
	switch( data_type ) {
        case zdf_null:
        	return(0);
        case zdf_int8:
        case zdf_uint8:
        	return(1);
        case zdf_int16:
        case zdf_uint16:
        	return(2);
        case zdf_int32:
        case zdf_uint32:
        case zdf_float32:
        	return(4);
        case zdf_int64:
        case zdf_uint64:
        case zdf_float64:
        	return(8);
	}
	return(0);
}

/* -----------------------------------------------------------------------------------------------
  Open / Close ZDF file
-------------------------------------------------------------------------------------------------- */

/**
 * Closes ZDF file
 * @param  zdf ZDF file to close
 * @return     Returns 1 on success, 0 on error
 */
int zdf_close_file( t_zdf_file* zdf ) {

	if ( fclose( zdf->fp ) ) {
		perror("(*error*) Unable to close ZDF file");
		return(0);
	}
	zdf -> fp = NULL;
	return(1);
}

/**
 * Opens ZDF file
 * @param  zdf      ZDF file to open
 * @param  filename Filename of the ZDF file to open, including path
 * @param  mode     Can be one of ZDF_CREATE, ZDF_READ, or ZDF_UPDATE
 * @return          Returns 1 on success, 0 on error
 */
int zdf_open_file( t_zdf_file* zdf, const char* filename, enum zdf_file_access_mode mode ){

	char test_magic[4];

	zdf -> mode = mode;

	switch( mode ) {
		case ZDF_CREATE :
			// Open file for writing
			// The "wb" mode must be used for compatibility with Windows
			if (!(zdf->fp = fopen( filename, "w+b"))) {
				perror("(*error*) Unable to open ZDF file for writing");
				return(0);
			}

			// Write magic number
			if ( fwrite( (void *) zdf_magic, sizeof(char), ZDF_MAGIC_LENGTH, zdf->fp )
			     != ZDF_MAGIC_LENGTH ) {
				fprintf(stderr, "(*error*) Unable to write magic number to ZDF file.\n");
				zdf_close_file( zdf );
				return(0);
			}

			break;

		case ZDF_READ :
			// Open file for reading
			if (!(zdf->fp = fopen( filename, "r"))) {
				perror("(*error*) Unable to open ZDF file for reading");
				return(0);
			}

			// Read magic number
			if ( fread( (void *) test_magic, sizeof(char), ZDF_MAGIC_LENGTH, zdf->fp )
				!= ZDF_MAGIC_LENGTH) {
				fprintf(stderr, "(*error*) Unable to read magic number from ZDF file.\n");
				zdf_close_file( zdf );
				return(0);
			}

			// Check magic number
			for( int i = 0; i < ZDF_MAGIC_LENGTH; i++) {
				if ( test_magic[i] != zdf_magic[i] ) {
					fprintf(stderr, "(*error*) Invalid magic number, file is not a proper ZDF file.\n");
					zdf_close_file( zdf );
					return(0);
				}
			}
			break;

		case ZDF_UPDATE :
			// Open file for reading and writing
			if (!(zdf->fp = fopen( filename, "r+"))) {
				perror("(*error*) Unable to open ZDF file for reading / writing.\n");
				return(0);
			}

			// Read magic number
			if (fread( (void *) test_magic, sizeof(char), ZDF_MAGIC_LENGTH, zdf->fp )
				!= ZDF_MAGIC_LENGTH) {
				fprintf(stderr, "(*error*) Unable to read magic number from ZDF file.\n");
				zdf_close_file( zdf );
				return(0);
			}

			// Check magic number
      		for( int i = 0; i < ZDF_MAGIC_LENGTH; i++) {
	    		if ( test_magic[i] != zdf_magic[i] ) {
					fprintf(stderr, "(*error*) Invalid magic number, file is not a proper ZDF file.\n");
					zdf_close_file( zdf );
					return(0);
				}
			}

			// Position the file pointer at the end of the file
			fseeko( zdf->fp, 0, SEEK_END );

			break;

		default:
			fprintf(stderr, "(*error*) zdf_open_file: unsupported mode.\n");
			return(0);
	}

	zdf -> ndatasets = 0;

	return(1);
}

/* -----------------------------------------------------------------------------------------------
  Elemental types
-------------------------------------------------------------------------------------------------- */


#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__

/**
 * Implementation for little endian systems (e.g. x86)
 *
 * For these systems just writing the data to disk is sufficient
 */

/**
 * Writes scalar int32 value to file
 * @param  zdf ZDF file descriptor
 * @param  i   int32_t value to write
 * @return     Returns number of bytes written on success, 0 on error
 */
size_t zdf_int32_write( t_zdf_file* zdf, const int32_t i ){
	if ( fwrite( (void *) &i, sizeof(int32_t), 1, zdf -> fp ) != 1 )
		return(0);
	return ( sizeof(int32_t) );
}

/**
 * Reads scalar int32 value from file
 * @param  zdf ZDF file descriptor
 * @param  i   int32_t value read
 * @return     Returns number of bytes written on success, 0 on error
 */
size_t zdf_int32_read( t_zdf_file* zdf, int32_t* i ){
	if ( fread( (void *) i, sizeof(int32_t), 1, zdf -> fp ) != 1 )
		return(0);
	return ( sizeof(int32_t) );
}

/**
 * Writes scalar uint32 value to file
 * @param  zdf ZDF file descriptor
 * @param  u   uint32_t value to write
 * @return     Returns number of bytes written on success, 0 on error
 */
size_t zdf_uint32_write( t_zdf_file* zdf, const uint32_t u ){
	if ( fwrite( (void *) &u, sizeof(uint32_t), 1, zdf -> fp ) != 1 )
		return(0);
	return(sizeof(uint32_t));
}

/**
 * Reads scalar uint32 value from file
 * @param  zdf ZDF file descriptor
 * @param  u   uint32_t value read
 * @return     Returns number of bytes written on success, 0 on error
 */
size_t zdf_uint32_read( t_zdf_file* zdf, uint32_t* u ){
	if ( fread( (void *) u, sizeof(uint32_t), 1, zdf -> fp ) != 1 )
		return(0);
	return(sizeof(uint32_t));
}

/**
 * Writes scalar int64 value to file
 * @param  zdf ZDF file descriptor
 * @param  i   int64_t value to write
 * @return     Returns number of bytes written on success, 0 on error
 */
size_t zdf_int64_write( t_zdf_file* zdf, const int64_t i ){
	if ( fwrite( (void *) &i, sizeof(int64_t), 1, zdf -> fp ) != 1 )
		return(0);
	return(sizeof(int64_t));
}

/**
 * Reads scalar int64 value from file
 * @param  zdf ZDF file descriptor
 * @param  i   int64_t value read
 * @return     Returns number of bytes written on success, 0 on error
 */
size_t zdf_int64_read( t_zdf_file* zdf, int64_t *i ){
	if ( fread( (void *) i, sizeof(int64_t), 1, zdf -> fp ) != 1 )
		return(0);
	return(sizeof(int64_t));
}

/**
 * Writes scalar uint64_t value to file
 * @param  zdf ZDF file descriptor
 * @param  u   uint64_t value to write
 * @return     Returns number of bytes written on success, 0 on error
 */
size_t zdf_uint64_write( t_zdf_file* zdf, const uint64_t u ){
	if ( fwrite( (void *) &u, sizeof(uint64_t), 1, zdf -> fp ) != 1 )
		return(0);
	return(sizeof(uint64_t));
}

/**
 * Reads scalar int64 value from file
 * @param  zdf ZDF file descriptor
 * @param  u   uint64_t value read
 * @return     Returns number of bytes written on success, 0 on error
 */
size_t zdf_uint64_read( t_zdf_file* zdf, uint64_t *u ){
	if ( fread( (void *) u, sizeof(int64_t), 1, zdf -> fp ) != 1 )
		return(0);
	return(sizeof(uint64_t));
}


/**
 * Writes scalar double (float64) value to file
 * @param  zdf ZDF file descriptor
 * @param  d   double value to write
 * @return     Returns number of bytes written on success, 0 on error
 */
size_t zdf_double_write( t_zdf_file* zdf, const double d ){
	if ( fwrite( (void *) &d, sizeof(double), 1, zdf -> fp ) != 1 )
		return(0);
	return(sizeof(double));
}

/**
 * Write arbitraty 16bit type vector to file
 * @param  zdf          ZDF file descriptor
 * @param  data         Pointer to data to write
 * @param  element_size Size of vector element
 * @param  len          Number of vector elements
 * @return              Returns number of bytes written on success, 0 on error
 */
size_t zdf_vector16_write( t_zdf_file* zdf, void const * const data, size_t len ) {
	if ( fwrite( data, sizeof(int16_t), len, zdf -> fp ) != len )
		return(0);
	return( len * sizeof(int16_t) );
}

/**
 * Write arbitraty 32bit type vector to file
 * @param  zdf          ZDF file descriptor
 * @param  data         Pointer to data to write
 * @param  element_size Size of vector element
 * @param  len          Number of vector elements
 * @return              Returns number of bytes written on success, 0 on error
 */
size_t zdf_vector32_write( t_zdf_file* zdf, void const * const data, size_t len ) {
	if ( fwrite( data, sizeof(int32_t), len, zdf -> fp ) != len )
		return(0);
	return( len * sizeof(int32_t) );
}

/**
 * Write arbitraty 64bit type vector to file
 * @param  zdf          ZDF file descriptor
 * @param  data         Pointer to data to write
 * @param  element_size Size of vector element
 * @param  len          Number of vector elements
 * @return              Returns number of bytes written on success, 0 on error
 */
size_t zdf_vector64_write( t_zdf_file* zdf, void const * const data, size_t len ) {
	if( fwrite( data, sizeof(int64_t), len, zdf -> fp ) != len )
		return(0);
	return( len * sizeof(int64_t) );
}

#elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__

/**
 * Implementation for big endian systems (e.g. PowerPC)
 *
 * For these systems we need to swap the bytes before writing. For vectors, conversion
 * is done in chunks of ENDIAN_CONV_BUF_SIZE values.
 *
 * The system may provide hardware optimized bswap_* routines, these are usually found in
 * the <byteswap.h> header, but they are missing on many platforms.
 */


/**
 * Swap bytes for 16 bit number
 * @param  x 16 bit number to convert
 * @return   byte-swapped version of x
 */
#define bswap_16(x) \
({ \
	uint16_t __x = (x); \
	((uint16_t)( \
		(((uint16_t)(__x) & (uint16_t)0x00ffU) << 8) | \
		(((uint16_t)(__x) & (uint16_t)0xff00U) >> 8) )); \
})

/**
 * Swap bytes for 32 bit number
 * @param  x 32 bit number to convert
 * @return   byte-swapped version of x
 */
#define bswap_32(x) \
({ \
	uint32_t __x = (x); \
	((uint32_t)( \
		(((uint32_t)(__x) & (uint32_t)0x000000ffUL) << 24) | \
		(((uint32_t)(__x) & (uint32_t)0x0000ff00UL) <<  8) | \
		(((uint32_t)(__x) & (uint32_t)0x00ff0000UL) >>  8) | \
		(((uint32_t)(__x) & (uint32_t)0xff000000UL) >> 24) )); \
})

/**
 * Swap bytes for 64 bit number
 * @param  x 64 bit number to convert
 * @return   byte-swapped version of x
 */
#define bswap_64(x) \
({ \
	uint64_t __x = (x); \
	((uint64_t)( \
		(uint64_t)(((uint64_t)(__x) & (uint64_t)0x00000000000000ffULL) << 56) | \
		(uint64_t)(((uint64_t)(__x) & (uint64_t)0x000000000000ff00ULL) << 40) | \
		(uint64_t)(((uint64_t)(__x) & (uint64_t)0x0000000000ff0000ULL) << 24) | \
		(uint64_t)(((uint64_t)(__x) & (uint64_t)0x00000000ff000000ULL) <<  8) | \
	    (uint64_t)(((uint64_t)(__x) & (uint64_t)0x000000ff00000000ULL) >>  8) | \
		(uint64_t)(((uint64_t)(__x) & (uint64_t)0x0000ff0000000000ULL) >> 24) | \
		(uint64_t)(((uint64_t)(__x) & (uint64_t)0x00ff000000000000ULL) >> 40) | \
		(uint64_t)(((uint64_t)(__x) & (uint64_t)0xff00000000000000ULL) >> 56) )); \
})

/**
 * Writes int32 value to file
 * @param  zdf ZDF file descriptor
 * @param  i   int32_t value to write
 * @return     Returns number of bytes written on success, 0 on error
 */
size_t zdf_int32_write( t_zdf_file* zdf, const int32_t i ){
	uint32_t tmp = bswap_32((uint32_t)i);
	if ( fwrite( (void *) &tmp, sizeof(uint32_t), 1, zdf -> fp ) != 1 )
		return(0);
	return( sizeof(uint32_t));
}

/**
 * Reads scalar int32 value from file
 * @param  zdf ZDF file descriptor
 * @param  i   int32_t value read
 * @return     Returns number of bytes read on success, 0 on error
 */
size_t zdf_int32_read( t_zdf_file* zdf, int32_t* i ){
	uint32_t tmp;
	if ( fread( (void *) &tmp, sizeof(uint32_t), 1, zdf -> fp ) != 1 )
		return(0);
	*i = (int32_t) bswap_32(tmp);
	return( sizeof(uint32_t));
}

/**
 * Writes uint32 value to file
 * @param  zdf ZDF file descriptor
 * @param  u   uint32_t value to write
 * @return     Returns number of bytes written on success, 0 on error
 */
size_t zdf_uint32_write( t_zdf_file* zdf, const uint32_t u ){
	uint32_t tmp = bswap_32(u);
	if ( fwrite( (void *) &tmp, sizeof(uint32_t), 1, zdf -> fp ) != 1 )
		return(0);
	return(sizeof(uint32_t));
}

/**
 * Reads scalar uint32 value from file
 * @param  zdf ZDF file descriptor
 * @param  i   uint32_t value read
 * @return     Returns number of bytes read on success, 0 on error
 */
size_t zdf_uint32_read( t_zdf_file* zdf, uint32_t* u ){
	uint32_t tmp;
	if ( fread( (void *) &tmp, sizeof(uint32_t), 1, zdf -> fp ) != 1 )
		return(0);
	*u = bswap_32(tmp);
	return(sizeof(uint32_t));
}

/**
 * Writes int64 value to file
 * @param  zdf ZDF file descriptor
 * @param  i   int64_t value to write
 * @return     Returns number of bytes written on success, 0 on error
 */
size_t zdf_int64_write( t_zdf_file* zdf, const int64_t i ){
	uint64_t tmp = bswap_64((uint64_t)i);
	if ( fwrite( (void *) &tmp, sizeof(uint64_t), 1, zdf -> fp ) != 1 )
		return(0);
	return(sizeof(uint64_t));
}

/**
 * Reads scalar int64 value from file
 * @param  zdf ZDF file descriptor
 * @param  i   int64_t value read
 * @return     Returns number of bytes read on success, 0 on error
 */
size_t zdf_int64_read( t_zdf_file* zdf, int64_t* i ){
	uint64_t tmp;
	if ( fread( (void *) &tmp, sizeof(uint64_t), 1, zdf -> fp ) != 1 )
		return(0);
	*i = (int64_t) bswap_64(tmp);
	return(sizeof(uint64_t));
}

/**
 * Writes uint64_t value to file
 * @param  zdf ZDF file descriptor
 * @param  u   uint64_t value to write
 * @return     Returns number of bytes read on success, 0 on error
 */
size_t zdf_uint64_write( t_zdf_file* zdf, const uint64_t u ){
	uint64_t tmp = bswap_64(u);
	if ( fwrite( (void *) &tmp, sizeof(uint64_t), 1, zdf -> fp ) != 1 )
		return(0);
	return(sizeof(uint64_t));
}

/**
 * Reads scalar uint64 value from file
 * @param  zdf ZDF file descriptor
 * @param  i   uint64_t value read
 * @return     Returns number of bytes read on success, 0 on error
 */
size_t zdf_uint64_read( t_zdf_file* zdf, uint64_t* i ){
	uint64_t tmp;
	if ( fread( (void *) &tmp, sizeof(uint64_t), 1, zdf -> fp ) != 1 )
		return(0);
	*i = bswap_64(tmp);
	return(sizeof(uint64_t));
}

/**
 * Writes double (float64) value to file
 * @param  zdf ZDF file descriptor
 * @param  d   double value to write
 * @return     Returns number of bytes written on success, 0 on error
 */
size_t zdf_double_write( t_zdf_file* zdf, const double d ){
	uint64_t tmp = bswap_64((uint64_t)d);
	if ( fwrite( (void *) &tmp, sizeof(uint64_t), 1, zdf -> fp ) != 1 )
		return(0);
	return(sizeof(uint64_t));
}

/**
 * Write arbitraty 16bit type vector to file
 * @param  zdf  ZDF file descriptor
 * @param  data Pointer to 16 bit data to write
 * @param  len  Number of vector elements
 * @return      Returns number of bytes written on success, 0 on error
 */
size_t zdf_vector16_write( t_zdf_file* zdf,  void const * const data, size_t len ) {

	uint16_t buffer[ENDIAN_CONV_BUF_SIZE];

	for( size_t offset = 0; offset < len; offset += ENDIAN_CONV_BUF_SIZE ) {

		// Number of values in chunk
		size_t chunk_len = (offset + ENDIAN_CONV_BUF_SIZE < len ) ? ENDIAN_CONV_BUF_SIZE : len - offset;

		// Convert chunk to little endian
		for( size_t i = 0; i < chunk_len; i++) buffer[i] = bswap_16( data[offset+i]);

		// Write chunk
		if ( fwrite( (void *) buffer, sizeof(uint16_t), chunk_len, zdf -> fp) != chunk_len ) 
			return(0);
    }

	return( len * sizeof(uint16_t) );
}

/**
 * Write arbitraty 32bit type vector to file
 * @param  zdf  ZDF file descriptor
 * @param  data Pointer to 32 bit data to write
 * @param  len  Number of vector elements
 * @return      Returns number of bytes written on success, 0 on error
 */
size_t zdf_vector32_write( t_zdf_file* zdf,  void const * const data, size_t len ) {

	uint32_t buffer[ENDIAN_CONV_BUF_SIZE];

	for( size_t offset = 0; offset < len; offset += ENDIAN_CONV_BUF_SIZE ) {

		// Number of values in chunk
		size_t chunk_len = (offset + ENDIAN_CONV_BUF_SIZE < len ) ? ENDIAN_CONV_BUF_SIZE : len - offset;

		// Convert chunk to little endian
		for( size_t i = 0; i < chunk_len; i++) buffer[i] = bswap_32( data[offset+i]);

		// Write chunk
		if ( fwrite( (void *) buffer, sizeof(uint32_t), chunk_len, zdf -> fp) != chunk_len ) 
			return(0);
    }

	return( len * sizeof(uint32_t) );
}

/**
 * Write arbitraty 64bit type vector to file
 * @param  zdf  ZDF file descriptor
 * @param  data Pointer to 64 bit data to write
 * @param  len  Number of vector elements
 * @return      Returns number of bytes written on success, 0 on error
 */
size_t zdf_vector64_write( t_zdf_file* zdf,  void const * const data, size_t len ) {
	uint64_t buffer[ENDIAN_CONV_BUF_SIZE];

	for( size_t offset = 0; offset < len; offset += ENDIAN_CONV_BUF_SIZE ) {

		// Number of values in chunk
		size_t chunk_len = (offset + ENDIAN_CONV_BUF_SIZE < len ) ? ENDIAN_CONV_BUF_SIZE : len - offset;

		// Convert chunk to little endian
		for( size_t i = 0; i < chunk_len; i++) buffer[i] = bswap_64( data[offset+i] );

		// Write chunk
		if ( fwrite( (void *) buffer, sizeof(uint64_t), chunk_len, zdf -> fp) != chunk_len ) 
			return(0);
    }

	return( len * sizeof(uint64_t) );
}


#else

#error "System is neither little endian nor big endian, aborting."

#endif

/**
 * Write arbitraty 8bit type vector to file. Adds padding at the end of the vector if
 * necessary to maitain alignment.
 * @param  zdf ZDF file descriptor
 * @param  u   Pointer to 8 bit data
 * @param  len Number of vector elements
 * @return     Returns number of bytes written (including padding), 0 on error
 */
size_t zdf_vector8_write( t_zdf_file* zdf, void const * const u, size_t len ){

	if ( fwrite( u, sizeof(uint8_t), len, zdf -> fp ) != len ) {
		return(0);
	}

	size_t npad = RNDUP(len) - len;

	if ( npad > 0 ) {
		const uint8_t pad[BYTES_PER_ZDF_UNIT] = {0};
		if ( fwrite( (void *) pad, sizeof(uint8_t), npad, zdf -> fp ) != npad ) {
			return(0);
		}
	}

	return( RNDUP(len) );
}

/**
 * Writes a vector of the specified datatype to file
 * @param  zdf       ZDF file descriptor
 * @param  data      Pointer to data
 * @param  data_type Data type descriptor
 * @param  len       Number of elements in vector
 * @return           Returns number of bytes written on success, 0 on error
 */
size_t zdf_vector_write( t_zdf_file* zdf, const void * data, enum zdf_data_type data_type, size_t len ){

 	switch ( data_type ) {
    	case zdf_int8:
    	case zdf_uint8:
    		return ( zdf_vector8_write( zdf, data, len ));
    	case zdf_int16:
    	case zdf_uint16:
    	    return ( zdf_vector16_write( zdf, data, len ) );
    	case zdf_int32:
    	case zdf_uint32:
    	case zdf_float32:
    	    return ( zdf_vector32_write( zdf, data, len ) );
    	case zdf_int64:
    	case zdf_uint64:
    	case zdf_float64:
    	    return ( zdf_vector64_write( zdf, data, len ) );
    	default:
    		fprintf(stderr,"(*error*) zdf_vector_write: Unsupported datatype.\n");
 	}
	return(0);
 }


/* -----------------------------------------------------------------------------------------------
  zdf_string
-------------------------------------------------------------------------------------------------- */

/**
 * Write a string to file
 * @param  zdf ZDF file descriptor
 * @param  str C string (char *) to write
 * @return     Returns number of bytes written on success, 0 on error
 */
size_t zdf_string_write( t_zdf_file* zdf, const char * str ){

	uint32_t len;
	int count;

	len = ( str ) ? strlen( str ) : 0;
	if ( ! zdf_uint32_write( zdf, len ) ) return(0);

	if ( len > 0 ) {
		if ( !(count = zdf_vector8_write( zdf, (void *) str, len ) ) )
			return(0);
	} else {
		count = 0;
	}

	return( sizeof(uint32_t) + count );
}

/**
 * Read a string from file
 * @param  zdf ZDF file descriptor
 * @param  str C string (char *) to read. The routine will allocate memory using
 *             malloc() for storing the string data.
 * @return     Returns number of bytes read on success, 0 on error
 */
size_t zdf_string_read( t_zdf_file* zdf, char * * str ){

	uint32_t len;
	if ( ! zdf_uint32_read( zdf, &len) ) return(0);

	uint32_t flen = RNDUP(len);

    char * buffer = (char *) malloc( (flen+1) * sizeof(uint8_t) );

	if ( len > 0 ) {
		if ( fread( (void *) buffer, sizeof(uint8_t), flen, zdf -> fp ) != flen ) return(0);
	}

	buffer[len] = 0;

	*str = buffer;

	return( sizeof(uint32_t) + flen );
}

/**
 * Calculate size of string record
 * @param  s C string (char *) to write
 * @return   The number of bytes required for writing the string data
 */
size_t size_zdf_string( const char *s )
{
	unsigned len = ( s ) ? strlen( s ) : 0;
	return size_zdf_uint32 + ((len > 0) ? RNDUP(len) : 0);
}


/* -----------------------------------------------------------------------------------------------
  zdf records
-------------------------------------------------------------------------------------------------- */

typedef struct {
	uint32_t id_version;
	char*    name;
	uint64_t length;
} t_zdf_record;

/**
 * Adds ZDF record header to file
 * @param  zdf ZDF file descriptor
 * @param  rec ZDF record
 * @return     Returns number of bytes written on success, 0 on error
 */
size_t zdf_record_write( t_zdf_file* zdf, const t_zdf_record* rec ){

  if ( ! zdf_uint32_write( zdf, rec -> id_version ) ) return(0);
  size_t len;	if ( ! (len=zdf_string_write( zdf, rec -> name )) ) return(0);
	if ( ! zdf_uint64_write( zdf, rec -> length )     ) return(0);

 	return( sizeof(uint32_t) + len + sizeof(uint64_t) );
}

/**
 * Reads ZDF record header from file
 * @param  zdf ZDF file descriptor
 * @param  rec ZDF record
 * @return     Returns number of bytes read on success, 0 on error
 */
size_t zdf_record_read( t_zdf_file* zdf, t_zdf_record* rec ){

  if ( ! zdf_uint32_read( zdf, & rec -> id_version ) ) return(0);
  size_t len; if ( ! (len=zdf_string_read( zdf, & rec -> name )) ) return(0);
	if ( ! zdf_uint64_read( zdf, & rec -> length )     ) return(0);

 	return( sizeof(uint32_t) + len + sizeof(uint64_t) );
}

/**
 * Skips to the end of ZDF record. File pointer is assumed to be after the record header.
 * @param  zdf ZDF file descriptor
 * @param  rec ZDF record
 * @return     Returns 1 on success, 0 on error
 */
int zdf_record_skip( t_zdf_file* zdf, const t_zdf_record* rec ){

    off_t offset = rec -> length;
    return ( ! fseeko( zdf -> fp, offset, SEEK_CUR ) ? 1 : 0 );
}

/* -----------------------------------------------------------------------------------------------
  zdf basic data tags
-------------------------------------------------------------------------------------------------- */

/**
 * Adds string element to ZDF file
 * @param  zdf  ZDF File handle
 * @param  name Element name
 * @param  str  String value
 * @return      Returns number of bytes written on success, 0 on error
 */
size_t zdf_add_string( t_zdf_file* zdf, const char* name, const char* str ){

	t_zdf_record rec = {
		.id_version = ZDF_STRING_ID,
		.name = (char *)name,
		.length = size_zdf_string( str )
	};

    size_t ok, len;
    len = ( ok = zdf_record_write( zdf, &rec ) );
    if ( !ok ) return(0);

    len += ( ok = zdf_string_write( zdf, str) );
    if ( !ok ) return(0);

    return(len);
}

/**
 * Adds int32 element to ZDF file
 * @param  zdf    ZDF File handle
 * @param  name   Element name
 * @param  value  int32 value
 * @return        Returns number of bytes written on success, 0 on error
 */
size_t zdf_add_int32( t_zdf_file* zdf, const char* name, const int32_t value ){

	t_zdf_record rec = {
		.id_version = ZDF_INT32_ID,
		.name = (char *)name,
		.length = size_zdf_int32
	};

    size_t ok, len;
    len = ( ok = zdf_record_write(zdf, &rec) );
    if (!ok) return(0);

    len += ( ok = zdf_int32_write(zdf, value) );
    if (!ok) return(0);

    return( len );
}

/**
 * Adds float64 element to ZDF file
 * @param  zdf    ZDF file handle
 * @param  name   Element name
 * @param  value  float64 value
 * @return        Returns number of bytes written on success, 0 on error
 */
size_t zdf_add_double( t_zdf_file* zdf, const char* name, const double value )
{
	t_zdf_record rec = {
		.id_version = ZDF_DOUBLE_ID,
		.name = (char *)name,
		.length = size_zdf_double
	};

    size_t ok, len;
    len = (ok = zdf_record_write(zdf, &rec) );    if ( !ok ) return(0);
    len += ( ok = zdf_double_write(zdf, value) ); if ( !ok ) return(0);

    return(len);
}


/* -----------------------------------------------------------------------------------------------
  zdf compound metadata tags
-------------------------------------------------------------------------------------------------- */

/**
 * Adds iteration metadata group to ZDF file
 * @param  zdf  ZDF file handle
 * @param  iter Iteration info
 * @return      Number of bytes written on success, 0 on error
 */
size_t zdf_add_iteration( t_zdf_file* zdf, const t_zdf_iteration* iter ){

	t_zdf_record rec = {
		.id_version = ZDF_ITERATION_ID,
		.name = (char *) iter -> name,
		.length = size_zdf_uint32 +
		             size_zdf_double +
		             size_zdf_string( iter -> time_units )
	};

    size_t ok, len;

    len = ( ok = zdf_record_write( zdf, &rec) );      if ( !ok ) return(0);
 	len += ( ok = zdf_int32_write( zdf, iter ->n ) ); if ( !ok ) return(0);
	len += ( ok = zdf_double_write( zdf, iter->t ) ); if ( !ok ) return(0);
 	len += ( ok = zdf_string_write( zdf, iter->time_units ) ); if ( !ok ) return(0);

 	return(len);
}

/**
 * Returns size of grid information metadata element
 * @param  grid Grid info
 * @return Metadata group size in bytes
 */
size_t size_zdf_grid_info(const t_zdf_grid_info* grid) {
	size_t size;

	size = size_zdf_uint32 + grid -> ndims * size_zdf_uint64 +
	        size_zdf_string(grid->label) +  size_zdf_string(grid->units);

   	// Includes axis information
   	size += size_zdf_int32;
    if ( grid -> axis ) {

    	for(unsigned i=0; i<grid -> ndims; i++)
    		size += size_zdf_string( grid->axis[i].name ) +
                    size_zdf_int32 +
    	            2 * size_zdf_double +
    	            size_zdf_string( grid -> axis[i].label ) +
    	            size_zdf_string( grid -> axis[i].units );
    }

	return size;
}

/**
 * Adds grid information metadata group to file
 * @param  zdf  File handle
 * @param  grid Grid information
 * @return      Number of bytes written on success, 0 on error
 */
size_t zdf_add_grid_info( t_zdf_file* zdf, const t_zdf_grid_info* grid ){

	size_t reclen;

	t_zdf_record rec = {
		.id_version = ZDF_GRID_INFO_ID,
		.name = (char *) grid -> name,
		.length = size_zdf_grid_info( grid )
	};

  if ( !( reclen = zdf_record_write( zdf, &rec ) )  ) return(0);
 	if ( ! zdf_uint32_write( zdf, grid -> ndims ) ) return(0);
    for( unsigned i=0; i < grid -> ndims; i++) {
	 	if ( ! zdf_uint64_write( zdf, grid -> count[i] ) ) return(0);
  }

 	if ( !zdf_string_write( zdf, grid->label ) ) return(0);
 	if ( !zdf_string_write( zdf, grid->units ) ) return(0);

  int32_t has_axis = ( grid -> axis != NULL);
 	if ( !zdf_int32_write( zdf, has_axis ) ) return(0);

 	if ( has_axis ) {

	    for( unsigned i=0; i < grid -> ndims; i++) {
            if ( !zdf_string_write( zdf, grid -> axis[i].name  ) ) return(0);
		    if ( !zdf_int32_write(  zdf, grid -> axis[i].type  ) ) return(0);
		    if ( !zdf_double_write( zdf, grid -> axis[i].min   ) ) return(0);
		    if ( !zdf_double_write( zdf, grid -> axis[i].max   ) ) return(0);
		    if ( !zdf_string_write( zdf, grid -> axis[i].label ) ) return(0);
		    if ( !zdf_string_write( zdf, grid -> axis[i].units ) ) return(0);
	    }
 	}

 	return( reclen + size_zdf_grid_info( grid ) );
}

/**
 * Returns size of particles information metadata group
 * @param  part Particles information
 * @return      Metadata group size in bytes
 */
size_t size_zdf_part_info(const t_zdf_part_info* part) {
	size_t size = size_zdf_string(part->label) +	// label
			size_zdf_uint64 + 						// np
			size_zdf_uint32; 						// nquants

	for( unsigned i = 0; i < part -> nquants; i++) {
		size += size_zdf_string( part -> quants[i] ) +
				size_zdf_string( part -> qlabels[i] ) +
				size_zdf_string( part -> qunits[i] );
	}

	return size;
}

/**
 * Adds particle information metadata group to file
 * @param  zdf  File handle
 * @param  part Particle information
 * @return      Number of bytes written on success, 0 on error
 */
size_t zdf_add_part_info( t_zdf_file* zdf, const t_zdf_part_info* part ){

	size_t reclen;

	t_zdf_record rec = {
		.id_version = ZDF_PART_INFO_ID,
		.name = (char *) part->name,
		.length = size_zdf_part_info( part )
	};

    if ( !( reclen = zdf_record_write( zdf, &rec ) ) ) return(0);
	if ( !zdf_string_write( zdf, part -> label )     ) return(0);
 	if ( !zdf_uint64_write( zdf, part -> np )        ) return(0);
 	if ( !zdf_uint32_write( zdf, part -> nquants )   ) return(0);

    for( unsigned i=0; i < part -> nquants; i++) {
		if ( !zdf_string_write( zdf, part -> quants[i] ) ) return(0);
    }
    for( unsigned i=0; i < part -> nquants; i++) {
		if ( !zdf_string_write( zdf, part -> qlabels[i] ) ) return(0);
    }
    for( unsigned i=0; i < part -> nquants; i++) {
		if ( !zdf_string_write( zdf, part -> qunits[i] )  ) return(0);
    }

 	return( reclen + size_zdf_part_info( part ) );
}

/**
 * Returns size of tracks information metadata group
 * @param  tracks Tracks information
 * @return        Metadata group size in bytes
 */
size_t size_zdf_track_info(const t_zdf_track_info* tracks) {
	size_t size = size_zdf_string(tracks->label) +	// label
			4 * size_zdf_uint32; // ntracks, ndump, niter, nquants

	for( unsigned i = 0; i < tracks -> nquants; i++) {
		size += size_zdf_string( tracks -> quants[i] ) +
				size_zdf_string( tracks -> qlabels[i] ) +
				size_zdf_string( tracks -> qunits[i] );
	}

	return size;
}

/**
 * Adds tracks information metadata group to file
 * @param  zdf    File handle
 * @param  tracks Tracks information
 * @return        Number of bytes written on success, 0 on error
 */
size_t zdf_add_track_info( t_zdf_file* zdf, const t_zdf_track_info* tracks ){

	int reclen;

	t_zdf_record rec = {
		.id_version = ZDF_TRACK_INFO_ID,
		.name = (char *)tracks -> name,
		.length = size_zdf_track_info( tracks )
	};

	if ( !( reclen = zdf_record_write( zdf, &rec ) ) ) return(0);
	if ( !zdf_string_write( zdf, tracks -> label )   ) return(0);
 	if ( !zdf_uint32_write( zdf, tracks -> ntracks ) ) return(0);
 	if ( !zdf_uint32_write( zdf, tracks -> ndump )   ) return(0);
 	if ( !zdf_uint32_write( zdf, tracks -> niter )   ) return(0);
 	if ( !zdf_uint32_write( zdf, tracks -> nquants ) ) return(0);

	for( unsigned i=0; i < tracks -> nquants; i++) {
		if ( !zdf_string_write( zdf, tracks -> quants[i] ) ) return(0);
	}
	for( unsigned i=0; i < tracks -> nquants; i++) {
		if ( !zdf_string_write( zdf, tracks -> qlabels[i] ) ) return(0);
	}
	for( unsigned i=0; i < tracks -> nquants; i++) {
		if ( !zdf_string_write( zdf, tracks -> qunits[i] ) ) return(0);
	}

 	return( reclen + size_zdf_track_info( tracks ) );
}
/* -----------------------------------------------------------------------------------------------
  zdf dataset
-------------------------------------------------------------------------------------------------- */

/**
 * [size_zdf_dataset_header description]
 * @param  dataset [description]
 * @return         [description]
 */
size_t size_zdf_dataset_header( const t_zdf_dataset* dataset ) {
	return size_zdf_uint32 + size_zdf_int32 + size_zdf_uint32 +
		dataset -> ndims * size_zdf_uint64;
}

/**
 * Write dataset header to file and modify the dataset object to include the
 * file position where the dataset is stored.
 *
 * @param  zdf     ZDF file handle
 * @param  dataset Dataset object
 * @return         Number of bytes written on success, 0 on error
 */
size_t zdf_dataset_header_write( t_zdf_file* zdf, t_zdf_dataset* dataset ) {

	// Get the file position for the dataset header and store it
	off_t offset = ftello( zdf -> fp );
	if ( offset < 0 ) return(0);
	dataset -> offset = offset;

 	// Version 0x0001
 	uint32_t id = dataset -> id;
 	if ( !zdf_uint32_write( zdf, id )                  ) return(0);
 	if ( !zdf_int32_write( zdf, dataset -> data_type ) ) return(0);
	if ( !zdf_uint32_write( zdf, dataset -> ndims )    ) return(0);

    for( unsigned i=0; i < dataset -> ndims; i++) {
	 	if ( !zdf_uint64_write( zdf, dataset -> count[i] ) ) return(0);
    }
    return( size_zdf_dataset_header( dataset ) );
}

/**
 * Read dataset header from file and store the file position where the dataset
 * header is stored.
 *
 * @param  zdf     ZDF file handle
 * @param  dataset Dataset object
 * @return         Number of bytes written on success, 0 on error
 */
size_t zdf_dataset_header_read( t_zdf_file* zdf, t_zdf_dataset* dataset ) {

	// Get the file position for the dataset header and store it
	off_t offset = ftello( zdf -> fp );
	if ( offset < 0 ) return(0);
	dataset -> offset = offset;

 	// Version 0x0001
 	uint32_t id;
	if ( !zdf_uint32_read( zdf, &id ) ) return(0);
	dataset -> id = id;

 	if ( !zdf_int32_read( zdf, (int32_t *) &dataset -> data_type ) ) return(0);
	if ( !zdf_uint32_read( zdf, &dataset -> ndims ) ) return(0);

    for( unsigned i=0; i < dataset -> ndims; i++) {
	 	if ( !zdf_uint64_read( zdf, &dataset -> count[i] ) ) return(0);
    }
    return( size_zdf_dataset_header( dataset ) );
}

/**
 * Adds dataset to file
 * @param  zdf     ZDF file handle
 * @param  dataset Dataset object (includes pointer to data)
 * @return         Number of bytes written on success, 0 on error. The dataset object is also
 *                 modified to include a unique id.
 */
size_t zdf_add_dataset( t_zdf_file* zdf, t_zdf_dataset* dataset ){

	size_t length;

    length = zdf_sizeof( dataset -> data_type );
    for( unsigned i = 0; i < dataset -> ndims; i++ ) {
    	length *= dataset -> count[i];
    }

	length += size_zdf_dataset_header( dataset ) ;

	t_zdf_record rec = {
		.id_version = ZDF_DATASET_ID,
		.name       = (char *) dataset -> name,
		.length     = length
    };

  size_t reclen;
  if ( !( reclen = zdf_record_write( zdf, &rec ) ) ) return(0);

 	dataset -> id = ++ zdf -> ndatasets;
  if ( !zdf_dataset_header_write( zdf, dataset ) ) return(0);

 	size_t count = 1;
  for( unsigned i=0; i < dataset -> ndims; i++) count *= dataset -> count[i];

 	if ( !zdf_vector_write( zdf, dataset -> data, dataset -> data_type, count ) ) return(0);

 	return( reclen + length );
}

/**
 * Adds a chunked dataset header to file
 * @param  zdf     ZDF file handle
 * @param  dataset Dataset object
 * @return         Number of bytes written on success, 0 on error. The dataset object is also
 *                 modified to include a unique id.
 */
size_t zdf_start_cdset( t_zdf_file* zdf, t_zdf_dataset* dataset ) {

	t_zdf_record rec = {
		.id_version = ZDF_CDSET_START_ID,
		.name       = (char *) dataset -> name,
		.length     = size_zdf_dataset_header( dataset )
	};

	size_t ok, len;

	len = ( ok =  zdf_record_write( zdf, &rec ) );
	if ( !ok ) return(0);

 	dataset -> id = ++ zdf -> ndatasets;
	len += (ok = zdf_dataset_header_write( zdf, dataset ) );
	if ( !ok ) return(0);

	return( len );
}

/**
 * Returns size of chunk header
 * @param  dataset Dataset object
 * @return         Chunk header size in bytes
 */
size_t size_zdf_chunk_header(const t_zdf_dataset* dataset) {

	return size_zdf_uint32 + (size_zdf_uint32 + 16) + size_zdf_uint64 + 
	       size_zdf_uint32 + 3 * dataset -> ndims * size_zdf_uint64;
}

/**
 * Adds a chunk header to file
 * @param  zdf     ZDF file handle
 * @param  dataset Dataset object
 * @param  chunk   Chunk object
 * @return         Number of bytes written on success, 0 on error.
 */
size_t zdf_write_chunk_header( t_zdf_file* zdf, t_zdf_dataset* dataset, t_zdf_chunk* chunk ){

	char name[16];
	snprintf( name, 16, "%08llx-chunk", dataset -> id);

	size_t length;
	length = zdf_sizeof( dataset -> data_type );
	for( unsigned i = 0; i < dataset -> ndims; i++ ) {
	    	length *= chunk -> count[i];
	    }
	length += size_zdf_uint32 + 3 * dataset -> ndims * size_zdf_uint64;

	t_zdf_record rec = {
		.id_version = ZDF_CDSET_CHUNK_ID,
		.name       = (char *)name,
		.length     = length
	};

	size_t reclen;
	if ( !( reclen = zdf_record_write( zdf, &rec ) ) ) return(0);

	// ID of dataset this chunk belongs to
	if ( !zdf_uint32_write( zdf, dataset -> id ) ) return(0);

	// Size of chunk
	for( unsigned i=0; i < dataset -> ndims; i++) {
	 	if ( !zdf_uint64_write( zdf, chunk -> count[i] ) ) return(0);
	}

	// Start position
	for( unsigned i = 0; i < dataset -> ndims; i++) {
		if ( !zdf_uint64_write( zdf, chunk -> start[i] ) ) return(0);
	}

	// Stride
	for( unsigned i = 0; i < dataset -> ndims; i++) {
		if ( !zdf_uint64_write( zdf, chunk -> stride[i] ) ) return(0);
	}

	return( reclen + size_zdf_uint32 + 3 * dataset -> ndims * size_zdf_uint64 );

}


/**
 * Adds a chunk of dataset data to file
 * @param  zdf     ZDF file handle
 * @param  dataset Dataset object
 * @param  chunk   Chunk object (includes pointer to data). Data is expected to be contiguos in memory
 * @return         Number of bytes written on success, 0 on error.
 */
size_t zdf_write_cdset( t_zdf_file* zdf, t_zdf_dataset* dataset, t_zdf_chunk* chunk ){

	size_t headerlen = zdf_write_chunk_header( zdf, dataset, chunk );
    if ( !headerlen ) return(0);

 	// Chunk data
 	size_t count = 1;
 	for( unsigned i = 0; i < dataset -> ndims; i++ ) count *= chunk -> count[i];

 	size_t vectorlen;
 	if ( !( vectorlen = zdf_vector_write( zdf, chunk -> data, dataset -> data_type, count ) ) )
 			return(0);

 	return( headerlen + vectorlen );
}

/**
 * Add an end marker for a chunked dataset to file
 * @param  zdf     ZDF file handle
 * @param  dataset Dataset object
 * @return         Number of bytes written on success, 0 on error.
 */
size_t zdf_end_cdset( t_zdf_file* zdf, t_zdf_dataset* dataset ){

	char name[16];
	snprintf( name, 16, "%08llx-end", dataset -> id);

	t_zdf_record rec = {
		.id_version = ZDF_CDSET_END_ID,
		.name       = name,
		.length     = 0
    };

    size_t reclen;
    if ( !( reclen = zdf_record_write( zdf, &rec ) ) ) return(0);

 	return(reclen);
}


/**
 * Open a dataset and read the header. The dataset is selected by name.
 * @param  zdf     ZDF file handle
 * @param  dataset Dataset object
 * @return         Returns 1 if successful, 0 on error.
 */
size_t zdf_open_dataset( t_zdf_file* zdf, t_zdf_dataset* dataset ) {

	// Set file position after the magic number
	if ( fseek( zdf -> fp, 4, SEEK_SET ) ) return(-1);

	// Locate the requested dataset
	// If the dataset is not found there will be an EOF error signaled
	// and the routine returns -1
	int found = 0;
	while(!feof(zdf -> fp)) {
		t_zdf_record rec = { .name = NULL };

		if ( !zdf_record_read( zdf, &rec ) ) return(0);

		if ( rec.id_version == ZDF_CDSET_START_ID || rec.id_version == ZDF_DATASET_ID ) {
			int diff = strcmp( dataset -> name, rec.name );
			free( rec.name );

			if ( diff == 0 ) {
				// Read the dataset header
				if ( !zdf_dataset_header_read( zdf, dataset ) ) return(0);
				found = 1;
				break;
			} else {
				// Skip the record
				if ( !zdf_record_skip( zdf, &rec ) ) return(0);
			}
		} else {
			free( rec.name );
			if ( !zdf_record_skip( zdf, &rec ) ) return(0);
		}
	}

	if ( !found ) {
		fprintf(stderr,"(*error*) Unable to find dataset %s\n", dataset -> name);
		return(-1);
	}

	// re-position the file pointer at the end of the file
	if ( fseeko( zdf->fp, 0, SEEK_END ) ) return(0);
	return(1);

}

/**
 * Extends dataset dimensions. After this call the file pointer is positioned at the
 * end of the file.
 *
 * @param  zdf       ZDF file handle
 * @param  dataset   Datasete object
 * @param  new_count New dimensions for the dataset, must be >= than the previous
 *                   dimensions.
 * @return           Returns 1 if successful, 0 otherwise.
 */
int zdf_extend_dataset( t_zdf_file* zdf, t_zdf_dataset* dataset, uint64_t* new_count ) {

	for( unsigned i = 0; i < dataset -> ndims; i++ ) {
		if ( new_count[i] < dataset -> count[i] ) {
 	  		fprintf(stderr,"(*error*) Invalid value for zdf_extend_dataset.\n");
 	  		fprintf(stderr,"(*error*) New size is smaller than original size.\n");
 			return(-1);
		}
		dataset -> count[i] = new_count[i];
	}

    off_t offset = dataset -> offset;
	if ( fseeko( zdf -> fp, offset , SEEK_SET ) ) return(0);
    if ( ! zdf_dataset_header_write( zdf, dataset ) ) return(0);
	if ( fseeko( zdf->fp, 0, SEEK_END ) ) return(0);

	return(1);
}

/* -----------------------------------------------------------------------------------------------
  zdf high level interface
-------------------------------------------------------------------------------------------------- */

/**
 * Opens ZDF file and adds TYPE, GRID, and ITERATION metadata
 * @param  zdf       File handle
 * @param  info      Grid information
 * @param  iteration Iteration information
 * @param  path      File path
 * @return           Returns 1 on success, 0 on error
 */
int zdf_open_grid_file( t_zdf_file *zdf, const t_zdf_grid_info *info,
	const t_zdf_iteration *iteration, char const path[] ){

	char filename[1024];

	// Ensure that the path is available
	create_path( path );

	// Build filename
	sprintf( filename, "%s/%s-%06u.zdf", path, info->label, (unsigned) iteration -> n );
    // printf("Saving filename %s\n", filename );

    // Create ZDF file
    if ( !zdf_open_file( zdf, filename, ZDF_CREATE ) ) {
    	fprintf(stderr,"(*error*) Unable to open ZDF file, aborting.\n");
    	return(-1);
    }

    // Add file type
    if ( !zdf_add_string( zdf, "TYPE", "grid") ) return(0);

    // Add grid info
    if ( !zdf_add_grid_info( zdf, info ) ) return(0);

    // Add iteration info
    if ( !zdf_add_iteration( zdf, iteration ) ) return(0);

    return( 1 );
}


/**
 * Saves a ZDF grid file
 * @param  data       Pointer to grid data
 * @param  data_type  ZDF Data type of grid data
 * @param  info       Grid information
 * @param  iteration  Iteration information
 * @param  path       File path
 * @return            Returns 1 on success, 0 on error
 */
int zdf_save_grid( const void * data, const enum zdf_data_type data_type, const t_zdf_grid_info *info,
	const t_zdf_iteration *iteration, char const path[] )
{

    t_zdf_file zdf;

    // Open grid file
    if ( !zdf_open_grid_file( &zdf, info, iteration, path ) ) return(0);

    // Add dataset
    t_zdf_dataset dataset = {
    	.name = info -> name,
        .data_type = data_type,
    	.ndims = info->ndims,
    	.data = (void *) data
    };
    for( unsigned i = 0; i < info->ndims; i ++) dataset.count[i] = info->count[i];

    if ( !zdf_add_dataset( &zdf, &dataset ) ) return(0);

    // Close ZDF file and return
    return( zdf_close_file( &zdf ) );
}

/**
 * Opens ZDF file and adds TYPE, PARTICLES, and ITERATION metadata
 * @param  zdf        File handle
 * @param  info       Particles information
 * @param  iteration  Iteration information
 * @param  path       File path
 * @return            1 on sucess, 0 on error
 */
int zdf_open_part_file( t_zdf_file *zdf, t_zdf_part_info *info,
	const t_zdf_iteration *iteration, char const path[] ){

	char filename[1024];

	// Ensure that the path is available
	create_path( path );

	// Build filename
	sprintf( filename, "%s/%s-%s-%06u.zdf", path, "particles", info->name, (unsigned) iteration -> n );
    //printf("Saving filename %s\n", filename );

    // Create ZDF file
    if ( !zdf_open_file( zdf, filename, ZDF_CREATE ) ) {
    	fprintf(stderr,"(*error*) Unable to open ZDF file, aborting.\n");
    	return(-1);
    }

    // Add file type
    if ( !zdf_add_string( zdf, "TYPE", "particles") ) return(0);

    // Add particle info
    if ( !zdf_add_part_info( zdf, info ) ) return(0);

    // Add iteration info
    if ( !zdf_add_iteration( zdf, iteration ) ) return(0);

    return(1);
}

/**
 * Adds individual particle quantity to file
 * @param  zdf  File handle
 * @param  name Quantity name
 * @param  data Quantity data (float32)
 * @param  np   Number of particles
 * @return      Number of bytes written on success, 0 on error
 */
int zdf_add_quant_part_file( t_zdf_file *zdf, const char *name, const float* data,
	const uint64_t np ) {

    t_zdf_dataset dataset = {
    	.name = (char *) name,
        .data_type = zdf_float32,
    	.ndims = 1,
    	.data = (void *) data
    };

    dataset.count[0] = np;

    return( zdf_add_dataset( zdf, &dataset ) );

}

#ifdef __TEST_ZDF__

#include <math.h>

int main (int argc, const char * argv[]) {

	const unsigned NX = 128;
	float buf[NX];

	for(unsigned i = 0; i < NX; i++) {
		float x = 8 * (M_PI/NX) * (i+1);
		buf[i] = sin(x)/x;
	}

    t_zdf_grid_axis axis[1];
    axis[0] = (t_zdf_grid_axis) {
    	.min = -1.0,
    	.max =  1.0,
    	.label = "axis label",
    	.units = "axis units"
    };

    t_zdf_grid_info info = {
    	.ndims = 1,
    	.label = "data label",
    	.units = "data units",
    	.axis = axis
    };

    info.count[0] = NX;

    t_zdf_iteration iter = {
    	.n = 123,
    	.t = 12.3,
    	.time_units = "time units"
    };

	zdf_save_grid( (void *) buf, zdf_float32, &info, &iter, "V1TEST" );

	return 0;
}


#if 0
int main (int argc, const char * argv[]) {

	const unsigned NX = 128;
	float buf[NX];

	for(unsigned i = 0; i < NX; i++) {
		float x = 8 * (M_PI/NX) * (i+1);
		buf[i] = sin(x)/x;
	}

    t_zdf_grid_axis axis[1];
    axis[0] = (t_zdf_grid_axis) {
    	.min = -1.0,
    	.max =  1.0,
    	.label = "axis label",
    	.units = "axis units"
    };

    t_zdf_grid_info info = {
    	.ndims = 1,
    	.label = "data label",
    	.units = "data units",
    	.axis = axis
    };

    info.count[0] = NX;

    t_zdf_iteration iter = {
    	.n = 123,
    	.t = 12.3,
    	.time_units = "time units"
    };

	// Open file
	t_zdf_file file;
	zdf_open_grid_file( &file, &info, &iter, "chunk_test" );

	// Write chunked dataset header
	t_zdf_dataset dset = {
		.data_type = zdf_float32,
		.ndims = 1
	};
	dset.count[0] = 128;
	zdf_start_cdset( &file, "DATA", &dset );

	// Write chunks
	t_zdf_chunk chunk;
	chunk.count[0] = 16;
	chunk.start[0] = 0;
	chunk.stride[0] = 1;
	for( int i = 0; i < 8; i ++) {
		chunk.data = &buf[i*16];
		zdf_write_cdset( &file, &dset, &chunk );
		chunk.start[0] += 16;
	}

	// End dataset (optional)
    zdf_end_cdset( &file, &dset );

    // Close file
	zdf_close_file( &file );

	return 0;
}

#endif

#endif

