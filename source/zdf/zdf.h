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

#ifndef __ZDF__
#define __ZDF__

#if !(defined _WIN32 || defined _WIN64)

  #if __SIZEOF_FLOAT__ != 4
  #error This version of ZDF only support systems were sizeof(float) = 4
  #endif

  #if __SIZEOF_DOUBLE__ != 8
  #error This version of ZDF only support systems were sizeof(double) = 8
  #endif

#endif

#include <stdint.h>
#include <stdio.h>

/**
 * Number of bytes in each ZDF file unit
 */
#define BYTES_PER_ZDF_UNIT	(4)

/**
 * Magic byte sequence identifying ZDF files
 */
#define ZDF_MAGIC_LENGTH BYTES_PER_ZDF_UNIT
static const char zdf_magic[ZDF_MAGIC_LENGTH] = {'Z','D','F','1'};

/**
 * Maximum array dimensions supported
 */
#define zdf_max_dims 3

/**
 * ZDF datatypes - Currently only float32 and float64 are supported
 */
enum zdf_data_type{ zdf_null,
	                zdf_int8,  zdf_uint8,  zdf_int16, zdf_uint16,
	                zdf_int32, zdf_uint32, zdf_int64, zdf_uint64,
	                zdf_float32, zdf_float64 };


enum zdf_file_access_mode {ZDF_CREATE, ZDF_READ, ZDF_UPDATE };

typedef struct {
	FILE *fp;
	enum zdf_file_access_mode mode;
	uint32_t ndatasets;
} t_zdf_file;

typedef struct {
	char* name;
 	enum zdf_data_type data_type;
	uint32_t ndims;
	uint64_t count[zdf_max_dims];
	void * data;
	uint64_t id;
	uint64_t offset;
} t_zdf_dataset;

typedef struct {
	uint64_t count[zdf_max_dims];
	uint64_t start[zdf_max_dims];
	uint64_t stride[zdf_max_dims];
	void * data;
} t_zdf_chunk;


enum zdf_axis_type{ zdf_linear, zdf_log10, zdf_log2 };

typedef struct {
	char* name;
	enum zdf_axis_type type;

	double min;
	double max;

	char* label;
	char* units;

} t_zdf_grid_axis;

typedef struct {
	char* name;
	uint32_t ndims;
	uint64_t count[zdf_max_dims];
	char* label;
	char* units;
    t_zdf_grid_axis *axis;
} t_zdf_grid_info;

typedef struct {
	char* name;
	int32_t n;
	double t;
	char* time_units;
} t_zdf_iteration;


typedef struct {
	char* name;
	char* label;
	uint64_t np;
	uint32_t nquants;
	char** quants;
	char** qlabels;
	char** qunits;

} t_zdf_part_info;

typedef struct {
	char* name;
	char* label;
	uint32_t ntracks;
	uint32_t ndump;
	uint32_t niter;
	uint32_t nquants;

	char** quants;
	char** qlabels;
	char** qunits;
} t_zdf_track_info;

// Low level interface

size_t zdf_sizeof( enum zdf_data_type data_type );

int zdf_open_file( t_zdf_file* zdf, const char* filename, enum zdf_file_access_mode mode );

int zdf_close_file( t_zdf_file* zdf );

size_t zdf_vector_write( t_zdf_file* zdf, const void * data, enum zdf_data_type data_type, size_t len );

size_t zdf_add_string( t_zdf_file* zdf, const char* name, const char* str );

size_t zdf_add_int32( t_zdf_file* zdf, const char* name, const int32_t value );

size_t zdf_add_double( t_zdf_file* zdf, const char* name, const double value );

size_t zdf_add_iteration( t_zdf_file* zdf, const t_zdf_iteration* iter );

size_t zdf_add_grid_info( t_zdf_file* zdf, const t_zdf_grid_info* grid );

size_t zdf_add_part_info( t_zdf_file* zdf, const t_zdf_part_info* part );

size_t zdf_add_track_info( t_zdf_file* zdf, const t_zdf_track_info* tracks );

size_t zdf_add_dataset( t_zdf_file* zdf, t_zdf_dataset* dataset );

// Chunked dataset interface

size_t zdf_start_cdset( t_zdf_file* zdf, t_zdf_dataset* dataset );

size_t size_zdf_chunk_header(const t_zdf_dataset* dataset);

size_t zdf_write_chunk_header( t_zdf_file* zdf, t_zdf_dataset* dataset, t_zdf_chunk* chunk );

size_t zdf_write_cdset( t_zdf_file* zdf, t_zdf_dataset* dataset, t_zdf_chunk* chunk );

size_t zdf_end_cdset( t_zdf_file* zdf, t_zdf_dataset* dataset );

size_t zdf_open_dataset( t_zdf_file* zdf, t_zdf_dataset* dataset );

int zdf_extend_dataset( t_zdf_file* zdf, t_zdf_dataset* dataset, uint64_t* new_count );


// High level interface

int zdf_open_grid_file( t_zdf_file *file, const t_zdf_grid_info *info,
	const t_zdf_iteration *iteration, char const path[] );

int zdf_save_grid( const void* data, enum zdf_data_type data_type, const t_zdf_grid_info *info,
	const t_zdf_iteration *iteration, char const path[] );

int zdf_open_part_file( t_zdf_file *file, t_zdf_part_info *info,
	const t_zdf_iteration *iteration, char const path[] );

int zdf_add_quant_part_file( t_zdf_file *zdf, const char *name, const float* data,
	const uint64_t np );


#endif


